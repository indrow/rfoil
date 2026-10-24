use std::panic::{self, AssertUnwindSafe};

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rfoil_core::analysis::IterationInfo;
use rfoil_core::{Analysis, OperatingPoint};

use crate::airfoil::PyAirfoil;
use crate::results::{PyBoundaryLayerSide, PyIteration, PyOperatingPoint, PySurface};
use crate::settings::{PyPaneling, PySettings};
use crate::to_py_err;

/// An airfoil analysis session: paneling, settings and the current solution.
///
/// The solver keeps its boundary layer between calls, so consecutive
/// operating points warm-start from each other.  Long solves
/// release the GIL, so they can run in a worker thread.
#[pyclass(name = "Analysis", module = "rfoil._core")]
pub struct PyAnalysis {
    inner: Analysis,
}

/// Calls a Python callback from a detached (GIL-free) solve and remembers the
/// first exception it raised.
struct Callback<'a> {
    func: Option<&'a Py<PyAny>>,
    error: Option<PyErr>,
}

impl<'a> Callback<'a> {
    fn new(func: Option<&'a Py<PyAny>>) -> Self {
        Self { func, error: None }
    }

    /// Returns `false` if the callback asked to stop (returned a falsy
    /// value other than `None`) or raised.
    fn call<T: for<'py> IntoPyObject<'py>>(&mut self, arg: T) -> bool {
        let Some(func) = self.func else { return true };
        if self.error.is_some() {
            return false;
        }
        Python::attach(|py| match func.call1(py, (arg,)) {
            Ok(result) => {
                let result = result.bind(py);
                result.is_none() || result.is_truthy().unwrap_or(true)
            }
            Err(e) => {
                self.error = Some(e);
                false
            }
        })
    }

    fn finish<T>(self, value: T) -> PyResult<T> {
        match self.error {
            Some(e) => Err(e),
            None => Ok(value),
        }
    }
}

#[pymethods]
impl PyAnalysis {
    #[new]
    #[pyo3(signature = (airfoil, paneling=None, settings=None))]
    fn new(
        airfoil: &PyAirfoil,
        paneling: Option<PyRef<'_, PyPaneling>>,
        settings: Option<PyRef<'_, PySettings>>,
    ) -> PyResult<Self> {
        let paneling = paneling.map(|p| p.to_core()).unwrap_or_default();
        let mut inner = Analysis::with_paneling(airfoil.inner.clone(), paneling).map_err(to_py_err)?;
        if let Some(s) = settings {
            inner.set_settings(s.to_core()?);
        }
        Ok(Self { inner })
    }

    #[getter]
    fn airfoil(&self) -> PyAirfoil {
        PyAirfoil { inner: self.inner.airfoil().clone() }
    }

    #[setter(airfoil)]
    fn set_airfoil(&mut self, airfoil: &PyAirfoil) -> PyResult<()> {
        self.inner.set_airfoil(airfoil.inner.clone()).map_err(to_py_err)
    }

    #[getter]
    fn settings(&self) -> PySettings {
        PySettings::from_core(self.inner.settings())
    }

    #[setter(settings)]
    fn set_settings(&mut self, settings: &PySettings) -> PyResult<()> {
        self.inner.set_settings(settings.to_core()?);
        Ok(())
    }

    #[getter]
    fn paneling(&self) -> PyPaneling {
        PyPaneling::from_core(self.inner.paneling())
    }

    #[setter(paneling)]
    fn set_paneling(&mut self, paneling: &PyPaneling) -> PyResult<()> {
        self.inner.set_paneling(paneling.to_core()).map_err(to_py_err)
    }

    /// Paneled airfoil node coordinates ``(x, y)``.
    #[getter]
    fn panels(&self) -> (Vec<f64>, Vec<f64>) {
        let p = self.inner.panels();
        (p.x[..p.n].to_vec(), p.y[..p.n].to_vec())
    }

    /// Wake node coordinates ``(x, y)`` (empty until a viscous solve).
    #[getter]
    fn wake(&self) -> (Vec<f64>, Vec<f64>) {
        let p = self.inner.panels();
        (p.x[p.n..].to_vec(), p.y[p.n..].to_vec())
    }

    /// ``(thickness, x_thickness, camber, x_camber)`` of the paneled airfoil.
    fn shape(&self) -> (f64, f64, f64, f64) {
        let s = self.inner.shape();
        (s.thickness, s.x_thickness, s.camber, s.x_camber)
    }

    /// Solves at angle of attack ``alpha`` (degrees).  ``progress`` is called
    /// with an :class:`Iteration` after each Newton step.
    #[pyo3(signature = (alpha, progress=None))]
    fn solve_alpha(&mut self, py: Python<'_>, alpha: f64, progress: Option<Py<PyAny>>) -> PyResult<PyOperatingPoint> {
        let mut cb = Callback::new(progress.as_ref());
        let point = self.guarded(py, |a| {
            a.solve_alpha_with(alpha, &mut |it: &IterationInfo| {
                cb.call(PyIteration::from(it));
            })
        })?;
        cb.finish(point.into())
    }

    /// Solves for lift coefficient ``cl``.
    #[pyo3(signature = (cl, progress=None))]
    fn solve_cl(&mut self, py: Python<'_>, cl: f64, progress: Option<Py<PyAny>>) -> PyResult<PyOperatingPoint> {
        let mut cb = Callback::new(progress.as_ref());
        let point = self.guarded(py, |a| {
            a.solve_cl_with(cl, &mut |it: &IterationInfo| {
                cb.call(PyIteration::from(it));
            })
        })?;
        cb.finish(point.into())
    }

    /// Sweeps α from ``start`` to ``end`` inclusive.  ``callback(point)`` is
    /// called for every point; returning ``False`` stops the sweep.
    #[pyo3(signature = (start, end, step, callback=None))]
    fn sweep_alpha(
        &mut self,
        py: Python<'_>,
        start: f64,
        end: f64,
        step: f64,
        callback: Option<Py<PyAny>>,
    ) -> PyResult<Vec<PyOperatingPoint>> {
        let mut cb = Callback::new(callback.as_ref());
        let points = self.guarded(py, |a| a.sweep_alpha(start, end, step, |p| cb.call(PyOperatingPoint::from(*p))))?;
        cb.finish(convert(points))
    }

    /// Sweeps CL from ``start`` to ``end`` inclusive (see :meth:`sweep_alpha`).
    #[pyo3(signature = (start, end, step, callback=None))]
    fn sweep_cl(
        &mut self,
        py: Python<'_>,
        start: f64,
        end: f64,
        step: f64,
        callback: Option<Py<PyAny>>,
    ) -> PyResult<Vec<PyOperatingPoint>> {
        let mut cb = Callback::new(callback.as_ref());
        let points = self.guarded(py, |a| a.sweep_cl(start, end, step, |p| cb.call(PyOperatingPoint::from(*p))))?;
        cb.finish(convert(points))
    }

    /// Forgets the boundary layer so the next point is initialised afresh.
    fn reset_boundary_layer(&mut self) {
        self.inner.reset_boundary_layer();
    }

    /// Pressure distribution of the current solution.
    fn surface(&self) -> PySurface {
        self.inner.surface().into()
    }

    /// ``(top, bottom)`` boundary-layer sides, or ``None`` if not viscous.
    fn boundary_layer(&self) -> Option<(PyBoundaryLayerSide, PyBoundaryLayerSide)> {
        self.inner.boundary_layer().map(|bl| (bl.top.into(), bl.bottom.into()))
    }

    fn __repr__(&self) -> String {
        let s = self.inner.settings();
        format!("Analysis({:?}, reynolds={}, mach={})", self.inner.airfoil().name, s.reynolds, s.mach)
    }
}

impl PyAnalysis {
    /// Runs a solve without the GIL.  A solver panic becomes a `RuntimeError`
    /// and discards the (possibly corrupted) boundary layer.
    fn guarded<T: Send>(&mut self, py: Python<'_>, f: impl FnOnce(&mut Analysis) -> T + Send) -> PyResult<T> {
        let inner = &mut self.inner;
        let result = py.detach(|| panic::catch_unwind(AssertUnwindSafe(|| f(inner))));
        result.map_err(|payload| {
            self.inner.reset_boundary_layer();
            let what = payload
                .downcast_ref::<String>()
                .map(String::as_str)
                .or_else(|| payload.downcast_ref::<&str>().copied())
                .unwrap_or("unknown error");
            PyRuntimeError::new_err(format!("solver failure: {what}"))
        })
    }
}

fn convert(points: Vec<OperatingPoint>) -> Vec<PyOperatingPoint> {
    points.into_iter().map(Into::into).collect()
}
