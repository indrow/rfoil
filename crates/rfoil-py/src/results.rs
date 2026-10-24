use pyo3::prelude::*;
use pyo3::types::PyDict;
use rfoil_core::analysis::{BoundaryLayerSide, IterationInfo, SurfaceDistribution};
use rfoil_core::OperatingPoint;

/// Integrated results of one operating point.
#[pyclass(name = "OperatingPoint", module = "rfoil._core", frozen, get_all)]
#[derive(Clone)]
pub struct PyOperatingPoint {
    pub alpha: f64,
    pub cl: f64,
    pub cm: f64,
    pub cd: f64,
    pub cdf: f64,
    pub cdp: f64,
    pub cp_min: f64,
    pub xtr_top: f64,
    pub xtr_bottom: f64,
    pub mach: f64,
    pub reynolds: f64,
    pub viscous: bool,
    pub converged: bool,
    pub iterations: usize,
    pub rms_change: f64,
}

impl From<OperatingPoint> for PyOperatingPoint {
    fn from(p: OperatingPoint) -> Self {
        Self {
            alpha: p.alpha,
            cl: p.cl,
            cm: p.cm,
            cd: p.cd,
            cdf: p.cdf,
            cdp: p.cdp,
            cp_min: p.cp_min,
            xtr_top: p.xtr_top,
            xtr_bottom: p.xtr_bottom,
            mach: p.mach,
            reynolds: p.reynolds,
            viscous: p.viscous,
            converged: p.converged,
            iterations: p.iterations,
            rms_change: p.rms_change,
        }
    }
}

#[pymethods]
impl PyOperatingPoint {
    /// Lift-to-drag ratio (``inf`` when inviscid).
    #[getter]
    fn ld(&self) -> f64 {
        if self.cd > 0.0 {
            self.cl / self.cd
        } else {
            f64::INFINITY
        }
    }

    fn as_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        for (k, v) in [
            ("alpha", self.alpha),
            ("cl", self.cl),
            ("cd", self.cd),
            ("cdp", self.cdp),
            ("cdf", self.cdf),
            ("cm", self.cm),
            ("cp_min", self.cp_min),
            ("xtr_top", self.xtr_top),
            ("xtr_bottom", self.xtr_bottom),
        ] {
            d.set_item(k, v)?;
        }
        d.set_item("converged", self.converged)?;
        d.set_item("iterations", self.iterations)?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!(
            "OperatingPoint(alpha={:.3}, cl={:.4}, cd={:.5}, cm={:.4}, converged={})",
            self.alpha,
            self.cl,
            self.cd,
            self.cm,
            if self.converged { "True" } else { "False" }
        )
    }
}

/// Progress report of one Newton iteration.
#[pyclass(name = "Iteration", module = "rfoil._core", frozen, get_all)]
#[derive(Clone)]
pub struct PyIteration {
    pub iteration: usize,
    pub rms: f64,
    pub max: f64,
    pub relaxation: f64,
    pub alpha: f64,
    pub cl: f64,
    pub cd: f64,
}

impl From<&IterationInfo> for PyIteration {
    fn from(i: &IterationInfo) -> Self {
        Self {
            iteration: i.iteration,
            rms: i.rms,
            max: i.max,
            relaxation: i.relaxation,
            alpha: i.alpha,
            cl: i.cl,
            cd: i.cd,
        }
    }
}

/// Pressure distribution; entries past ``airfoil_nodes`` belong to the wake.
#[pyclass(name = "Surface", module = "rfoil._core", frozen, get_all)]
pub struct PySurface {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub cp_inviscid: Vec<f64>,
    pub cp_viscous: Vec<f64>,
    pub airfoil_nodes: usize,
}

impl From<SurfaceDistribution> for PySurface {
    fn from(s: SurfaceDistribution) -> Self {
        Self { x: s.x, y: s.y, cp_inviscid: s.cp_inviscid, cp_viscous: s.cp_viscous, airfoil_nodes: s.airfoil_nodes }
    }
}

/// Boundary-layer variables along one side; stations past ``te`` are the wake.
#[pyclass(name = "BoundaryLayerSide", module = "rfoil._core", frozen, get_all)]
pub struct PyBoundaryLayerSide {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub xi: Vec<f64>,
    pub ue: Vec<f64>,
    pub dstar: Vec<f64>,
    pub theta: Vec<f64>,
    pub h: Vec<f64>,
    pub cf: Vec<f64>,
    pub n: Vec<f64>,
    pub ctau: Vec<f64>,
    pub te: usize,
}

impl From<BoundaryLayerSide> for PyBoundaryLayerSide {
    fn from(s: BoundaryLayerSide) -> Self {
        Self {
            x: s.x,
            y: s.y,
            xi: s.xi,
            ue: s.ue,
            dstar: s.dstar,
            theta: s.theta,
            h: s.h,
            cf: s.cf,
            n: s.n,
            ctau: s.ctau,
            te: s.te,
        }
    }
}
