use pyo3::prelude::*;
use rfoil_core::Airfoil;

use crate::to_py_err;

/// Airfoil coordinates ordered
/// from the TE over the upper surface to the LE, then back along the lower surface to the TE.
#[pyclass(name = "Airfoil", module = "rfoil._core", frozen)]
#[derive(Clone)]
pub struct PyAirfoil {
    pub inner: Airfoil,
}

#[pymethods]
impl PyAirfoil {
    /// Builds an airfoil from coordinates (clockwise input is reversed).
    #[new]
    fn new(name: String, x: Vec<f64>, y: Vec<f64>) -> PyResult<Self> {
        Ok(Self { inner: Airfoil::new(name, x, y).map_err(to_py_err)? })
    }

    /// NACA 4- or 5-digit section, e.g. ``Airfoil.naca("2412")``.
    #[staticmethod]
    fn naca(designation: &str) -> PyResult<Self> {
        Ok(Self { inner: Airfoil::naca(designation).map_err(to_py_err)? })
    }

    /// Parses coordinate-file text (plain, labelled, ISES or Lednicer format).
    #[staticmethod]
    fn parse(text: &str) -> PyResult<Self> {
        Ok(Self { inner: Airfoil::parse(text).map_err(to_py_err)? })
    }

    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn x(&self) -> Vec<f64> {
        self.inner.x.clone()
    }

    #[getter]
    fn y(&self) -> Vec<f64> {
        self.inner.y.clone()
    }

    /// A copy with a new name.
    fn renamed(&self, name: String) -> Self {
        Self { inner: Airfoil { name, ..self.inner.clone() } }
    }

    /// Coordinates as a name line followed by ``x y`` pairs.
    fn to_dat(&self) -> String {
        self.inner.to_dat()
    }

    /// ``(thickness, x_thickness, camber, x_camber)`` of the raw coordinates.
    fn shape(&self) -> (f64, f64, f64, f64) {
        let s = self.inner.shape();
        (s.thickness, s.x_thickness, s.camber, s.x_camber)
    }

    fn __len__(&self) -> usize {
        self.inner.x.len()
    }

    fn __repr__(&self) -> String {
        format!("Airfoil({:?}, {} points)", self.inner.name, self.inner.x.len())
    }
}
