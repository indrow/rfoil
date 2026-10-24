//! Python bindings for `rfoil-core`, exposed as the `rfoil._core` module.

mod airfoil;
mod analysis;
mod results;
mod settings;

use pyo3::prelude::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<airfoil::PyAirfoil>()?;
    m.add_class::<settings::PySettings>()?;
    m.add_class::<settings::PyPaneling>()?;
    m.add_class::<analysis::PyAnalysis>()?;
    m.add_class::<results::PyOperatingPoint>()?;
    m.add_class::<results::PyIteration>()?;
    m.add_class::<results::PySurface>()?;
    m.add_class::<results::PyBoundaryLayerSide>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}

/// Converts library errors into Python `ValueError`s.
pub(crate) fn to_py_err(e: rfoil_core::Error) -> PyErr {
    pyo3::exceptions::PyValueError::new_err(e.to_string())
}
