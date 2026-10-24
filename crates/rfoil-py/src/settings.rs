use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rfoil_core::viscous::Amplification;
use rfoil_core::{PanelingOptions, Scaling, Settings};

/// Operating conditions and solver options.
///
/// ``scaling`` is one of ``"fixed"`` (type 1), ``"lift"`` (type 2: fixed
/// M·√CL and Re·√CL) or ``"lift_speed"`` (type 3: fixed M and Re·CL).
/// ``amplification`` is ``"envelope"`` or ``"modified"``.
#[pyclass(name = "Settings", module = "rfoil._core", get_all, set_all)]
#[derive(Clone)]
pub struct PySettings {
    /// Reynolds number; 0 for inviscid analysis.
    pub reynolds: f64,
    pub mach: f64,
    pub ncrit: f64,
    pub xtr_top: f64,
    pub xtr_bottom: f64,
    pub max_iterations: usize,
    pub vaccel: f64,
    pub wake_length: f64,
    pub scaling: String,
    pub amplification: String,
}

#[pymethods]
impl PySettings {
    #[new]
    #[pyo3(signature = (reynolds=0.0, mach=0.0, ncrit=9.0, xtr_top=1.0, xtr_bottom=1.0, max_iterations=None, vaccel=None, wake_length=None, scaling="fixed".to_string(), amplification="envelope".to_string()))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        reynolds: f64,
        mach: f64,
        ncrit: f64,
        xtr_top: f64,
        xtr_bottom: f64,
        max_iterations: Option<usize>,
        vaccel: Option<f64>,
        wake_length: Option<f64>,
        scaling: String,
        amplification: String,
    ) -> Self {
        let d = Settings::default();
        Self {
            reynolds,
            mach,
            ncrit,
            xtr_top,
            xtr_bottom,
            max_iterations: max_iterations.unwrap_or(d.max_iterations),
            vaccel: vaccel.unwrap_or(d.vaccel),
            wake_length: wake_length.unwrap_or(d.wake_length),
            scaling,
            amplification,
        }
    }

    #[getter]
    fn viscous(&self) -> bool {
        self.reynolds > 0.0
    }

    fn __repr__(&self) -> String {
        format!(
            "Settings(reynolds={}, mach={}, ncrit={}, xtr_top={}, xtr_bottom={}, scaling={:?})",
            self.reynolds, self.mach, self.ncrit, self.xtr_top, self.xtr_bottom, self.scaling
        )
    }
}

impl PySettings {
    pub fn to_core(&self) -> PyResult<Settings> {
        let scaling = match self.scaling.as_str() {
            "fixed" => Scaling::Fixed,
            "lift" => Scaling::FixedLift,
            "lift_speed" => Scaling::FixedLiftAndSpeed,
            other => return Err(PyValueError::new_err(format!("unknown scaling {other:?}"))),
        };
        let amplification = match self.amplification.as_str() {
            "envelope" => Amplification::Envelope,
            "modified" => Amplification::Modified,
            other => return Err(PyValueError::new_err(format!("unknown amplification model {other:?}"))),
        };
        if !(0.0..1.0).contains(&self.mach) {
            return Err(PyValueError::new_err("Mach number must be in [0, 1)"));
        }
        Ok(Settings {
            reynolds: self.reynolds,
            mach: self.mach,
            scaling,
            ncrit: self.ncrit,
            xtr_top: self.xtr_top,
            xtr_bottom: self.xtr_bottom,
            max_iterations: self.max_iterations,
            vaccel: self.vaccel,
            wake_length: self.wake_length,
            amplification,
            ..Settings::default()
        })
    }

    pub fn from_core(s: &Settings) -> Self {
        Self {
            reynolds: s.reynolds,
            mach: s.mach,
            ncrit: s.ncrit,
            xtr_top: s.xtr_top,
            xtr_bottom: s.xtr_bottom,
            max_iterations: s.max_iterations,
            vaccel: s.vaccel,
            wake_length: s.wake_length,
            scaling: match s.scaling {
                Scaling::Fixed => "fixed",
                Scaling::FixedLift => "lift",
                Scaling::FixedLiftAndSpeed => "lift_speed",
            }
            .into(),
            amplification: match s.amplification {
                Amplification::Envelope => "envelope",
                Amplification::Modified => "modified",
            }
            .into(),
        }
    }
}

/// Panel distribution options.
#[pyclass(name = "Paneling", module = "rfoil._core", get_all, set_all)]
#[derive(Clone)]
pub struct PyPaneling {
    /// Number of panel nodes.
    pub nodes: usize,
    /// Curvature attraction (0 = uniform, ~1 typical).
    pub bunching: f64,
    /// TE / LE panel density ratio.
    pub te_le_ratio: f64,
    /// Refined-area / LE panel density ratio.
    pub refine_ratio: f64,
    pub top_refine: (f64, f64),
    pub bottom_refine: (f64, f64),
}

#[pymethods]
impl PyPaneling {
    #[new]
    #[pyo3(signature = (nodes=160, bunching=1.0, te_le_ratio=0.15, refine_ratio=0.2, top_refine=(1.0, 1.0), bottom_refine=(1.0, 1.0)))]
    fn new(
        nodes: usize,
        bunching: f64,
        te_le_ratio: f64,
        refine_ratio: f64,
        top_refine: (f64, f64),
        bottom_refine: (f64, f64),
    ) -> Self {
        Self { nodes, bunching, te_le_ratio, refine_ratio, top_refine, bottom_refine }
    }

    fn __repr__(&self) -> String {
        format!("Paneling(nodes={}, bunching={}, te_le_ratio={})", self.nodes, self.bunching, self.te_le_ratio)
    }
}

impl PyPaneling {
    pub fn to_core(&self) -> PanelingOptions {
        PanelingOptions {
            nodes: self.nodes,
            bunching: self.bunching,
            te_le_ratio: self.te_le_ratio,
            refine_ratio: self.refine_ratio,
            top_refine: self.top_refine,
            bottom_refine: self.bottom_refine,
        }
    }

    pub fn from_core(p: &PanelingOptions) -> Self {
        Self {
            nodes: p.nodes,
            bunching: p.bunching,
            te_le_ratio: p.te_le_ratio,
            refine_ratio: p.refine_ratio,
            top_refine: p.top_refine,
            bottom_refine: p.bottom_refine,
        }
    }
}
