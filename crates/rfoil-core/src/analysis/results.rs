//! Analysis results.

/// Integrated results of one operating point.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct OperatingPoint {
    /// Angle of attack in degrees.
    pub alpha: f64,
    pub cl: f64,
    pub cm: f64,
    /// Total drag (Squire–Young at the end of the wake); zero if inviscid.
    pub cd: f64,
    /// Skin-friction drag.
    pub cdf: f64,
    /// Pressure drag `cd − cdf`.
    pub cdp: f64,
    /// Minimum pressure coefficient on the surface.
    pub cp_min: f64,
    /// Transition x/c on the top and bottom surfaces.
    pub xtr_top: f64,
    pub xtr_bottom: f64,
    pub mach: f64,
    pub reynolds: f64,
    pub viscous: bool,
    pub converged: bool,
    /// Newton iterations used.
    pub iterations: usize,
    /// Final RMS of the normalised Newton changes.
    pub rms_change: f64,
}

/// Pressure distribution over the airfoil (and wake, if viscous).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SurfaceDistribution {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Inviscid Cp.
    pub cp_inviscid: Vec<f64>,
    /// Viscous Cp (equal to inviscid Cp when inviscid).
    pub cp_viscous: Vec<f64>,
    /// Number of airfoil nodes; any further entries are wake nodes.
    pub airfoil_nodes: usize,
}

/// Boundary-layer quantities along one side.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BoundaryLayerSide {
    /// Node coordinates.
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Arc length from the stagnation point.
    pub xi: Vec<f64>,
    /// Edge speed `Ue/V∞`.
    pub ue: Vec<f64>,
    pub dstar: Vec<f64>,
    pub theta: Vec<f64>,
    /// Kinematic shape parameter `H = δ*/θ`.
    pub h: Vec<f64>,
    /// Skin friction coefficient `τw/(½ρV∞²)`.
    pub cf: Vec<f64>,
    /// Amplification ratio `N` (laminar stations, else zero).
    pub n: Vec<f64>,
    /// Shear stress `√Cτ` (turbulent stations, else zero).
    pub ctau: Vec<f64>,
    /// Index of the trailing-edge station in these arrays.
    pub te: usize,
}

/// Boundary-layer distributions of the last viscous solution.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BoundaryLayerResult {
    pub top: BoundaryLayerSide,
    /// Bottom surface followed by the wake.
    pub bottom: BoundaryLayerSide,
}

/// Progress of a Newton iteration, for live reporting.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IterationInfo {
    pub iteration: usize,
    pub rms: f64,
    pub max: f64,
    pub relaxation: f64,
    pub alpha: f64,
    pub cl: f64,
    pub cd: f64,
}
