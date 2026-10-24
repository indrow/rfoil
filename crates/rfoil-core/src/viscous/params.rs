//! Calibration constants and free-stream parameters of the boundary layer.

use super::closure::Amplification;
use super::dual::{slot, D};

/// Turbulence-model calibration constants.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Calibration {
    /// Shear-lag constant.
    pub shear_lag: f64,
    /// G–β equilibrium locus `G = A·√(1 + B·β) + C/(H·Rθ·√(Cf/2))`.
    pub g_beta_a: f64,
    pub g_beta_b: f64,
    pub g_beta_c: f64,
    /// Wake/wall dissipation length ratio.
    pub wake_dissipation_length: f64,
    /// Initial shear stress at transition, `Cτ = scale·exp(−exponent/(Hk−1))·Cτ_eq`.
    pub transition_ctau_scale: f64,
    pub transition_ctau_exponent: f64,
    /// Weight of the equilibrium `dUe/dx` term in the lag equation.
    pub equilibrium_dudx_weight: f64,
    /// Turbulent skin-friction scale factor.
    pub cf_scale: f64,
}

impl Default for Calibration {
    fn default() -> Self {
        Self {
            shear_lag: 5.6,
            g_beta_a: 6.70,
            g_beta_b: 0.75,
            g_beta_c: 18.0,
            wake_dissipation_length: 0.9,
            transition_ctau_scale: 1.8,
            transition_ctau_exponent: 3.3,
            equilibrium_dudx_weight: 1.0,
            cf_scale: 1.0,
        }
    }
}

impl Calibration {
    /// Scale of the equilibrium shear coefficient implied by the G–β locus.
    pub fn ctau_eq_scale(&self) -> f64 {
        0.5 / (self.g_beta_a * self.g_beta_a * self.g_beta_b)
    }
}

/// `γ − 1` for air.
pub(crate) const GAMMA_MINUS_1: f64 = 0.4;
/// Sutherland's-law temperature ratio; zero means viscosity ∝ temperature^1.5.
pub(crate) const SUTHERLAND_RATIO: f64 = 0.0;

/// Free-stream parameters as dual numbers in the `MACH_SQ` and `REYNOLDS` slots.
#[derive(Clone, Copy, Debug)]
pub struct Freestream {
    pub speed: f64,
    /// Kármán–Tsien parameter λ.
    pub karman_tsien: D,
    pub stagnation_density: D,
    /// Reciprocal of the stagnation enthalpy.
    pub inv_stagnation_enthalpy: D,
    /// Reynolds number based on free-stream conditions.
    pub reynolds: D,
}

impl Freestream {
    pub fn new(mach: f64, reynolds: f64, speed: f64) -> Self {
        let msq = D::var(mach * mach, slot::MACH_SQ);
        let re = D::var(reynolds, slot::REYNOLDS);
        let beta = (1.0 - msq).sqrt();
        let stagnation_temperature_ratio = 1.0 + 0.5 * GAMMA_MINUS_1 * msq;
        let inv_stagnation_enthalpy = GAMMA_MINUS_1 * msq / (speed * speed) / stagnation_temperature_ratio;
        let enthalpy_ratio = 1.0 - 0.5 * speed * speed * inv_stagnation_enthalpy;
        Self {
            speed,
            karman_tsien: msq / ((1.0 + beta) * (1.0 + beta)),
            stagnation_density: stagnation_temperature_ratio.powf(1.0 / GAMMA_MINUS_1),
            inv_stagnation_enthalpy,
            reynolds: re * enthalpy_ratio.powf(1.5) * (1.0 + SUTHERLAND_RATIO) / (enthalpy_ratio + SUTHERLAND_RATIO),
        }
    }

    /// Edge Mach number squared at an incompressible-equivalent speed `ue` (value only).
    pub fn edge_msq(&self, ue: f64) -> f64 {
        let h = self.inv_stagnation_enthalpy.value;
        ue * ue * h / (GAMMA_MINUS_1 * (1.0 - 0.5 * ue * ue * h))
    }
}

/// Everything a boundary-layer interval needs besides its two stations.
#[derive(Clone, Copy, Debug)]
pub struct Context {
    pub freestream: Freestream,
    pub calibration: Calibration,
    /// Critical amplification exponent of the eᴺ method.
    pub ncrit: f64,
    pub amplification: Amplification,
    /// Arc length of the forced-transition trip on the current side.
    pub xi_trip: f64,
}
