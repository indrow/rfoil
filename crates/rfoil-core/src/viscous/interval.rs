//! Discretised boundary-layer equations over one interval.
//!
//! Every interval contributes three residuals:
//!
//! 1. amplification (laminar) or shear-lag (turbulent / wake) equation,
//! 2. von Kármán momentum equation,
//! 3. kinetic-energy shape-parameter equation,
//!
//! written in logarithmic differences between stations 1 and 2.

use super::closure;
use super::dual::{ceil_value, slot, var, D};
use super::params::Context;
use super::station::{low_reynolds_hk_excess, Regime, Station};
use super::transition::{mean_growth, GrowthPoint, Onset};

/// Pressure-gradient parameter `ξ/Ue·dUe/dξ` assumed at the similarity station.
const SIMILARITY_PRESSURE_GRADIENT: f64 = 1.0;

/// Residuals of one interval, as duals over the interval variables.
pub type Residuals = [D; 3];

/// What governs the interval ending at a station.
#[derive(Clone, Copy, Debug)]
pub enum Segment {
    /// First station after the stagnation point (self-similar flow).
    Similarity,
    Laminar,
    /// Laminar flow up to the onset, turbulent after it.
    Transition(Onset),
    Turbulent,
    Wake,
    /// First wake interval, starting from the combined trailing-edge values.
    TrailingEdge {
        shear: f64,
        theta: f64,
        dstar: f64,
    },
}

impl Segment {
    /// Chooses the interval type from the marching state.
    pub fn select(similarity: bool, wake: bool, turbulent: bool, onset: Option<Onset>) -> Self {
        match onset {
            Some(o) => Segment::Transition(o),
            None if similarity => Segment::Similarity,
            None if !turbulent => Segment::Laminar,
            None if wake => Segment::Wake,
            None => Segment::Turbulent,
        }
    }

    /// Closure regime of the downstream station.
    pub fn regime(&self) -> Regime {
        match self {
            Segment::Similarity | Segment::Laminar => Regime::Laminar,
            Segment::Transition(_) | Segment::Turbulent => Regime::Turbulent,
            Segment::Wake | Segment::TrailingEdge { .. } => Regime::Wake,
        }
    }

    /// Closes the downstream station and evaluates the interval residuals.
    pub fn equations(&self, s1: &Station, s2: &mut Station, ctx: &Context) -> Residuals {
        s2.close(self.regime(), &ctx.calibration);
        match *self {
            Segment::Similarity => similarity(s2, ctx),
            Segment::Laminar => residuals(Kind::Laminar, s1, s2, ctx),
            Segment::Turbulent => residuals(Kind::Turbulent, s1, s2, ctx),
            Segment::Wake => residuals(Kind::Wake, s1, s2, ctx),
            Segment::Transition(ref onset) => transition(onset, s1, s2, ctx),
            Segment::TrailingEdge { shear, theta, dstar } => trailing_edge(shear, theta, dstar, s2),
        }
    }
}

/// Type of an ordinary interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Laminar,
    Turbulent,
    Wake,
}

/// Logarithmic differences `ln(q2/q1)` across the interval.
struct LogRatios {
    xi: D,
    ue: D,
    theta: D,
    h_energy: D,
}

impl LogRatios {
    fn between(s1: &Station, s2: &Station) -> Self {
        Self {
            xi: (s2.xi / s1.xi).ln(),
            ue: (s2.ue / s1.ue).ln(),
            theta: (s2.theta / s1.theta).ln(),
            h_energy: (s2.h_energy / s1.h_energy).ln(),
        }
    }

    /// The self-similar flow at the stagnation point.
    fn similarity() -> Self {
        let m = SIMILARITY_PRESSURE_GRADIENT;
        Self { xi: 1.0.into(), ue: m.into(), theta: (0.5 * (1.0 - m)).into(), h_energy: 0.0.into() }
    }
}

/// Similarity station: the station is its own upstream end.
fn similarity(s: &Station, ctx: &Context) -> Residuals {
    let logs = LogRatios::similarity();
    let upwind = upwinding(s, s, false);
    [s.ampl, momentum(Regime::Laminar, s, s, &logs, ctx), shape(s, s, upwind, &logs)]
}

/// Residuals of an ordinary interval.
fn residuals(kind: Kind, s1: &Station, s2: &Station, ctx: &Context) -> Residuals {
    let logs = LogRatios::between(s1, s2);
    let upwind = upwinding(s1, s2, kind == Kind::Wake);
    let first = match kind {
        Kind::Laminar => {
            let rate = mean_growth(GrowthPoint::of(s1), GrowthPoint::of(s2), ctx);
            s2.ampl - s1.ampl - rate * (s2.xi - s1.xi)
        }
        Kind::Turbulent | Kind::Wake => shear_lag(kind == Kind::Wake, s1, s2, upwind, &logs, ctx),
    };
    let regime = match kind {
        Kind::Laminar => Regime::Laminar,
        Kind::Turbulent => Regime::Turbulent,
        Kind::Wake => Regime::Wake,
    };
    [first, momentum(regime, s1, s2, &logs, ctx), shape(s1, s2, upwind, &logs)]
}

fn mean(a: D, b: D) -> D {
    0.5 * (a + b)
}

/// Upwinding weight: 0.5 is trapezoidal, 1.0 backward Euler.  It kicks in
/// where `ln(Hk − 1)` changes rapidly, mainly at transition.
fn upwinding(s1: &Station, s2: &Station, wake: bool) -> D {
    let sharpness = if wake { 1.0 } else { 5.0 } / (s2.hk * s2.hk);
    let log_change = ((s2.hk - 1.0) / (s1.hk - 1.0)).ln();
    let log_change_sq = ceil_value(log_change * log_change, 15.0);
    1.0 - 0.5 * (-log_change_sq * sharpness).exp()
}

/// Lag equation for the maximum shear coefficient `√Cτ`.
fn shear_lag(wake: bool, s1: &Station, s2: &Station, upwind: D, logs: &LogRatios, ctx: &Context) -> D {
    let cal = &ctx.calibration;
    let shear = s1.shear.lerp(s2.shear, upwind);
    let shear_eq = s1.shear_eq.lerp(s2.shear_eq, upwind);
    let cf = s1.cf.lerp(s2.cf, upwind);
    let hk = s1.hk.lerp(s2.hk, upwind);
    let slip = mean(s1.slip, s2.slip);
    let re_theta = mean(s1.re_theta, s2.re_theta);
    let delta = mean(s1.delta, s2.delta);
    let dstar = mean(s1.dstar, s2.dstar);

    // the wake has a longer dissipation length
    let length_ratio = if wake { cal.wake_dissipation_length } else { 1.0 };

    // equilibrium 1/Ue·dUe/dξ
    let hk_excess = if wake { hk - 1.0 } else { low_reynolds_hk_excess(hk, re_theta, cal) };
    let hr = hk_excess / (cal.g_beta_a * length_ratio * hk);
    let equilibrium_dudx = (0.5 * cf - hr * hr) / (cal.g_beta_b * dstar);

    let lag_rate = cal.shear_lag * 1.333 / (1.0 + slip);
    let dxi = s2.xi - s1.xi;
    let log_shear = (s2.shear / s1.shear).ln();
    lag_rate * (shear_eq - shear * length_ratio) * dxi - delta * 2.0 * log_shear
        + delta * 2.0 * (equilibrium_dudx * dxi - logs.ue) * cal.equilibrium_dudx_weight
}

/// Integral momentum equation; the skin-friction term uses the midpoint Cf for accuracy.
fn momentum(regime: Regime, s1: &Station, s2: &Station, logs: &LogRatios, ctx: &Context) -> D {
    let cf_mid = midpoint_cf(regime, s1, s2, ctx.calibration.cf_scale);
    let xi = mean(s1.xi, s2.xi);
    let theta = mean(s1.theta, s2.theta);
    let friction = 0.5 * cf_mid * xi / theta + 0.25 * (s1.cf * s1.xi / s1.theta + s2.cf * s2.xi / s2.theta);
    let shape_factor = mean(s1.h, s2.h) + 2.0 - mean(s1.msq, s2.msq) + wake_gap_ratio(s1, s2);
    logs.theta + shape_factor * logs.ue - logs.xi * 0.5 * friction
}

/// Kinetic-energy shape-parameter equation.
fn shape(s1: &Station, s2: &Station, upwind: D, logs: &LogRatios) -> D {
    let xi_over_theta1 = s1.xi / s1.theta;
    let xi_over_theta2 = s2.xi / s2.theta;
    let dissipation = (s1.dissipation * xi_over_theta1).lerp(s2.dissipation * xi_over_theta2, upwind);
    let friction = (s1.cf * xi_over_theta1).lerp(s2.cf * xi_over_theta2, upwind);
    let shape_factor = 2.0 * mean(s1.h_density, s2.h_density) / mean(s1.h_energy, s2.h_energy) + 1.0
        - mean(s1.h, s2.h)
        - wake_gap_ratio(s1, s2);
    logs.h_energy + shape_factor * logs.ue + logs.xi * (0.5 * friction - dissipation)
}

/// Mean of `δ_wake/θ`, zero on the airfoil.
fn wake_gap_ratio(s1: &Station, s2: &Station) -> D {
    0.5 * (s1.wake_gap / s1.theta + s2.wake_gap / s2.theta)
}

/// Skin friction at the interval midpoint.
fn midpoint_cf(regime: Regime, s1: &Station, s2: &Station, cf_scale: f64) -> D {
    let (hk, re_theta, msq) = (mean(s1.hk, s2.hk), mean(s1.re_theta, s2.re_theta), mean(s1.msq, s2.msq));
    let fit = match regime {
        Regime::Wake => return D::constant(0.0),
        Regime::Laminar => closure::skin_friction_laminar(hk.value, re_theta.value),
        Regime::Turbulent => closure::skin_friction_turbulent(hk.value, re_theta.value, msq.value, cf_scale)
            .max(closure::skin_friction_laminar(hk.value, re_theta.value)),
    };
    fit.lift(hk, re_theta, msq)
}

/// Dummy interval between the trailing edge and the first wake station:
/// the wake starts with the combined TE values `(√Cτ, θ, δ*)`.
fn trailing_edge(shear: f64, theta: f64, dstar: f64, s2: &Station) -> Residuals {
    let shear1 = D::var(shear, slot::upstream(var::AMPL_OR_SHEAR));
    let theta1 = D::var(theta, slot::upstream(var::THETA));
    let dstar1 = D::var(dstar, slot::upstream(var::DSTAR));
    [s2.shear - shear1, s2.theta - theta1, s2.dstar + s2.wake_gap - dstar1]
}

/// Interval containing transition: the laminar part `ξ1..ξt` and the
/// turbulent part `ξt..ξ2` are discretised separately and summed.
fn transition(onset: &Onset, s1: &Station, s2: &Station, ctx: &Context) -> Residuals {
    let w = (onset.xi - s1.xi) / (s2.xi - s1.xi);
    let mut at_onset = Station {
        xi: onset.xi,
        ue: s1.ue.lerp(s2.ue, w),
        theta: s1.theta.lerp(s2.theta, w),
        dstar: s1.dstar.lerp(s2.dstar, w),
        ampl: D::constant(ctx.ncrit),
        wake_gap: s2.wake_gap,
        ..Station::default()
    }
    .with_kinematics(&ctx.freestream);

    at_onset.close(Regime::Laminar, &ctx.calibration);
    let laminar = residuals(Kind::Laminar, s1, &at_onset, ctx);

    // initial shear stress at transition from the equilibrium value
    let cal = &ctx.calibration;
    at_onset.close(Regime::Turbulent, cal);
    let ratio = cal.transition_ctau_scale * (-cal.transition_ctau_exponent / (at_onset.hk - 1.0)).exp();
    at_onset.shear = ratio * at_onset.shear_eq;
    at_onset.close(Regime::Turbulent, cal);
    let turbulent = residuals(Kind::Turbulent, &at_onset, s2, ctx);

    [turbulent[0], laminar[1] + turbulent[1], laminar[2] + turbulent[2]]
}

/// Newton linearisation `J₁·Δ₁ + J₂·Δ₂ = rhs` of the three interval residuals.
#[derive(Clone, Copy, Debug, Default)]
pub struct Linearization {
    /// Negated residuals.
    pub rhs: [f64; 3],
    /// Partials w.r.t. the upstream and downstream primaries, indexed by [`var`].
    pub upstream: [[f64; var::COUNT]; 3],
    pub downstream: [[f64; var::COUNT]; 3],
    /// Partials w.r.t. M∞², Re and the trip location.
    pub d_msq: [f64; 3],
    pub d_reynolds: [f64; 3],
    pub d_xi_trip: [f64; 3],
}

impl From<Residuals> for Linearization {
    fn from(residuals: Residuals) -> Self {
        let up = slot::upstream(0)..slot::upstream(var::COUNT);
        let down = slot::downstream(0)..slot::downstream(var::COUNT);
        let mut lin = Self::default();
        for (k, r) in residuals.iter().enumerate() {
            lin.rhs[k] = -r.value;
            lin.upstream[k].copy_from_slice(&r.partials[up.clone()]);
            lin.downstream[k].copy_from_slice(&r.partials[down.clone()]);
            lin.d_msq[k] = r.partials[slot::MACH_SQ];
            lin.d_reynolds[k] = r.partials[slot::REYNOLDS];
            lin.d_xi_trip[k] = r.partials[slot::XI_TRIP];
        }
        lin
    }
}
