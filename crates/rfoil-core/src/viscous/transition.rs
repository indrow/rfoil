//! eᴺ transition prediction: amplification growth and transition location.

use super::closure;
use super::dual::{ceil_value, slot, D};
use super::params::Context;
use super::station::Station;

/// Inputs of the amplification rate at one point.
#[derive(Clone, Copy, Debug)]
pub struct GrowthPoint {
    pub hk: D,
    pub theta: D,
    pub re_theta: D,
    pub ampl: D,
}

impl GrowthPoint {
    pub fn of(st: &Station) -> Self {
        Self { hk: st.hk, theta: st.theta, re_theta: st.re_theta, ampl: st.ampl }
    }
}

/// Mean amplification rate `dN/dξ` over an interval (RMS of the end values),
/// plus a small term keeping `dN/dξ > 0` near `N = Ncrit`.
pub fn mean_growth(p1: GrowthPoint, p2: GrowthPoint, ctx: &Context) -> D {
    let rate = |p: GrowthPoint| {
        let r = closure::amplification_rate(ctx.amplification, p.hk.value, p.theta.value, p.re_theta.value);
        D::chain(r.value, &[(r.d_hk, p.hk), (r.d_theta, p.theta), (r.d_re_theta, p.re_theta)])
    };
    let (rate1, rate2) = (rate(p1), rate(p2));
    let mean_sq = 0.5 * (rate1 * rate1 + rate2 * rate2);
    let rms = if mean_sq.value <= 0.0 { D::constant(0.0) } else { mean_sq.sqrt() };

    let arg = ceil_value(20.0 * (ctx.ncrit - 0.5 * (p1.ampl + p2.ampl)), 20.0);
    let near_critical = if arg.value <= 0.0 { D::constant(1.0) } else { (-arg).exp() };
    rms + near_critical * 0.002 / (p1.theta + p2.theta)
}

/// Where transition happens inside an interval.
#[derive(Clone, Copy, Debug)]
pub struct Onset {
    /// Transition arc length, with partials w.r.t. the interval variables.
    pub xi: D,
    /// True if set by the trip rather than by amplification.
    pub forced: bool,
}

/// Result of integrating the amplification equation across an interval.
#[derive(Clone, Copy, Debug)]
pub struct TransitionCheck {
    /// Amplification at the downstream station.
    pub downstream_ampl: f64,
    pub onset: Option<Onset>,
}

/// Solves the implicit amplification equation from `s1` to `s2`
///
/// ```text
///   (N2 − N1) / (ξ2 − ξ1) = ½·[N'(ξt, Nt) + N'(ξ1, N1)]
/// ```
///
/// where `(ξt, Nt)` is station 2 if `N2 < Ncrit`, or the interpolated point
/// where `N = Ncrit` otherwise, and reports any free or forced transition.
pub fn check(s1: &Station, s2: &Station, ctx: &Context) -> TransitionCheck {
    const TOLERANCE: f64 = 5.0e-5;
    let (xi1, xi2) = (s1.xi.value, s2.xi.value);
    let dxi = xi2 - xi1;

    // Newton iteration for N2, starting from an explicit step
    let rate = mean_growth(GrowthPoint::of(s1), GrowthPoint::of(s2), ctx);
    let mut ampl2 = s1.ampl.value + rate.value * dxi;
    for _ in 0..30 {
        let trial = Trial::new(s1, s2, ampl2, ctx);
        if trial.rate <= 0.0 {
            break;
        }
        let step = -trial.residual.value / trial.residual.partials[slot::AMPL_TRIAL];
        let xi_shift = trial.xi.partials[slot::AMPL_TRIAL] * step;
        // limit the step so that neither N nor the transition point jump far
        let mut relaxation = 1.0_f64;
        if relaxation * (xi_shift / dxi).abs() > 0.05 {
            relaxation = 0.05 * (dxi / xi_shift).abs();
        }
        if relaxation * step.abs() > 1.0 {
            relaxation = 1.0 / step.abs();
        }
        if step.abs() < TOLERANCE {
            break;
        }
        // land exactly on Ncrit when crossing it
        let next = ampl2 + relaxation * step;
        let crosses = (ampl2 > ctx.ncrit) != (next > ctx.ncrit) && ampl2 != ctx.ncrit;
        ampl2 = if crosses { ctx.ncrit } else { next };
    }

    let mut free = ampl2 >= ctx.ncrit;
    let mut forced = ctx.xi_trip > xi1 && ctx.xi_trip <= xi2;
    if !free && !forced {
        return TransitionCheck { downstream_ampl: ampl2, onset: None };
    }
    let trial = Trial::new(s1, s2, ampl2, ctx);
    if free && forced {
        // whichever comes first wins
        forced = ctx.xi_trip < trial.xi.value;
        free = !forced;
    }

    let xi = if forced {
        D::var(ctx.xi_trip, slot::XI_TRIP)
    } else {
        // eliminate N2 using the stationarity of the amplification residual
        let ratio = trial.xi.partials[slot::AMPL_TRIAL] / trial.residual.partials[slot::AMPL_TRIAL];
        let mut xi = trial.xi;
        for (x, r) in xi.partials.iter_mut().zip(trial.residual.partials) {
            *x -= ratio * r;
        }
        xi.partials[slot::AMPL_TRIAL] = 0.0;
        xi
    };
    debug_assert!(forced || free);
    TransitionCheck { downstream_ampl: ampl2, onset: Some(Onset { xi, forced }) }
}

/// The amplification equation evaluated for a trial downstream amplification.
struct Trial {
    /// End point of the laminar part of the interval.
    xi: D,
    residual: D,
    rate: f64,
}

impl Trial {
    fn new(s1: &Station, s2: &Station, ampl2: f64, ctx: &Context) -> Self {
        let a1 = s1.ampl;
        let a2 = D::var(ampl2, slot::AMPL_TRIAL);
        let (xi1, xi2) = (s1.xi, s2.xi);

        // fraction of the interval before free transition...
        let (free_fraction, ampl_end) = if a2.value <= ctx.ncrit {
            (D::constant(1.0), a2)
        } else {
            ((ctx.ncrit - a1) / (a2 - a1), D::constant(ctx.ncrit))
        };
        // ...or before the trip
        let trip_fraction = if ctx.xi_trip < xi2.value {
            (D::var(ctx.xi_trip, slot::XI_TRIP) - xi1) / (xi2 - xi1)
        } else {
            D::constant(1.0)
        };
        let w = if free_fraction.value < trip_fraction.value { free_fraction } else { trip_fraction };

        let xi = xi1.lerp(xi2, w);
        let end = Station {
            xi,
            ue: s1.ue.lerp(s2.ue, w),
            theta: s1.theta.lerp(s2.theta, w),
            dstar: s1.dstar.lerp(s2.dstar, w),
            ampl: ampl_end,
            ..Station::default()
        }
        .with_kinematics(&ctx.freestream);
        let rate = mean_growth(GrowthPoint::of(s1), GrowthPoint::of(&end), ctx);
        Self { xi, residual: a2 - a1 - rate * (xi2 - xi1), rate: rate.value }
    }
}
