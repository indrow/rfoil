//! Independent variables of one boundary-layer interval.
//!
//! Each interval couples an upstream station "1" and a downstream station
//! "2".  Every quantity is a [`D`] carrying partials with respect to the
//! primary variables of both stations plus a few global parameters, so the
//! Newton Jacobian falls out of ordinary arithmetic.

use crate::numeric::Dual;

/// Column of each primary variable within one station's Jacobian block.
pub mod var {
    /// Amplification `N` while laminar, shear-stress coefficient `√Cτ` once turbulent.
    pub const AMPL_OR_SHEAR: usize = 0;
    pub const THETA: usize = 1;
    pub const DSTAR: usize = 2;
    pub const UE: usize = 3;
    pub const XI: usize = 4;
    /// Number of primary variables per station.
    pub const COUNT: usize = 5;
}

/// Partial-derivative slots of [`D`].
pub mod slot {
    use super::var;

    /// Slot of primary variable `v` of the upstream station.
    pub const fn upstream(v: usize) -> usize {
        v
    }

    /// Slot of primary variable `v` of the downstream station.
    pub const fn downstream(v: usize) -> usize {
        var::COUNT + v
    }

    /// Free-stream Mach number squared.
    pub const MACH_SQ: usize = 10;
    pub const REYNOLDS: usize = 11;
    /// Arc length of the forced-transition trip.
    pub const XI_TRIP: usize = 12;
    /// Trial downstream amplification, used only while locating transition.
    pub const AMPL_TRIAL: usize = 13;
    pub const COUNT: usize = 14;
}

/// Dual number over the interval variables.
pub type D = Dual<{ slot::COUNT }>;

/// Moves downstream-station partials into the upstream-station slots.
pub fn as_upstream(mut x: D) -> D {
    for v in 0..var::COUNT {
        x.partials[slot::upstream(v)] = x.partials[slot::downstream(v)];
        x.partials[slot::downstream(v)] = 0.0;
    }
    x
}

/// Value-only clamp from below that keeps the partials, so the Newton
/// linearisation still sees the unclamped trend.
pub fn floor_value(mut x: D, min: f64) -> D {
    x.value = x.value.max(min);
    x
}

/// Value-only clamp from above that keeps the partials.
pub fn ceil_value(mut x: D, max: f64) -> D {
    x.value = x.value.min(max);
    x
}
