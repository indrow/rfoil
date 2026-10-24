//! # rfoil-core
//!
//! Viscous/inviscid airfoil analysis: a linear-vorticity panel method coupled
//! to an integral boundary-layer formulation with eᴺ transition prediction,
//! solved simultaneously by a global Newton method.

// Numerical kernels index several arrays in lock-step; index loops read clearer there.
#![allow(clippy::needless_range_loop)]

pub mod analysis;
pub mod error;
pub mod geometry;
pub mod inviscid;
pub mod numeric;
pub mod viscous;

pub use analysis::{Analysis, OperatingPoint, Scaling, Settings};
pub use error::{Error, Result};
pub use geometry::{Airfoil, PanelingOptions};
