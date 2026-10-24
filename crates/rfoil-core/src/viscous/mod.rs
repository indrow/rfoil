//! Integral boundary layer with eᴺ transition, strongly coupled to the panel method.

pub mod closure;
pub mod dual;
pub mod interval;
pub mod layer;
pub mod march;
pub mod newton;
pub mod params;
pub mod station;
pub mod transition;

pub use closure::Amplification;
pub use layer::{BoundaryLayer, Node, Side, TransitionPoint, LOWER, UPPER};
pub use params::Calibration;
