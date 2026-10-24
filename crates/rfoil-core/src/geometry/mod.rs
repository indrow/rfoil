//! Airfoil geometry: coordinates, splined contours, NACA sections, paneling.

mod airfoil;
mod curve;
mod naca;
mod paneling;

pub use airfoil::Airfoil;
pub use curve::{ContourPoint, Curve, SectionShape};
pub use paneling::{panel, PanelingOptions};
