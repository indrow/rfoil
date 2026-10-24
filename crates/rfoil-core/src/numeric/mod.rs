//! Numerical building blocks: splines, dense linear algebra, dual numbers.

pub mod dual;
pub mod linalg;
pub mod spline;

pub use dual::Dual;
pub use linalg::{solve_small, solve_tridiagonal, Lu, Matrix};
pub use spline::Spline;
