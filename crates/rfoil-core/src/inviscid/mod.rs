//! Linear-vorticity panel method with source coupling for the boundary layer.

pub mod forces;
pub mod influence;
pub mod panels;
pub mod system;
pub mod wake;

pub use forces::{Compressibility, Condition, Loads};
pub use influence::UniformFlow;
pub use panels::{Panels, TrailingEdge};
pub use system::InviscidSystem;
