//! The global viscous/inviscid Newton system.
//!
//! Unknowns per station are `(ΔS, Δθ, Δm)`, where `S` is the amplification or
//! shear variable and `m = Ue·δ*` the mass defect.  Edge speeds follow from
//! the masses through the source-influence matrix, `Ue = Ue_inv + D·m`,
//! which couples every station to every other one:
//!
//! ```text
//!   A |  |  .  |  |  .  |    d       R       S
//!   B  A |  .  |  |  .  |    d       R       S
//!   |  B  A .  |  |  .  |    d       R       S
//!   .  .  .  .  |  |  .  |   d   =   R − Δp·S
//!   |  Z  |  |  B  A .  |    d       R       S
//!   |  |  |  |  |  |  B A    d       R       S
//! ```
//!
//! `A` (diagonal) and `B` (lower) are 3×2 blocks from the local interval
//! equations in `(S, θ)`; the dense `|` columns hold the mass coupling; `Z`
//! couples the first wake station to the upper-surface trailing edge.  The
//! second right-hand side `S` is the sensitivity to the floating global
//! parameter `p` (CL when α is prescribed, α when CL is prescribed).

mod assemble;
mod solve;
mod update;

pub use assemble::{assemble, Assembly};
pub use update::{apply, UpdateInputs, UpdateReport};

/// Three equations by two columns.
type Block = [[f64; 2]; 3];

/// Right-hand-side column of the residuals.
const RESIDUAL: usize = 0;
/// Right-hand-side column of the sensitivity to the floating parameter.
const PARAMETER: usize = 1;

/// Which global parameter floats while the other is held.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Constraint {
    /// α prescribed; CL floats.  Carries `dRe/dCL` and `dM²/dCL`.
    Alpha { reynolds_cl: f64, msq_cl: f64 },
    /// CL prescribed; α floats.
    Lift,
}

/// The assembled Newton system.
pub struct System {
    rows: usize,
    /// `A`: each row's equations w.r.t. its own `(S, θ)`.
    diag: Vec<Block>,
    /// `B`: each row's equations w.r.t. the previous station's `(S, θ)`.
    lower: Vec<Block>,
    /// Dense mass-coupling columns: `mass[row * rows + col][equation]`.
    mass: Vec<[f64; 3]>,
    /// Residual and parameter right-hand sides, replaced by the solution.
    rhs: Vec<Block>,
    /// `Z`: the first wake row w.r.t. the upper-surface TE station's `(S, θ)`.
    wake_coupling: Block,
    te_upper_row: usize,
    first_wake_row: usize,
}

impl System {
    fn zeros(rows: usize) -> Self {
        Self {
            rows,
            diag: vec![[[0.0; 2]; 3]; rows],
            lower: vec![[[0.0; 2]; 3]; rows],
            mass: vec![[0.0; 3]; rows * rows],
            rhs: vec![[[0.0; 2]; 3]; rows],
            wake_coupling: [[0.0; 2]; 3],
            te_upper_row: 0,
            first_wake_row: 0,
        }
    }

    #[inline]
    fn mass(&self, row: usize, col: usize) -> &[f64; 3] {
        &self.mass[row * self.rows + col]
    }

    #[inline]
    fn mass_mut(&mut self, row: usize, col: usize) -> &mut [f64; 3] {
        &mut self.mass[row * self.rows + col]
    }

    /// Solution `(ΔS, Δθ, Δm)` of `row` for the residual right-hand side, and
    /// its sensitivity to the floating parameter.  Valid after [`System::solve`].
    pub fn delta(&self, row: usize) -> ([f64; 3], [f64; 3]) {
        let d = &self.rhs[row];
        ([d[0][RESIDUAL], d[1][RESIDUAL], d[2][RESIDUAL]], [d[0][PARAMETER], d[1][PARAMETER], d[2][PARAMETER]])
    }
}
