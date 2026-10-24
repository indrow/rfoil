//! The inviscid panel system and the source-influence matrix.
//!
//! The airfoil streamfunction must equal a constant `ψ₀` at every node, with a
//! Kutta condition closing the system.  Solving once for onset flows at
//! α = 0° and 90° gives two unit vorticity distributions that superpose to
//! any angle of attack.

use super::influence::{self, Probe, UniformFlow};
use super::panels::Panels;
use crate::numeric::{Lu, Matrix};

/// Control point distance ahead of a sharp TE, as a fraction of the shorter TE panel.
const BISECTOR_OFFSET: f64 = 0.1;

/// Factored inviscid system of an airfoil paneling.
#[derive(Clone, Debug)]
pub struct InviscidSystem {
    lu: Lu,
    /// `−dψ/dσ` of the airfoil sources in each equation (`n+1 × n`).
    source_terms: Matrix,
    /// Unit vorticity distributions for α = 0°, 90°; entry `n` is `ψ₀`.
    pub unit_gamma: [Vec<f64>; 2],
}

impl InviscidSystem {
    pub fn new(p: &Panels, speed: f64) -> Self {
        let n = p.n;
        let zero = vec![0.0; n];
        let no_vorticity = [zero.clone(), zero.clone()];
        let onset = UniformFlow { speed, alpha: 0.0 };

        let mut matrix = Matrix::zeros(n + 1, n + 1);
        let mut source_terms = Matrix::zeros(n + 1, n);
        let mut unit_rhs = [vec![0.0; n + 1], vec![0.0; n + 1]];

        // ψ(γ, σ) − ψ₀ = 0 at every node
        for i in 0..n {
            let (v, sources) = influence::airfoil(p, &zero, &no_vorticity, onset, &Probe::at_node(p, i), true);
            let sources = sources.expect("source influence requested");
            matrix.row_mut(i)[..n].copy_from_slice(&v.dpsi_dgamma);
            matrix[(i, n)] = -1.0;
            for (term, d) in source_terms.row_mut(i).iter_mut().zip(&sources.dpsi_dsigma) {
                *term = -d;
            }
            unit_rhs[0][i] = -speed * p.y[i];
            unit_rhs[1][i] = speed * p.x[i];
        }

        // Kutta condition: γ_first + γ_last = 0
        matrix[(n, 0)] = 1.0;
        matrix[(n, n - 1)] = 1.0;

        // a sharp TE instead has zero internal speed along its bisector
        if p.te.sharp {
            let probe = bisector_probe(p);
            let (v, sources) = influence::airfoil(p, &zero, &no_vorticity, onset, &probe, true);
            let sources = sources.expect("source influence requested");
            let last = n - 1;
            matrix.row_mut(last)[..n].copy_from_slice(&v.dq_dgamma);
            matrix[(last, n)] = 0.0;
            for (term, d) in source_terms.row_mut(last).iter_mut().zip(&sources.dq_dsigma) {
                *term = -d;
            }
            unit_rhs[0][last] = -probe.ny;
            unit_rhs[1][last] = probe.nx;
        }

        let lu = Lu::factor(matrix);
        for rhs in &mut unit_rhs {
            lu.solve(rhs);
        }
        Self { lu, source_terms, unit_gamma: unit_rhs }
    }

    /// Vorticity and its α-derivative on the `n` airfoil nodes at angle of attack `alpha`.
    pub fn vorticity(&self, n: usize, alpha: f64) -> (Vec<f64>, Vec<f64>) {
        let (sa, ca) = alpha.sin_cos();
        let [g0, g90] = &self.unit_gamma;
        let gamma = (0..n).map(|i| ca * g0[i] + sa * g90[i]).collect();
        let gamma_alpha = (0..n).map(|i| -sa * g0[i] + ca * g90[i]).collect();
        (gamma, gamma_alpha)
    }

    /// Unit solutions restricted to the `n` airfoil nodes.
    pub fn unit_vorticity(&self, n: usize) -> [Vec<f64>; 2] {
        [self.unit_gamma[0][..n].to_vec(), self.unit_gamma[1][..n].to_vec()]
    }

    /// Matrix `D` with `dq_i/dσ_j` over all airfoil and wake nodes, including
    /// the indirect effect of the sources on the airfoil vorticity.
    pub fn source_influence(&self, p: &Panels) -> Matrix {
        let (n, total) = (p.n, p.len());
        let mut d = Matrix::zeros(total, total);
        let mut column = vec![0.0; n + 1];

        // airfoil vorticity (= surface speed) per unit source j, from the system's rhs
        let mut solve_column = |d: &mut Matrix, j: usize, rhs: &dyn Fn(usize) -> f64| {
            column.iter_mut().enumerate().for_each(|(i, c)| *c = rhs(i));
            self.lu.solve(&mut column);
            for i in 0..n {
                d[(i, j)] = column[i];
            }
        };

        // how the airfoil sources change the airfoil vorticity
        for j in 0..n {
            solve_column(&mut d, j, &|i| self.source_terms[(i, j)]);
        }

        // how the wake sources change the airfoil vorticity
        let wake_rows: Vec<Vec<f64>> =
            (0..n).map(|i| influence::wake_sources(p, &Probe::at_node(p, i)).dpsi_dsigma).collect();
        // (neither the Kutta row nor the sharp-TE closure feel sources directly)
        let unaffected = |i: usize| i == n || (p.te.sharp && i == n - 1);
        for j in n..total {
            solve_column(&mut d, j, &|i| if unaffected(i) { 0.0 } else { -wake_rows[i][j] });
        }

        // how all sources change the wake speeds, directly and through the airfoil vorticity
        let zero = vec![0.0; n];
        let no_vorticity = [zero.clone(), zero.clone()];
        let still = UniformFlow { speed: 0.0, alpha: 0.0 };
        for i in n..total {
            let probe = Probe::at_node(p, i);
            let (v, sources) = influence::airfoil(p, &zero, &no_vorticity, still, &probe, true);
            let sources = sources.expect("source influence requested");
            let wake = influence::wake_sources(p, &probe);
            for j in 0..total {
                let direct = if j < n { sources.dq_dsigma[j] } else { wake.dq_dsigma[j] };
                let via_vorticity: f64 = (0..n).map(|k| v.dq_dgamma[k] * d[(k, j)]).sum();
                d[(i, j)] = direct + via_vorticity;
            }
        }

        // the first wake point coincides with the TE
        if total > n {
            let te_row = d.row(n - 1).to_vec();
            d.row_mut(n).copy_from_slice(&te_row);
        }
        d
    }
}

/// Probe on the TE bisector just ahead of a sharp trailing edge, pointing along it.
fn bisector_probe(p: &Panels) -> Probe {
    let last = p.n - 1;
    let upper_angle = (-p.dy_ds[0]).atan2(-p.dx_ds[0]);
    let lower_angle = atan_continuous(p.dy_ds[last], p.dx_ds[last], upper_angle);
    let (sin_bis, cos_bis) = (0.5 * (upper_angle + lower_angle)).sin_cos();
    let offset = BISECTOR_OFFSET * p.distance(0, 1).min(p.distance(last, last - 1));
    let (xte, yte) = p.te_midpoint;
    Probe { x: xte - offset * cos_bis, y: yte - offset * sin_bis, nx: -sin_bis, ny: cos_bis, node: None }
}

/// `atan2(y, x)` unwrapped to lie within ±π of `reference`.
pub fn atan_continuous(y: f64, x: f64, reference: f64) -> f64 {
    use std::f64::consts::{PI, TAU};
    let delta = y.atan2(x) - reference;
    reference + delta - TAU * ((delta + PI.copysign(delta)) / TAU).trunc()
}
