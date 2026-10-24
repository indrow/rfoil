//! Assembly of the Newton system from the interval linearisations.

use super::{Constraint, System, PARAMETER, RESIDUAL};
use crate::inviscid::Panels;
use crate::numeric::Matrix;
use crate::viscous::dual::var;
use crate::viscous::interval::{Linearization, Segment};
use crate::viscous::layer::{BoundaryLayer, TransitionPoint, LOWER, UPPER};
use crate::viscous::params::Context;
use crate::viscous::station::{Primary, Station};
use crate::viscous::transition::{self, Onset};

/// Linear dependence of a station quantity on every mass defect and on the
/// floating parameter, plus its "forced" change: the mismatch between its
/// current value and the one implied by `Ue_inv + D·m`.
#[derive(Clone)]
struct Sensitivity {
    mass: Vec<f64>,
    parameter: f64,
    forced: f64,
}

impl Sensitivity {
    fn zeros(rows: usize) -> Self {
        Self { mass: vec![0.0; rows], parameter: 0.0, forced: 0.0 }
    }
}

/// Edge-speed and displacement-thickness sensitivities of one interval end.
#[derive(Clone)]
struct Coupling {
    ue: Sensitivity,
    dstar: Sensitivity,
}

impl Coupling {
    fn zeros(rows: usize) -> Self {
        Self { ue: Sensitivity::zeros(rows), dstar: Sensitivity::zeros(rows) }
    }
}

/// How the stagnation point (and with it every ξ) moves with the edge
/// speeds of the first station on each side.
struct StagnationCoupling {
    /// Edge-speed sensitivities of the first station on each side.
    first_ue: [Sensitivity; 2],
    /// Change of this side's ξ per unit change of each of those speeds.
    xi_per_ue: [f64; 2],
}

/// Inputs shared by the whole assembly.
pub struct Assembly<'a> {
    pub panels: &'a Panels,
    pub dij: &'a Matrix,
    pub constraint: Constraint,
    pub ctx: [Context; 2],
}

/// Linearises the boundary-layer equations about the current state.
/// Also refreshes the transition locations and station outputs in `bl`.
pub fn assemble(bl: &mut BoundaryLayer, a: &Assembly) -> System {
    let rows = bl.rows;
    let mut sys = System::zeros(rows);
    let p = a.panels;

    // edge speeds implied by the current masses; differences from the actual
    // speeds are carried as forced changes
    let implied_ue = bl.displacement_speeds(a.dij);

    let ue_sensitivity = |bl: &BoundaryLayer, side: usize, i: usize| -> Sensitivity {
        let node = bl.sides[side].nodes[i];
        let mut mass = vec![0.0; rows];
        for (_, _, other) in bl.stations() {
            mass[other.row] = BoundaryLayer::ue_mass_coupling(&node, other, a.dij);
        }
        Sensitivity { mass, parameter: node.ue_inv_alpha, forced: node.ue - implied_ue[side][i] }
    };
    let first_ue = [ue_sensitivity(bl, UPPER, 1), ue_sensitivity(bl, LOWER, 1)];
    let te_ue = [UPPER, LOWER].map(|side| ue_sensitivity(bl, side, bl.sides[side].te).mass);

    sys.te_upper_row = bl.sides[UPPER].nodes[bl.sides[UPPER].te].row;
    sys.first_wake_row = bl.sides[LOWER].nodes[bl.sides[LOWER].te + 1].row;

    for side in 0..2 {
        let ctx = &a.ctx[side];
        let (te, first_turbulent) = (bl.sides[side].te, bl.sides[side].first_turbulent);
        record_transition_at_te(bl, p, side);

        let (ds_upper, ds_lower) = (bl.stagnation.ds_dgamma_upper, bl.stagnation.ds_dgamma_lower);
        let stagnation = StagnationCoupling {
            first_ue: first_ue.clone(),
            xi_per_ue: if side == UPPER { [ds_upper, -ds_lower] } else { [-ds_upper, ds_lower] },
        };

        let mut upstream_coupling = Coupling::zeros(rows);
        let mut upstream = Station::default();
        for i in 1..bl.sides[side].len() {
            let node = bl.sides[side].nodes[i];
            let wake = i > te;
            let laminar = i < first_turbulent && !wake;

            // δ* follows from the mass defect and the edge speed
            let dstar = node.mass / node.ue;
            let dstar_ue = -dstar / node.ue;
            let ue = ue_sensitivity(bl, side, i);
            let mut dstar_mass: Vec<f64> = ue.mass.iter().map(|u| dstar_ue * u).collect();
            dstar_mass[node.row] += 1.0 / node.ue;
            let dstar_sensitivity =
                Sensitivity { mass: dstar_mass, parameter: dstar_ue * ue.parameter, forced: dstar_ue * ue.forced };
            let coupling = Coupling { ue, dstar: dstar_sensitivity };

            let mut primary = Primary { dstar, ..node.primary(laminar) };
            if !laminar {
                primary.ampl = upstream.ampl.value;
            }
            let mut s2 = Station::new(&primary, &ctx.freestream);
            let mut onset: Option<Onset> = None;
            if i == first_turbulent && !wake {
                let check = transition::check(&upstream, &s2, ctx);
                s2.ampl.value = check.downstream_ampl;
                onset = check.onset;
            }

            let segment = if i == te + 1 {
                let (shear, theta, dstar) = bl.trailing_edge_values(p.te.normal_gap);
                upstream_coupling = trailing_edge_coupling(bl, &te_ue, &implied_ue, rows);
                Segment::TrailingEdge { shear, theta, dstar }
            } else {
                Segment::select(i == 1, wake, i > first_turbulent || wake, onset)
            };
            let lin = Linearization::from(segment.equations(&upstream, &mut s2, ctx));
            bl.sides[side].nodes[i].outputs = s2.outputs();

            fill_rows(&mut sys, node.row, &lin, [&upstream_coupling, &coupling], &stagnation, a.constraint);
            if let Segment::TrailingEdge { theta, .. } = segment {
                couple_wake_to_both_te(&mut sys, bl, node.row, &lin, theta);
            }
            if let Some(onset) = onset {
                record_transition(bl, p, side, &onset);
            }

            upstream = s2.as_upstream();
            upstream_coupling = coupling;
        }
    }
    sys
}

/// Upstream coupling of the first wake interval, whose upstream end is the
/// combination of both trailing-edge stations.
fn trailing_edge_coupling(
    bl: &BoundaryLayer,
    te_ue_mass: &[Vec<f64>; 2],
    implied_ue: &[Vec<f64>; 2],
    rows: usize,
) -> Coupling {
    let mut c = Coupling::zeros(rows);
    for (side, s) in bl.sides.iter().enumerate() {
        let node = s.nodes[s.te];
        let dstar_ue = -node.dstar / node.ue;
        for (d, u) in c.dstar.mass.iter_mut().zip(&te_ue_mass[side]) {
            *d += dstar_ue * u;
        }
        c.dstar.mass[node.row] += 1.0 / node.ue;
        c.dstar.parameter += dstar_ue * node.ue_inv_alpha;
        c.dstar.forced += dstar_ue * (node.ue - implied_ue[side][s.te]);
    }
    c
}

/// Writes the three Newton rows of station `row` from the interval
/// linearisation and the couplings of its upstream and downstream ends.
fn fill_rows(
    sys: &mut System,
    row: usize,
    lin: &Linearization,
    [c1, c2]: [&Coupling; 2],
    stagnation: &StagnationCoupling,
    constraint: Constraint,
) {
    use var::{AMPL_OR_SHEAR as S, DSTAR as D, THETA as T, UE as U, XI as X};
    let [le_upper, le_lower] = &stagnation.first_ue;
    let [xi_upper, xi_lower] = stagnation.xi_per_ue;

    for k in 0..3 {
        let (j1, j2) = (&lin.upstream[k], &lin.downstream[k]);
        // ξ enters through the stagnation-point position, and through the trip
        let dxi = j1[X] + j2[X] + lin.d_xi_trip[k];

        for col in 0..sys.rows {
            sys.mass_mut(row, col)[k] = j1[D] * c1.dstar.mass[col]
                + j1[U] * c1.ue.mass[col]
                + j2[D] * c2.dstar.mass[col]
                + j2[U] * c2.ue.mass[col]
                + dxi * (xi_upper * le_upper.mass[col] + xi_lower * le_lower.mass[col]);
        }
        sys.lower[row][k] = [j1[S], j1[T]];
        sys.diag[row][k] = [j2[S], j2[T]];

        sys.rhs[row][k][PARAMETER] = match constraint {
            Constraint::Alpha { reynolds_cl, msq_cl } => lin.d_reynolds[k] * reynolds_cl + lin.d_msq[k] * msq_cl,
            Constraint::Lift => {
                j1[U] * c1.ue.parameter
                    + j1[D] * c1.dstar.parameter
                    + j2[U] * c2.ue.parameter
                    + j2[D] * c2.dstar.parameter
                    + dxi * (xi_upper * le_upper.parameter + xi_lower * le_lower.parameter)
            }
        };
        sys.rhs[row][k][RESIDUAL] = lin.rhs[k]
            + j1[U] * c1.ue.forced
            + j1[D] * c1.dstar.forced
            + j2[U] * c2.ue.forced
            + j2[D] * c2.dstar.forced
            + dxi * (xi_upper * le_upper.forced + xi_lower * le_lower.forced);
    }
}

/// The first wake interval starts from the θ-weighted mean shear and the
/// summed θ of both TE stations; splits its upstream `(S, θ)` columns
/// between the upper TE station (`Z`) and the lower TE station (`B`).
fn couple_wake_to_both_te(sys: &mut System, bl: &BoundaryLayer, row: usize, lin: &Linearization, theta_te: f64) {
    let [u, l] = bl.te_nodes();
    let shear_te = (u.ampl_or_shear * u.theta + l.ampl_or_shear * l.theta) / theta_te;
    // ∂S_te/∂S and ∂S_te/∂θ of each TE station
    let (shear_d_shear_upper, shear_d_shear_lower) = (u.theta / theta_te, l.theta / theta_te);
    let (shear_d_theta_upper, shear_d_theta_lower) =
        ((u.ampl_or_shear - shear_te) / theta_te, (l.ampl_or_shear - shear_te) / theta_te);
    for k in 0..3 {
        let (d_shear, d_theta) = (lin.upstream[k][var::AMPL_OR_SHEAR], lin.upstream[k][var::THETA]);
        sys.wake_coupling[k] = [d_shear * shear_d_shear_upper, d_shear * shear_d_theta_upper + d_theta];
        sys.lower[row][k] = [d_shear * shear_d_shear_lower, d_shear * shear_d_theta_lower + d_theta];
    }
}

/// Default while a side stays laminar: transition at the trailing edge.
fn record_transition_at_te(bl: &mut BoundaryLayer, p: &Panels, side: usize) {
    let s = &mut bl.sides[side];
    let node = s.nodes[s.te];
    let (x_over_c, y_over_c) = p.chord_coordinates(p.x[node.panel], p.y[node.panel]);
    s.transition = TransitionPoint { xi: node.xi, x_over_c, y_over_c, forced: false };
}

fn record_transition(bl: &mut BoundaryLayer, p: &Panels, side: usize, onset: &Onset) {
    let s = if side == UPPER { bl.stagnation.s - onset.xi.value } else { bl.stagnation.s + onset.xi.value };
    let contour = p.curve();
    let (x_over_c, y_over_c) = p.chord_coordinates(contour.x_at(s), contour.y_at(s));
    bl.sides[side].transition = TransitionPoint { xi: onset.xi.value, x_over_c, y_over_c, forced: onset.forced };
}
