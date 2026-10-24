//! Under-relaxed application of the Newton step.

use super::{Constraint, System};
use crate::inviscid::{forces, Compressibility, Panels};
use crate::numeric::Matrix;
use crate::viscous::layer::{hk_floor, limit_dstar, BoundaryLayer};
use crate::viscous::params::Freestream;

/// Largest α change per iteration when CL is prescribed (0.5°).
const MAX_ALPHA_STEP: f64 = 0.5 * std::f64::consts::PI / 180.0;
/// Largest CL increase per iteration when α is prescribed.
const MAX_CL_STEP: f64 = 0.5;
/// Limits of the normalised change of any boundary-layer variable.
const MAX_GROWTH: f64 = 1.5;
const MAX_SHRINK: f64 = -0.5;
/// Edge-speed change counted as a unit normalised change.
const UE_CHANGE_SCALE: f64 = 0.25;

/// Summary of one Newton update.
#[derive(Clone, Copy, Debug, Default)]
pub struct UpdateReport {
    /// RMS of the normalised changes.
    pub rms: f64,
    /// Largest normalised change.
    pub max: f64,
    /// Under-relaxation factor applied.
    pub relaxation: f64,
    /// Relaxed change of the floating global parameter (CL or α).
    pub parameter_change: f64,
}

/// What [`apply`] needs to compute the new lift.
pub struct UpdateInputs<'a> {
    pub panels: &'a Panels,
    pub dij: &'a Matrix,
    pub constraint: Constraint,
    /// Current CL (α prescribed) or target CL (CL prescribed).
    pub cl: f64,
    pub alpha: f64,
    pub speed: f64,
    pub compressibility: Compressibility,
    pub mach_cl: f64,
    /// Lower limit on the CL change per step (relevant when Mach varies with CL).
    pub min_cl_step: f64,
    pub freestream: Freestream,
}

/// Unrelaxed Newton change of one station.
#[derive(Clone, Copy, Debug, Default)]
struct StationChange {
    ampl_or_shear: f64,
    theta: f64,
    dstar: f64,
    ue: f64,
}

/// Unrelaxed edge speeds after the step, and their sensitivity to the
/// floating parameter, indexed `[side][station]`.
struct EdgeSpeeds {
    ue: [Vec<f64>; 2],
    ue_parameter: [Vec<f64>; 2],
}

/// Applies the solved Newton changes with under-relaxation.
pub fn apply(bl: &mut BoundaryLayer, sys: &System, inp: &UpdateInputs) -> UpdateReport {
    let speeds = new_edge_speeds(bl, sys, inp);
    let (parameter_step, mut relaxation) = parameter_change(bl, &speeds, inp);
    let changes = station_changes(bl, sys, &speeds, parameter_step);

    // normalised changes: statistics, and the relaxation keeping them in bounds
    let (mut sum_sq, mut max, mut count) = (0.0_f64, 0.0_f64, 0);
    for (side, side_changes) in bl.sides.iter().zip(&changes) {
        count += side.len();
        for ((i, node), change) in side.stations().zip(&side_changes[1..]) {
            let shear_scale = if i < side.first_turbulent { 10.0 } else { node.ampl_or_shear };
            let normalised = [
                change.ampl_or_shear / shear_scale,
                change.theta / node.theta,
                change.dstar / node.dstar,
                change.ue.abs() / UE_CHANGE_SCALE,
            ];
            for d in normalised {
                sum_sq += d * d;
                if d.abs() > max.abs() {
                    max = d;
                }
                if relaxation * d > MAX_GROWTH {
                    relaxation = MAX_GROWTH / d;
                }
                if relaxation * d < MAX_SHRINK {
                    relaxation = MAX_SHRINK / d;
                }
            }
        }
    }

    for (side, side_changes) in bl.sides.iter_mut().zip(&changes) {
        let (te, first_turbulent) = (side.te, side.first_turbulent);
        for (i, change) in side_changes.iter().enumerate().skip(1) {
            let node = &mut side.nodes[i];
            node.ampl_or_shear += relaxation * change.ampl_or_shear;
            node.theta += relaxation * change.theta;
            node.dstar += relaxation * change.dstar;
            node.ue += relaxation * change.ue;
            if i >= first_turbulent {
                node.ampl_or_shear = node.ampl_or_shear.min(0.25);
            }
            let msq = inp.freestream.edge_msq(node.ue);
            node.dstar = limit_dstar(node.dstar - node.wake_gap, node.theta, msq, hk_floor(i > te)) + node.wake_gap;
            node.mass = node.dstar * node.ue;
        }
        // no isolated patches of reversed edge speed on the surface
        for i in 2..=te {
            if side.nodes[i - 1].ue > 0.0 && side.nodes[i].ue <= 0.0 {
                side.nodes[i].ue = side.nodes[i - 1].ue;
                side.nodes[i].mass = side.nodes[i].dstar * side.nodes[i].ue;
            }
        }
    }

    let rms = (sum_sq / (4.0 * count as f64)).sqrt();
    UpdateReport { rms, max, relaxation, parameter_change: relaxation * parameter_step }
}

/// Full (unrelaxed) changes of every station, indexed `[side][station]`
/// with a zero entry for the virtual station 0.
fn station_changes(
    bl: &BoundaryLayer,
    sys: &System,
    speeds: &EdgeSpeeds,
    parameter_step: f64,
) -> [Vec<StationChange>; 2] {
    let mut changes = [Vec::new(), Vec::new()];
    for (side, s) in bl.sides.iter().enumerate() {
        changes[side] = s
            .nodes
            .iter()
            .enumerate()
            .map(|(i, node)| {
                if i == 0 {
                    return StationChange::default();
                }
                let (residual, parameter) = sys.delta(node.row);
                let dmass = residual[2] - parameter_step * parameter[2];
                let due = speeds.ue[side][i] + parameter_step * speeds.ue_parameter[side][i] - node.ue;
                StationChange {
                    ampl_or_shear: residual[0] - parameter_step * parameter[0],
                    theta: residual[1] - parameter_step * parameter[1],
                    dstar: (dmass - node.dstar * due) / node.ue,
                    ue: due,
                }
            })
            .collect();
    }
    changes
}

/// Edge speeds implied by the solved mass changes.
fn new_edge_speeds(bl: &BoundaryLayer, sys: &System, inp: &UpdateInputs) -> EdgeSpeeds {
    let alpha_floats = inp.constraint == Constraint::Lift;
    let stations: Vec<_> = bl.stations().map(|(_, _, node)| *node).collect();
    let mut speeds = EdgeSpeeds { ue: [Vec::new(), Vec::new()], ue_parameter: [Vec::new(), Vec::new()] };
    for (side, s) in bl.sides.iter().enumerate() {
        for (i, node) in s.nodes.iter().enumerate() {
            let (mut du, mut du_parameter) = (0.0, 0.0);
            if i > 0 {
                for other in &stations {
                    let (residual, parameter) = sys.delta(other.row);
                    let ue_mass = BoundaryLayer::ue_mass_coupling(node, other, inp.dij);
                    du += ue_mass * (other.mass + residual[2]);
                    du_parameter -= ue_mass * parameter[2];
                }
            }
            speeds.ue[side].push(node.ue_inv + du);
            speeds.ue_parameter[side].push(if alpha_floats { node.ue_inv_alpha } else { 0.0 } + du_parameter);
        }
    }
    speeds
}

/// Newton change of the floating parameter (CL or α) implied by the lift of
/// the new surface speeds, and the relaxation needed to limit it.
fn parameter_change(bl: &BoundaryLayer, speeds: &EdgeSpeeds, inp: &UpdateInputs) -> (f64, f64) {
    let p = inp.panels;
    let (mut q, mut q_parameter) = (vec![0.0; p.n], vec![0.0; p.n]);
    for (side, s) in bl.sides.iter().enumerate() {
        for i in 1..=s.te {
            let node = s.nodes[i];
            q[node.panel] = node.sign * speeds.ue[side][i];
            q_parameter[node.panel] = node.sign * speeds.ue_parameter[side][i];
        }
    }
    let cond = forces::Condition { alpha: inp.alpha, speed: inp.speed, compressibility: inp.compressibility };
    // only the lift is needed, so the moment reference is irrelevant
    let loads = forces::integrate(&p.x[..p.n], &p.y[..p.n], &q, &q_parameter, &cond, (0.25, 0.0));

    let (step, lo, hi) = match inp.constraint {
        Constraint::Alpha { .. } => {
            let dcl_dmsq = loads.cl_msq * 2.0 * inp.compressibility.mach * inp.mach_cl;
            ((loads.cl - inp.cl) / (1.0 - loads.cl_speed - dcl_dmsq), inp.min_cl_step, MAX_CL_STEP)
        }
        Constraint::Lift => {
            ((loads.cl - inp.cl) / (-loads.cl_speed - loads.cl_geometry), -MAX_ALPHA_STEP, MAX_ALPHA_STEP)
        }
    };
    let mut relaxation = 1.0_f64;
    if relaxation * step > hi {
        relaxation = hi / step;
    }
    if relaxation * step < lo {
        relaxation = lo / step;
    }
    (step, relaxation)
}
