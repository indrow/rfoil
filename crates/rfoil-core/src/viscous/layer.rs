//! Boundary-layer state and its mapping onto the panel nodes.
//!
//! The layer is split at the stagnation point into two sides.  The upper
//! side runs over the upper surface to the TE; the lower side runs over the
//! lower surface, to the TE and on along the wake.  Station 0 of each side is
//! a virtual point at the stagnation point itself.

use super::closure::kinematic_shape;
use super::station::{Primary, StationOutputs};
use crate::inviscid::Panels;
use crate::numeric::spline::{slopes, End, Spline};
use crate::numeric::Matrix;

pub const UPPER: usize = 0;
pub const LOWER: usize = 1;

/// Trailing-edge "flap" length relative to the TE gap, shaping the wake gap.
const WAKE_GAP_LENGTH_RATIO: f64 = 2.5;
/// Smallest edge speed allowed near the stagnation point.
const MIN_UE: f64 = 1.0e-7;

/// One boundary-layer station.
#[derive(Clone, Copy, Debug, Default)]
pub struct Node {
    /// Panel node this station sits on.
    pub panel: usize,
    /// +1 on the upper side, −1 on the lower side (turns panel speed into edge speed).
    pub sign: f64,
    /// Arc length from the stagnation point.
    pub xi: f64,
    /// Wake gap thickness (zero on the airfoil).
    pub wake_gap: f64,
    /// Viscous and inviscid edge speeds, and `dUe_inv/dα`.
    pub ue: f64,
    pub ue_inv: f64,
    pub ue_inv_alpha: f64,
    /// Amplification `N` while laminar, `√Cτ` once turbulent.
    pub ampl_or_shear: f64,
    pub theta: f64,
    pub dstar: f64,
    /// Mass defect `Ue·δ*`.
    pub mass: f64,
    pub outputs: StationOutputs,
    /// Row of this station in the global Newton system.
    pub row: usize,
}

impl Node {
    /// Primary variables for building a [`super::station::Station`].
    pub fn primary(&self, laminar: bool) -> Primary {
        Primary {
            xi: self.xi,
            ampl: if laminar { self.ampl_or_shear } else { 0.0 },
            shear: if laminar { 0.0 } else { self.ampl_or_shear },
            theta: self.theta,
            dstar: self.dstar,
            wake_gap: self.wake_gap,
            ue: self.ue,
        }
    }
}

/// Transition point found on a side.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TransitionPoint {
    pub xi: f64,
    pub x_over_c: f64,
    pub y_over_c: f64,
    pub forced: bool,
}

/// Stations along one side.
#[derive(Clone, Debug, Default)]
pub struct Side {
    pub nodes: Vec<Node>,
    /// Index of the trailing-edge station.
    pub te: usize,
    /// Index of the first turbulent station.
    pub first_turbulent: usize,
    pub transition: TransitionPoint,
}

impl Side {
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    pub fn is_wake(&self, i: usize) -> bool {
        i > self.te
    }

    pub fn is_laminar(&self, i: usize) -> bool {
        i < self.first_turbulent && i <= self.te
    }

    /// The real stations (all but the virtual stagnation station 0).
    pub fn stations(&self) -> impl Iterator<Item = (usize, &Node)> {
        self.nodes.iter().enumerate().skip(1)
    }
}

/// Stagnation point location on the panel contour.
#[derive(Clone, Copy, Debug, Default)]
pub struct Stagnation {
    /// The stagnation point lies between panel nodes `panel` and `panel + 1`.
    pub panel: usize,
    pub s: f64,
    /// Sensitivity of `s` to the vorticity at the two bracketing nodes.
    pub ds_dgamma_upper: f64,
    pub ds_dgamma_lower: f64,
}

impl Stagnation {
    /// Locates the zero crossing of the surface vorticity.
    pub fn find(gamma: &[f64], s: &[f64]) -> Self {
        let n = gamma.len();
        let i = (0..n - 1).find(|&i| gamma[i] >= 0.0 && gamma[i + 1] < 0.0).unwrap_or(n / 2);
        let dgamma = gamma[i + 1] - gamma[i];
        let ds = s[i + 1] - s[i];
        // interpolate from the node with the smaller |γ| to minimise round-off
        let s_stag = if gamma[i] < -gamma[i + 1] {
            s[i] - ds * (gamma[i] / dgamma)
        } else {
            s[i + 1] - ds * (gamma[i + 1] / dgamma)
        }
        .clamp(s[i] + 1.0e-7, s[i + 1] - 1.0e-7);
        Self {
            panel: i,
            s: s_stag,
            ds_dgamma_upper: (s_stag - s[i + 1]) / dgamma,
            ds_dgamma_lower: (s[i] - s_stag) / dgamma,
        }
    }
}

/// The complete boundary layer.
#[derive(Clone, Debug, Default)]
pub struct BoundaryLayer {
    pub sides: [Side; 2],
    pub stagnation: Stagnation,
    /// Number of rows in the global Newton system.
    pub rows: usize,
}

impl BoundaryLayer {
    /// Lays out stations on both sides of the stagnation point.
    pub fn new(p: &Panels, stagnation: Stagnation) -> Self {
        let n = p.n;
        let station = |panel, sign| Node { panel, sign, ..Node::default() };
        let virtual_station = std::iter::once(Node::default());

        let upper: Vec<Node> =
            virtual_station.clone().chain((0..=stagnation.panel).rev().map(|i| station(i, 1.0))).collect();
        let lower: Vec<Node> = virtual_station
            .chain((stagnation.panel + 1..n).map(|i| station(i, -1.0)))
            .chain((n..p.len()).map(|i| station(i, -1.0)))
            .collect();

        let te = [upper.len() - 1, n - stagnation.panel - 1];
        let side = |nodes, te| Side { nodes, te, first_turbulent: te, ..Side::default() };
        let mut bl = Self { sides: [side(upper, te[UPPER]), side(lower, te[LOWER])], stagnation, rows: 0 };
        bl.number_rows();
        bl.set_arc_lengths(p);
        bl
    }

    /// All real stations as `(side, index, node)`.
    pub fn stations(&self) -> impl Iterator<Item = (usize, usize, &Node)> {
        self.sides.iter().enumerate().flat_map(|(side, s)| s.stations().map(move |(i, node)| (side, i, node)))
    }

    pub fn stations_mut(&mut self) -> impl Iterator<Item = &mut Node> {
        self.sides.iter_mut().flat_map(|side| side.nodes.iter_mut().skip(1))
    }

    /// The trailing-edge station of each side.
    pub fn te_nodes(&self) -> [Node; 2] {
        self.sides.each_ref().map(|s| s.nodes[s.te])
    }

    fn number_rows(&mut self) {
        let mut row = 0;
        for node in self.stations_mut() {
            node.row = row;
            row += 1;
        }
        self.rows = row;
    }

    /// Arc lengths from the stagnation point and wake-gap thickness.
    pub fn set_arc_lengths(&mut self, p: &Panels) {
        let xi_min = 1.0e-7 * (p.s[p.n - 1] - p.s[0]);
        let s_stag = self.stagnation.s;
        for (side_index, side) in self.sides.iter_mut().enumerate() {
            let direction = if side_index == UPPER { -1.0 } else { 1.0 };
            for node in &mut side.nodes[1..=side.te] {
                node.xi = (direction * (p.s[node.panel] - s_stag)).max(xi_min);
            }
        }

        // the wake continues from the lower-surface TE station
        let lower = &mut self.sides[LOWER];
        let te = lower.te;
        for i in te + 1..lower.len() {
            let panel = lower.nodes[i].panel;
            let step = if i == te + 1 { 0.0 } else { p.distance(panel, panel - 1) };
            lower.nodes[i].xi = lower.nodes[i - 1].xi + step;
        }
        set_wake_gap(lower, p);
    }

    /// Combined trailing-edge `(√Cτ, θ, δ*)` that start the wake.  The TE base
    /// thickness `normal_gap` is added to the displacement thickness.
    pub fn trailing_edge_values(&self, normal_gap: f64) -> (f64, f64, f64) {
        let [u, l] = self.te_nodes();
        let theta = u.theta + l.theta;
        let shear = (u.ampl_or_shear * u.theta + l.ampl_or_shear * l.theta) / theta;
        (shear, theta, u.dstar + l.dstar + normal_gap)
    }

    /// Arc length from the stagnation point of a trip at chordwise station
    /// `x_over_c` on `side` (`x_over_c >= 1` means no trip: the TE).
    pub fn trip_xi(&self, p: &Panels, side: usize, x_over_c: f64) -> f64 {
        let s = &self.sides[side];
        let xi_te = s.nodes[s.te].xi;
        if x_over_c >= 1.0 {
            return xi_te;
        }
        let n = p.n;
        let arc = &p.s[..n];
        let chordwise: Vec<f64> = (0..n).map(|i| p.chord_coordinates(p.x[i], p.y[i]).0).collect();
        let chordwise_slope = slopes(&chordwise, arc, End::ZeroJerk, End::ZeroJerk);
        let (s_le, s_end) = (p.le.s, if side == UPPER { arc[0] } else { arc[n - 1] });
        let guess = s_le + (s_end - s_le) * x_over_c;
        let s_trip = Spline::new(&chordwise, &chordwise_slope, arc).invert(x_over_c, guess);
        let xi = if side == UPPER { self.stagnation.s - s_trip } else { s_trip - self.stagnation.s };
        if xi < 0.0 {
            xi_te
        } else {
            xi.min(xi_te)
        }
    }

    /// Inviscid edge speeds from the panel speeds.
    pub fn set_inviscid_speeds(&mut self, q_inv: &[f64], q_inv_alpha: &[f64]) {
        for node in self.stations_mut() {
            node.ue_inv = node.sign * q_inv[node.panel];
            node.ue_inv_alpha = node.sign * q_inv_alpha[node.panel];
        }
    }

    /// Panel speeds from the viscous edge speeds.
    pub fn write_panel_speeds(&self, q: &mut [f64]) {
        for (_, _, node) in self.stations() {
            q[node.panel] = node.sign * node.ue;
        }
    }

    /// Edge-speed sensitivity `dUe/dm` of `node` to the mass defect of `other`.
    pub fn ue_mass_coupling(node: &Node, other: &Node, dij: &Matrix) -> f64 {
        -node.sign * other.sign * dij[(node.panel, other.panel)]
    }

    /// Edge speed of every station implied by the inviscid speed plus the
    /// displacement effect of all mass defects (`Ue = Ue_inv + D·m`),
    /// indexed `[side][station]`.
    pub fn displacement_speeds(&self, dij: &Matrix) -> [Vec<f64>; 2] {
        self.sides.each_ref().map(|side| {
            side.nodes
                .iter()
                .enumerate()
                .map(|(i, node)| {
                    if i == 0 {
                        return 0.0;
                    }
                    let induced: f64 = self
                        .stations()
                        .map(|(_, _, other)| Self::ue_mass_coupling(node, other, dij) * other.mass)
                        .sum();
                    node.ue_inv + induced
                })
                .collect()
        })
    }

    /// Relocates the stagnation point after the vorticity `gamma` has changed,
    /// shifting the station variables if it moved to a different panel.
    /// Speeds below a small threshold are raised, also in `q` and `gamma`.
    pub fn move_stagnation(
        &mut self,
        p: &Panels,
        gamma: &mut [f64],
        q: &mut [f64],
        q_inv: &[f64],
        q_inv_alpha: &[f64],
    ) {
        let old_panel = self.stagnation.panel;
        let stagnation = Stagnation::find(&gamma[..p.n], &p.s[..p.n]);
        if stagnation.panel == old_panel {
            self.stagnation = stagnation;
            self.set_arc_lengths(p);
        } else {
            let mut moved = Self::new(p, stagnation);
            moved.set_inviscid_speeds(q_inv, q_inv_alpha);
            let shift = stagnation.panel as isize - old_panel as isize;
            // the side gaining stations grows at its start, the other loses them
            let (grow, shrink) = if shift > 0 { (UPPER, LOWER) } else { (LOWER, UPPER) };
            let k = shift.unsigned_abs();
            moved.sides[grow].first_turbulent = self.sides[grow].first_turbulent + k;
            moved.sides[shrink].first_turbulent = self.sides[shrink].first_turbulent.saturating_sub(k);
            moved.sides[grow].transition = self.sides[grow].transition;
            moved.sides[shrink].transition = self.sides[shrink].transition;

            copy_shifted(&self.sides[grow], &mut moved.sides[grow], |i| i.checked_sub(k).filter(|&j| j >= 1));
            extrapolate_to_stagnation(&mut moved.sides[grow], k);
            copy_shifted(&self.sides[shrink], &mut moved.sides[shrink], |i| Some(i + k));

            for node in moved.stations_mut() {
                if node.ue <= MIN_UE {
                    node.ue = MIN_UE;
                    q[node.panel] = node.sign * MIN_UE;
                    gamma[node.panel] = node.sign * MIN_UE;
                }
            }
            *self = moved;
        }
        for node in self.stations_mut() {
            node.mass = node.dstar * node.ue;
        }
    }
}

/// Copies station variables from `from` into `to`, station `i` taking
/// the values of `from` station `source(i)` when that exists.
fn copy_shifted(from: &Side, to: &mut Side, source: impl Fn(usize) -> Option<usize>) {
    for i in 1..to.len() {
        if let Some(src) = source(i).and_then(|j| from.nodes.get(j)) {
            let node = &mut to.nodes[i];
            (node.ampl_or_shear, node.theta, node.dstar, node.ue) = (src.ampl_or_shear, src.theta, src.dstar, src.ue);
        }
    }
}

/// Fills the `k` stations nearest the stagnation point from the first shifted station.
fn extrapolate_to_stagnation(side: &mut Side, k: usize) {
    let Some(&anchor) = side.nodes.get(k + 1) else { return };
    let due_dxi = anchor.ue / anchor.xi;
    for node in &mut side.nodes[1..=k] {
        (node.ampl_or_shear, node.theta, node.dstar) = (anchor.ampl_or_shear, anchor.theta, anchor.dstar);
        node.ue = due_dxi * node.xi;
    }
}

/// Wake gap: a cubic "dead air" thickness decaying behind a blunt TE.
fn set_wake_gap(side: &mut Side, p: &Panels) {
    let te = side.te;
    let last = p.n - 1;
    // sine of the angle between the two TE tangents sets the initial gap slope
    let (t0, t1) = ((p.dx_ds[0], p.dy_ds[0]), (p.dx_ds[last], p.dy_ds[last]));
    let cross = (t0.0 * t1.1 - t0.1 * t1.0) / ((t0.0.powi(2) + t0.1.powi(2)) * (t1.0.powi(2) + t1.1.powi(2))).sqrt();
    let limit = 3.0 / WAKE_GAP_LENGTH_RATIO;
    let slope = (cross / (1.0 - cross * cross).sqrt()).clamp(-limit, limit);
    let a = 3.0 + WAKE_GAP_LENGTH_RATIO * slope;
    let b = -2.0 - WAKE_GAP_LENGTH_RATIO * slope;
    let base = p.te.normal_gap;
    let xi_te = side.nodes[te].xi;

    for node in side.nodes.iter_mut().skip(te + 1) {
        node.wake_gap = if p.te.sharp {
            0.0
        } else {
            let z = 1.0 - (node.xi - xi_te) / (WAKE_GAP_LENGTH_RATIO * base);
            if z >= 0.0 {
                base * (a + b * z) * z * z
            } else {
                0.0
            }
        };
    }
}

/// Raises `δ*` so that the kinematic shape parameter is at least `hk_min`.
pub fn limit_dstar(dstar: f64, theta: f64, msq: f64, hk_min: f64) -> f64 {
    let (hk, hk_h, _) = kinematic_shape(dstar / theta, msq);
    dstar + (hk_min - hk).max(0.0) / hk_h * theta
}

/// Minimum kinematic shape parameter on the airfoil and in the wake.
pub fn hk_floor(wake: bool) -> f64 {
    if wake {
        1.00005
    } else {
        1.02
    }
}
