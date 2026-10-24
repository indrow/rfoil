//! Station-by-station marching of the boundary layer along each side.
//!
//! Each station is converged by a local 4×4 Newton iteration on
//! `(N or √Cτ, θ, δ*, Ue)`: the three interval equations plus one extra
//! constraint, which depends on the marching mode.

use super::closure::kinematic_shape;
use super::dual::{slot, var};
use super::interval::{Linearization, Segment};
use super::layer::{hk_floor, limit_dstar, BoundaryLayer, Node};
use super::params::Context;
use super::station::{Primary, Regime, Station};
use super::transition::{self, Onset};
use crate::numeric::solve_small;

/// How the edge speed is treated while marching.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// Ue is prescribed (switching to a prescribed Hk if the flow separates).
    /// Used to initialise the boundary layer from the inviscid solution.
    PrescribedUe,
    /// Ue and Hk move along a line normal to the local Ue–Hk characteristic,
    /// which never meets the Goldstein/Levy–Lees separation singularity.
    Mixed,
}

const MAX_ITERATIONS: usize = 25;
/// Weight of the Ue–Hk characteristic in mixed mode.
const SENSITIVITY_WEIGHT: f64 = 1000.0;
/// Largest Hk tolerated by the direct (prescribed-Ue) mode.
const HK_MAX_LAMINAR: f64 = 3.8;
const HK_MAX_TURBULENT: f64 = 2.5;
/// Initial `√Cτ` of a freshly turbulent station.
const INITIAL_SHEAR: f64 = 0.03;

/// Marches both sides of the boundary layer.  `context_for` supplies the
/// context (with that side's trip location) for each side.
pub fn march(bl: &mut BoundaryLayer, te_gap: f64, mode: Mode, context_for: impl Fn(&BoundaryLayer, usize) -> Context) {
    for side in 0..2 {
        let ctx = context_for(bl, side);
        Marcher::new(bl, side, ctx, te_gap, mode).run();
    }
}

/// Running guesses for the primary variables, carried between stations.
#[derive(Clone, Copy, Debug, Default)]
struct Guess {
    ampl: f64,
    shear: f64,
    theta: f64,
    dstar: f64,
    ue: f64,
}

/// Edge speed and shape parameter the mixed-mode constraint is anchored to.
#[derive(Clone, Copy, Debug, Default)]
struct Reference {
    ue: f64,
    hk: f64,
}

/// The local Newton system of one station in the downstream variables
/// `(N or √Cτ, θ, δ*, Ue)`: three interval equations plus one constraint.
#[derive(Clone, Copy)]
struct LocalSystem {
    a: [[f64; 4]; 4],
    rhs: [f64; 4],
}

impl LocalSystem {
    fn from_interval(lin: &Linearization) -> Self {
        let mut local = Self { a: [[0.0; 4]; 4], rhs: [0.0; 4] };
        for k in 0..3 {
            local.a[k].copy_from_slice(&lin.downstream[k][..4]);
            local.rhs[k] = lin.rhs[k];
        }
        local
    }

    fn set_constraint(&mut self, (row, rhs): ([f64; 4], f64)) {
        (self.a[3], self.rhs[3]) = (row, rhs);
    }

    /// The Newton changes of the four variables.
    fn solve(mut self) -> [f64; 4] {
        solve_small(&mut self.a, &mut self.rhs);
        self.rhs
    }
}

struct Marcher<'a> {
    bl: &'a mut BoundaryLayer,
    side: usize,
    ctx: Context,
    te_gap: f64,
    mode: Mode,
    guess: Guess,
    /// The converged previous station, as the upstream end of the next interval.
    upstream: Station,
    turbulent: bool,
    /// First turbulent station before this march started.
    previous_first_turbulent: usize,
    /// Smoothed `dUe/dHk` along the local characteristic (mixed mode).
    ue_hk_sensitivity: f64,
}

impl<'a> Marcher<'a> {
    fn new(bl: &'a mut BoundaryLayer, side: usize, ctx: Context, te_gap: f64, mode: Mode) -> Self {
        let s = &mut bl.sides[side];
        let previous_first_turbulent = s.first_turbulent;
        s.first_turbulent = s.te;
        Self {
            bl,
            side,
            ctx,
            te_gap,
            mode,
            guess: Guess::default(),
            upstream: Station::default(),
            turbulent: false,
            previous_first_turbulent,
            ue_hk_sensitivity: 0.0,
        }
    }

    fn run(mut self) {
        if self.mode == Mode::PrescribedUe {
            self.guess = self.thwaites_start();
        }
        for i in 1..self.bl.sides[self.side].len() {
            self.start_guess(i);
            let onset = self.converge(i);
            self.store(i, onset);
        }
    }

    fn node(&self, i: usize) -> Node {
        self.bl.sides[self.side].nodes[i]
    }

    fn te(&self) -> usize {
        self.bl.sides[self.side].te
    }

    fn first_turbulent(&self) -> usize {
        self.bl.sides[self.side].first_turbulent
    }

    fn is_wake(&self, i: usize) -> bool {
        self.bl.sides[self.side].is_wake(i)
    }

    /// Similarity-station θ from Thwaites' formula.
    fn thwaites_start(&self) -> Guess {
        let node = self.node(1);
        let due_dxi = node.ue / node.xi;
        let theta = (0.45 / (due_dxi * 6.0 * self.ctx.freestream.reynolds.value)).sqrt();
        Guess { ampl: 0.0, shear: INITIAL_SHEAR, theta, dstar: 2.2 * theta, ue: node.ue }
    }

    /// Mixed mode starts from the stored solution; direct mode carries the
    /// previous station's values forward.
    fn start_guess(&mut self, i: usize) {
        let node = self.node(i);
        self.guess.ue = node.ue;
        if self.mode == Mode::Mixed {
            self.guess.theta = node.theta;
            self.guess.dstar = node.dstar;
            if i < self.previous_first_turbulent {
                self.guess.ampl = node.ampl_or_shear;
                self.guess.shear = INITIAL_SHEAR;
            } else {
                self.guess.shear = if node.ampl_or_shear <= 0.0 { INITIAL_SHEAR } else { node.ampl_or_shear };
            }
            let floor = hk_floor(self.is_wake(i));
            self.guess.dstar = (self.guess.dstar - node.wake_gap).max(floor * self.guess.theta) + node.wake_gap;
        }
    }

    fn primary(&self, i: usize) -> Primary {
        let node = self.node(i);
        Primary {
            xi: node.xi,
            ampl: self.guess.ampl,
            shear: self.guess.shear,
            theta: self.guess.theta,
            dstar: self.guess.dstar,
            wake_gap: node.wake_gap,
            ue: self.guess.ue,
        }
    }

    /// Builds station `i` from the current guess and checks for transition.
    fn evaluate(&mut self, i: usize) -> (Station, Option<Onset>) {
        let mut s2 = Station::new(&self.primary(i), &self.ctx.freestream);
        let mut onset = None;
        if i > 1 && !self.turbulent {
            let check = transition::check(&self.upstream, &s2, &self.ctx);
            self.guess.ampl = check.downstream_ampl;
            s2.ampl.value = check.downstream_ampl;
            onset = check.onset;
            self.bl.sides[self.side].first_turbulent = if onset.is_some() { i } else { i + 2 };
        }
        (s2, onset)
    }

    fn segment(&self, i: usize, onset: Option<Onset>) -> Segment {
        if i == self.te() + 1 {
            let (shear, theta, dstar) = self.bl.trailing_edge_values(self.te_gap);
            Segment::TrailingEdge { shear, theta, dstar }
        } else {
            Segment::select(i == 1, self.is_wake(i), self.turbulent, onset)
        }
    }

    /// Converges station `i`, returning the transition onset if it lies in the interval.
    fn converge(&mut self, i: usize) -> Option<Onset> {
        let tolerance = if self.mode == Mode::Mixed { 5.0e-6 } else { 1.0e-5 };
        let mut inverse_target: Option<f64> = None;
        let mut reference = Reference::default();
        let mut max_change = 0.0;

        for iteration in 0..MAX_ITERATIONS {
            let (mut s2, onset) = self.evaluate(i);
            self.seed_shear(i, iteration, onset, &mut s2);

            let lin = Linearization::from(self.segment(i, onset).equations(&self.upstream, &mut s2, &self.ctx));
            if self.mode == Mode::Mixed && iteration == 0 {
                reference = self.reference(i, &s2);
            }

            let mut local = LocalSystem::from_interval(&lin);
            local.set_constraint(self.constraint(i, iteration, &s2, inverse_target, reference, &local));
            let delta = local.solve();

            let direct = self.mode == Mode::PrescribedUe && inverse_target.is_none();
            let turbulent_here = i >= self.first_turbulent();
            max_change = self.max_relative_change(&delta, turbulent_here, direct);
            let relaxation = if max_change > 0.3 { 0.3 / max_change } else { 1.0 };

            if direct && i != self.te() + 1 {
                // switch to inverse mode if the prescribed Ue would separate the flow
                let hk_max = if turbulent_here { HK_MAX_TURBULENT } else { HK_MAX_LAMINAR };
                if self.would_separate(&delta, relaxation, hk_max) {
                    inverse_target = Some(self.target_hk(i, &s2, onset, hk_max));
                    continue;
                }
            }

            if turbulent_here {
                self.guess.shear += relaxation * delta[0];
            } else if self.mode == Mode::Mixed {
                self.guess.ampl += relaxation * delta[0];
            }
            self.guess.theta += relaxation * delta[1];
            self.guess.dstar += relaxation * delta[2];
            if !direct {
                self.guess.ue += relaxation * delta[3];
            }

            self.apply_limits(i);
            if max_change <= tolerance {
                return onset;
            }
        }

        // not converged: keep a nearly converged answer, else extrapolate from upstream
        if max_change > 0.1 && i > 2 {
            self.extrapolate(i);
        }
        self.evaluate(i).1
    }

    /// Gives a station that just became turbulent a starting shear stress.
    fn seed_shear(&mut self, i: usize, iteration: usize, onset: Option<Onset>, s2: &mut Station) {
        let newly_turbulent = match self.mode {
            Mode::PrescribedUe => onset.is_some() && self.guess.shear <= 0.0,
            // a formerly laminar station that is now turbulent needs a fresh Cτ
            Mode::Mixed => iteration == 0 && i < self.previous_first_turbulent && (onset.is_some() || self.turbulent),
        };
        if newly_turbulent {
            self.guess.shear = match (self.mode, onset) {
                (Mode::Mixed, None) => self.node(i - 1).ampl_or_shear,
                _ => INITIAL_SHEAR,
            };
            s2.shear.value = self.guess.shear;
        }
    }

    /// Ue and Hk at the start of the mixed-mode iteration.  At a station that
    /// was laminar before but is newly in the transition region, Hk is taken
    /// from the upstream station instead.
    fn reference(&self, i: usize, s2: &Station) -> Reference {
        let mut hk = s2.hk.value;
        if i < self.first_turbulent() && i >= self.previous_first_turbulent {
            let prev = self.node(i - 1);
            hk = kinematic_shape(prev.dstar / prev.theta, self.ctx.freestream.edge_msq(prev.ue)).0;
        }
        Reference { ue: s2.ue.value, hk }
    }

    /// The fourth row `(coefficients, rhs)` of the local system, whose first
    /// three rows hold the interval equations.
    fn constraint(
        &mut self,
        i: usize,
        iteration: usize,
        s2: &Station,
        inverse_target: Option<f64>,
        reference: Reference,
        local: &LocalSystem,
    ) -> ([f64; 4], f64) {
        let down = |v| slot::downstream(v);
        let hk_partials =
            [0.0, s2.hk.partials[down(var::THETA)], s2.hk.partials[down(var::DSTAR)], s2.hk.partials[down(var::UE)]];
        let ue_partial = s2.ue.partials[down(var::UE)];

        match (self.mode, inverse_target) {
            (Mode::PrescribedUe, None) => ([0.0, 0.0, 0.0, 1.0], 0.0),
            (Mode::PrescribedUe, Some(hk_target)) => (hk_partials, hk_target - s2.hk.value),
            // the similarity and first wake stations keep their edge speed
            (Mode::Mixed, _) if i == 1 || i == self.te() + 1 => {
                ([0.0, 0.0, 0.0, ue_partial], reference.ue - s2.ue.value)
            }
            (Mode::Mixed, _) => {
                // dUe response to a unit dHk along the local characteristic
                let mut probe = *local;
                probe.set_constraint((hk_partials, 1.0));
                let response = probe.solve();
                let sensitivity = SENSITIVITY_WEIGHT * response[3] * reference.hk / reference.ue;
                if iteration < 5 {
                    self.ue_hk_sensitivity = sensitivity;
                } else if iteration < 15 {
                    self.ue_hk_sensitivity = 0.5 * (self.ue_hk_sensitivity + sensitivity);
                }

                let (s, hk_ref, ue_ref) = (self.ue_hk_sensitivity, reference.hk, reference.ue);
                (
                    [
                        0.0,
                        hk_partials[1] * hk_ref,
                        hk_partials[2] * hk_ref,
                        hk_partials[3] * hk_ref + s / ue_ref * ue_partial,
                    ],
                    -hk_ref * hk_ref * (s2.hk.value / hk_ref - 1.0) - s * (s2.ue.value / ue_ref - 1.0),
                )
            }
        }
    }

    /// Largest Newton change relative to the current guess.
    fn max_relative_change(&self, delta: &[f64; 4], turbulent: bool, direct: bool) -> f64 {
        let g = self.guess;
        let mut max = (delta[1] / g.theta).abs().max((delta[2] / g.dstar).abs());
        if !direct {
            max = max.max((delta[3] / g.ue).abs());
        }
        max.max(match (turbulent, self.mode) {
            (true, Mode::Mixed) => (delta[0] / (10.0 * g.shear)).abs(),
            (true, _) => (delta[0] / g.shear).abs(),
            (false, Mode::PrescribedUe) if direct => (delta[0] / 10.0).abs(),
            (false, _) => 0.0,
        })
    }

    /// Whether the relaxed step would push Hk past `hk_max`.
    fn would_separate(&self, delta: &[f64; 4], relaxation: f64, hk_max: f64) -> bool {
        let g = self.guess;
        let msq = self.ctx.freestream.edge_msq(g.ue);
        let h = (g.dstar + relaxation * delta[2]) / (g.theta + relaxation * delta[1]);
        kinematic_shape(h, msq).0 >= hk_max
    }

    /// Hk prescribed in inverse mode, following plausible separated-flow trends.
    fn target_hk(&self, i: usize, s2: &Station, onset: Option<Onset>, hk_max: f64) -> f64 {
        let s1 = &self.upstream;
        let (xi1, xi2, theta1, hk1) = (s1.xi.value, s2.xi.value, s1.theta.value, s1.hk.value);
        let wake = self.is_wake(i);
        let first_turbulent = self.first_turbulent();
        let target = if i < first_turbulent {
            hk1 + 0.03 * (xi2 - xi1) / theta1
        } else if i == first_turbulent {
            let xi_t = onset.map_or(xi2, |o| o.xi.value);
            hk1 + (0.03 * (xi_t - xi1) - 0.15 * (xi2 - xi_t)) / theta1
        } else if wake {
            // asymptotic wake decay, by backward-Euler Newton steps
            let c = 0.03 * (xi2 - xi1) / theta1;
            (0..3).fold(hk1, |hk, _| hk - (hk + c * (hk - 1.0).powi(3) - hk1) / (1.0 + 3.0 * c * (hk - 1.0).powi(2)))
        } else {
            hk1 - 0.15 * (xi2 - xi1) / theta1
        };
        if wake {
            target.max(1.01)
        } else {
            target.max(hk_max)
        }
    }

    fn apply_limits(&mut self, i: usize) {
        if i >= self.first_turbulent() {
            self.guess.shear = self.guess.shear.clamp(1.0e-7, 0.30);
        }
        let node = self.node(i);
        let msq = self.ctx.freestream.edge_msq(self.guess.ue);
        let floor = hk_floor(self.is_wake(i));
        self.guess.dstar = limit_dstar(self.guess.dstar - node.wake_gap, self.guess.theta, msq, floor) + node.wake_gap;
    }

    /// Replaces a diverged station by values extrapolated from upstream.
    fn extrapolate(&mut self, i: usize) {
        let te = self.te();
        let (prev, node) = (self.node(i - 1), self.node(i));
        if i <= te {
            let growth = (node.xi / prev.xi).sqrt();
            self.guess.theta = prev.theta * growth;
            self.guess.dstar = prev.dstar * growth;
        } else if i == te + 1 {
            (self.guess.shear, self.guess.theta, self.guess.dstar) = self.bl.trailing_edge_values(self.te_gap);
        } else {
            self.guess.theta = prev.theta;
            let length_ratio = (node.xi - prev.xi) / (10.0 * prev.dstar);
            self.guess.dstar = (prev.dstar + self.guess.theta * length_ratio) / (1.0 + length_ratio);
        }
        if i == self.first_turbulent() {
            self.guess.shear = 0.05;
        }
        if i > self.first_turbulent() {
            self.guess.shear = prev.ampl_or_shear;
        }
        self.guess.ue = match self.mode {
            Mode::Mixed => prev.ue,
            Mode::PrescribedUe => {
                let side = &self.bl.sides[self.side];
                if i > 1 && i + 1 < side.len() {
                    0.5 * (prev.ue + side.nodes[i + 1].ue)
                } else {
                    node.ue
                }
            }
        };
    }

    /// Stores the converged station and prepares it as the next upstream station.
    fn store(&mut self, i: usize, onset: Option<Onset>) {
        let first_turbulent = self.first_turbulent();
        let te = self.te();
        let regime = if self.is_wake(i) {
            Regime::Wake
        } else if i >= first_turbulent {
            Regime::Turbulent
        } else {
            Regime::Laminar
        };
        let mut station = Station::new(&self.primary(i), &self.ctx.freestream);
        station.ampl.value = self.guess.ampl;
        station.close(regime, &self.ctx.calibration);

        let g = self.guess;
        let node = &mut self.bl.sides[self.side].nodes[i];
        node.ampl_or_shear = if i < first_turbulent { g.ampl } else { g.shear };
        node.theta = g.theta;
        node.dstar = g.dstar;
        node.ue = g.ue;
        node.mass = g.dstar * g.ue;
        node.outputs = station.outputs();

        if let Some(onset) = onset {
            let transition = &mut self.bl.sides[self.side].transition;
            transition.xi = onset.xi.value;
            transition.forced = onset.forced;
        }
        if onset.is_some() || i == te {
            self.turbulent = true;
        }
        if i == te && self.mode == Mode::PrescribedUe {
            // the wake starts from the combined TE values
            let (_, theta, dstar) = self.bl.trailing_edge_values(self.te_gap);
            (self.guess.theta, self.guess.dstar) = (theta, dstar);
        }
        self.upstream = station.as_upstream();
    }
}
