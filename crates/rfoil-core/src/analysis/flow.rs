//! Solution state and the inviscid / viscous solve sequence.

use super::results::{BoundaryLayerResult, BoundaryLayerSide, IterationInfo, OperatingPoint, SurfaceDistribution};
use super::settings::{FlowNumbers, Scaling, Settings};
use crate::inviscid::{forces, wake, Compressibility, Condition, InviscidSystem, Loads, Panels, UniformFlow};
use crate::numeric::Matrix;
use crate::viscous::layer::Stagnation;
use crate::viscous::march::{self, Mode};
use crate::viscous::newton::{self, Assembly, Constraint, UpdateInputs, UpdateReport};
use crate::viscous::params::{Context, Freestream};
use crate::viscous::{BoundaryLayer, LOWER, UPPER};

/// Free-stream speed; all speeds are normalised by it.
const SPEED: f64 = 1.0;
/// Newton convergence tolerance on the RMS normalised change.
const TOLERANCE: f64 = 1.0e-4;
/// Changes of α or M below this keep the wake and a converged solution valid.
const SAME_CONDITION: f64 = 1.0e-5;

/// What an operating point prescribes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Target {
    /// Angle of attack in degrees.
    Alpha(f64),
    Lift(f64),
}

/// Viscous solution state, kept between operating points as a warm start.
#[derive(Default)]
struct Viscous {
    /// α at which the wake trajectory was traced.
    wake_alpha: Option<f64>,
    /// Source influence matrix `dq/dσ` for the current wake.
    dij: Option<Matrix>,
    bl: Option<BoundaryLayer>,
    /// Whether `bl` holds a marched solution.
    initialized: bool,
    /// `(α, M)` of the last converged solution.
    converged_at: Option<(f64, f64)>,
}

/// Panels, inviscid system and the current flow solution.
pub struct Flow {
    pub panels: Panels,
    system: InviscidSystem,
    alpha: f64,
    loads: Loads,
    cd: f64,
    cd_friction: f64,
    numbers: FlowNumbers,
    /// Surface vorticity (= surface speed) and its α-derivative on the airfoil.
    gamma: Vec<f64>,
    gamma_alpha: Vec<f64>,
    /// Speeds of the α = 0°, 90° unit solutions over airfoil and wake.
    q_unit: [Vec<f64>; 2],
    /// Inviscid speed at the current α, its α-derivative, and the viscous speed.
    q_inviscid: Vec<f64>,
    q_inviscid_alpha: Vec<f64>,
    q_viscous: Vec<f64>,
    viscous: Viscous,
    last: OperatingPoint,
}

impl Flow {
    pub fn new(panels: Panels) -> Self {
        let system = InviscidSystem::new(&panels, SPEED);
        let n = panels.n;
        let q_unit = system.unit_vorticity(n);
        let numbers = FlowNumbers { mach: 0.0, mach_cl: 0.0, reynolds: 0.0, reynolds_cl: 0.0 };
        Self {
            panels,
            system,
            alpha: 0.0,
            loads: Loads::default(),
            cd: 0.0,
            cd_friction: 0.0,
            numbers,
            gamma: vec![0.0; n],
            gamma_alpha: vec![0.0; n],
            q_unit,
            q_inviscid: vec![0.0; n],
            q_inviscid_alpha: vec![0.0; n],
            q_viscous: vec![0.0; n],
            viscous: Viscous::default(),
            last: OperatingPoint::default(),
        }
    }

    pub fn invalidate_convergence(&mut self) {
        self.viscous.converged_at = None;
    }

    pub fn invalidate_wake(&mut self) {
        self.viscous.wake_alpha = None;
        self.viscous.converged_at = None;
    }

    pub fn reset_boundary_layer(&mut self) {
        self.viscous.bl = None;
        self.viscous.initialized = false;
        self.viscous.converged_at = None;
    }

    fn compressibility(&self) -> Compressibility {
        Compressibility::new(self.numbers.mach)
    }

    fn integrate(&self, settings: &Settings) -> Loads {
        let n = self.panels.n;
        let (x, y) = (&self.panels.x[..n], &self.panels.y[..n]);
        let cond = Condition { alpha: self.alpha, speed: SPEED, compressibility: self.compressibility() };
        forces::integrate(x, y, &self.gamma, &self.gamma_alpha, &cond, settings.moment_ref)
    }

    /// Solves one operating point.
    pub fn solve(
        &mut self,
        target: Target,
        settings: &Settings,
        iterations: usize,
        progress: &mut dyn FnMut(&IterationInfo),
    ) -> OperatingPoint {
        match target {
            Target::Alpha(degrees) => {
                self.alpha = degrees.to_radians();
                self.solve_inviscid_alpha(settings);
            }
            Target::Lift(cl) => {
                self.alpha = 0.0;
                self.solve_inviscid_cl(cl, settings);
            }
        }
        self.cd = 0.0;
        self.cd_friction = 0.0;
        self.forget_stale_viscous_state();

        let (converged, iterations, rms) = if settings.is_viscous() {
            let cl_target = match target {
                Target::Alpha(_) => None,
                Target::Lift(cl) => Some(cl),
            };
            self.solve_viscous(cl_target, settings, iterations, progress)
        } else {
            (true, 0, 0.0)
        };
        self.last = self.operating_point(settings, converged, iterations, rms);
        if let Target::Alpha(degrees) = target {
            self.last.alpha = degrees; // report the prescribed value exactly
        }
        self.last
    }

    /// Drops the wake and the convergence mark if α or M changed.
    fn forget_stale_viscous_state(&mut self) {
        let (alpha, mach) = (self.alpha, self.numbers.mach);
        let v = &mut self.viscous;
        if v.wake_alpha.is_some_and(|a| (a - alpha).abs() > SAME_CONDITION) {
            v.wake_alpha = None;
        }
        if v.converged_at.is_some_and(|(a, m)| (a - alpha).abs() > SAME_CONDITION || (m - mach).abs() > SAME_CONDITION)
        {
            v.converged_at = None;
        }
    }

    /// Superimposes the unit solutions for the current α, iterating CL ↔ M(CL)
    /// when the Mach number depends on lift.
    fn solve_inviscid_alpha(&mut self, settings: &Settings) {
        self.set_inviscid_vorticity();
        let mut cl_mach = 1.0;
        self.numbers = settings.flow_numbers(cl_mach);
        self.loads = self.integrate(settings);
        for _ in 0..20 {
            let msq_cl = 2.0 * self.numbers.mach * self.numbers.mach_cl;
            let step = (self.loads.cl - cl_mach) / (1.0 - self.loads.cl_msq * msq_cl);
            let start = cl_mach;
            let mut relaxation = 1.0;
            // halve the step while it would push M(CL) onto its clamp
            for _ in 0..12 {
                cl_mach = start + relaxation * step;
                self.numbers = settings.flow_numbers(cl_mach);
                let mach_fixed = settings.scaling != Scaling::FixedLift;
                if mach_fixed || self.numbers.mach == 0.0 || self.numbers.mach_cl != 0.0 {
                    break;
                }
                relaxation *= 0.5;
            }
            self.loads = self.integrate(settings);
            if step.abs() <= 1.0e-6 {
                break;
            }
        }
        self.numbers = settings.flow_numbers(self.loads.cl);
        self.loads = self.integrate(settings);
    }

    /// Finds the α giving the inviscid lift `cl` at the Mach number of that CL.
    fn solve_inviscid_cl(&mut self, cl: f64, settings: &Settings) {
        self.numbers = settings.flow_numbers(cl);
        self.set_inviscid_vorticity();
        self.loads = self.integrate(settings);
        for _ in 0..20 {
            let dalpha = (cl - self.loads.cl) / self.loads.cl_alpha();
            self.alpha += dalpha;
            self.set_inviscid_vorticity();
            self.loads = self.integrate(settings);
            if dalpha.abs() <= 1.0e-6 {
                break;
            }
        }
    }

    fn set_inviscid_vorticity(&mut self) {
        (self.gamma, self.gamma_alpha) = self.system.vorticity(self.panels.n, self.alpha);
        self.set_inviscid_speeds();
    }

    /// Inviscid speeds at the current α over airfoil and wake.
    fn set_inviscid_speeds(&mut self) {
        let (sa, ca) = self.alpha.sin_cos();
        let [q0, q90] = &self.q_unit;
        self.q_inviscid = q0.iter().zip(q90).map(|(a, b)| ca * a + sa * b).collect();
        self.q_inviscid_alpha = q0.iter().zip(q90).map(|(a, b)| -sa * a + ca * b).collect();
        self.q_viscous.resize(self.q_inviscid.len(), 0.0);
    }

    /// Couples the boundary layer to the inviscid flow by global Newton iteration.
    /// Returns `(converged, iterations, rms)`.
    fn solve_viscous(
        &mut self,
        cl_target: Option<f64>,
        settings: &Settings,
        iterations: usize,
        progress: &mut dyn FnMut(&IterationInfo),
    ) -> (bool, usize, f64) {
        self.prepare_viscous(settings);
        let mut rms = f64::INFINITY;

        for iteration in 1..=iterations {
            // Mach and Re follow the current (or prescribed) lift
            self.numbers = settings.flow_numbers(cl_target.unwrap_or(self.loads.cl));
            let constraint = match cl_target {
                None => Constraint::Alpha {
                    reynolds_cl: self.numbers.reynolds_cl,
                    msq_cl: 2.0 * self.numbers.mach * self.numbers.mach_cl,
                },
                Some(_) => Constraint::Lift,
            };

            let report = self.newton_step(settings, constraint, cl_target);
            rms = report.rms;

            match constraint {
                Constraint::Alpha { .. } => {
                    self.loads.cl += report.parameter_change;
                    self.numbers = settings.flow_numbers(self.loads.cl);
                }
                Constraint::Lift => {
                    self.alpha += report.parameter_change;
                    self.set_inviscid_speeds();
                    let bl = self.viscous.bl.as_mut().expect("boundary layer exists");
                    bl.set_inviscid_speeds(&self.q_inviscid, &self.q_inviscid_alpha);
                }
            }
            self.update_surface_from_layer();
            self.loads = self.integrate(settings);
            self.compute_drag();

            progress(&IterationInfo {
                iteration,
                rms,
                max: report.max,
                relaxation: report.relaxation,
                alpha: self.alpha.to_degrees(),
                cl: self.loads.cl,
                cd: self.cd,
            });
            if rms < TOLERANCE {
                self.viscous.converged_at = Some((self.alpha, self.numbers.mach));
                return (true, iteration, rms);
            }
        }
        (false, iterations, rms)
    }

    /// Wake, source influences and boundary-layer topology for the current α.
    fn prepare_viscous(&mut self, settings: &Settings) {
        if self.viscous.wake_alpha.is_none() {
            let onset = UniformFlow { speed: SPEED, alpha: self.alpha };
            wake::trace(&mut self.panels, &self.gamma, onset, settings.wake_length);
            self.viscous.wake_alpha = Some(self.alpha);
            self.viscous.dij = None;
        }
        self.q_unit = wake::unit_speeds(&self.panels, &self.system.unit_gamma, SPEED);
        self.set_inviscid_speeds();

        let n = self.panels.n;
        if self.viscous.bl.is_none() {
            if self.viscous.initialized {
                self.gamma.copy_from_slice(&self.q_viscous[..n]);
            }
            let stagnation = Stagnation::find(&self.gamma, &self.panels.s[..n]);
            self.viscous.bl = Some(BoundaryLayer::new(&self.panels, stagnation));
        }
        let bl = self.viscous.bl.as_mut().expect("boundary layer exists");
        bl.set_inviscid_speeds(&self.q_inviscid, &self.q_inviscid_alpha);
        if !self.viscous.initialized {
            for node in bl.stations_mut() {
                node.ue = node.ue_inv;
            }
        }
        if self.viscous.converged_at.is_some() {
            // resuming a converged point: restore its viscous lift
            self.update_surface_from_layer();
            self.loads = self.integrate(settings);
        }
        if self.viscous.dij.is_none() {
            self.viscous.dij = Some(self.system.source_influence(&self.panels));
        }
    }

    /// Marches, assembles, solves and applies one Newton step.
    fn newton_step(&mut self, settings: &Settings, constraint: Constraint, cl_target: Option<f64>) -> UpdateReport {
        let freestream = Freestream::new(self.numbers.mach, self.numbers.reynolds, SPEED);
        let min_cl_step = if settings.scaling == Scaling::FixedLift { (-0.9 * self.loads.cl).max(-0.5) } else { -0.5 };
        let cl = cl_target.unwrap_or(self.loads.cl);
        let panels = &self.panels;
        let viscous = &mut self.viscous;
        let bl = viscous.bl.as_mut().expect("boundary layer exists");
        let dij = viscous.dij.as_ref().expect("source influence exists");
        let trips = [settings.xtr_top, settings.xtr_bottom];
        let context = |bl: &BoundaryLayer, side: usize| Context {
            freestream,
            calibration: settings.calibration,
            ncrit: settings.ncrit,
            amplification: settings.amplification,
            xi_trip: bl.trip_xi(panels, side, trips[side]),
        };

        let te_gap = panels.te.normal_gap;
        if !viscous.initialized {
            march::march(bl, te_gap, Mode::PrescribedUe, context);
            viscous.initialized = true;
        }
        march::march(bl, te_gap, Mode::Mixed, context);

        let ctx = [context(bl, UPPER), context(bl, LOWER)];
        let mut sys = newton::assemble(bl, &Assembly { panels, dij, constraint, ctx });
        sys.solve(settings.vaccel, panels.s[panels.n - 1] - panels.s[0]);

        let inputs = UpdateInputs {
            panels,
            dij,
            constraint,
            cl,
            alpha: self.alpha,
            speed: SPEED,
            compressibility: Compressibility::new(self.numbers.mach),
            mach_cl: self.numbers.mach_cl,
            min_cl_step,
            freestream,
        };
        newton::apply(bl, &sys, &inputs)
    }

    /// Surface speeds and vorticity from the boundary-layer edge speeds,
    /// then relocates the stagnation point.
    fn update_surface_from_layer(&mut self) {
        let n = self.panels.n;
        let bl = self.viscous.bl.as_mut().expect("boundary layer exists");
        bl.write_panel_speeds(&mut self.q_viscous);
        self.gamma.copy_from_slice(&self.q_viscous[..n]);
        self.gamma_alpha.copy_from_slice(&self.q_inviscid_alpha[..n]);
        bl.move_stagnation(
            &self.panels,
            &mut self.gamma,
            &mut self.q_viscous,
            &self.q_inviscid,
            &self.q_inviscid_alpha,
        );
    }

    /// Total drag by Squire–Young at the wake end; friction drag by integrating τw.
    fn compute_drag(&mut self) {
        let Some(bl) = &self.viscous.bl else { return };
        let kt = self.compressibility().karman_tsien;
        let wake_end = bl.sides[LOWER].nodes.last().expect("wake stations");
        let speed_ratio = wake_end.ue / SPEED;
        let ue_wake = wake_end.ue * (1.0 - kt) / (1.0 - kt * speed_ratio * speed_ratio);
        let h_wake = wake_end.dstar / wake_end.theta;
        self.cd = 2.0 * wake_end.theta * (ue_wake / SPEED).powf(0.5 * (5.0 + h_wake));

        let (sa, ca) = self.alpha.sin_cos();
        let p = &self.panels;
        self.cd_friction = bl
            .sides
            .iter()
            .flat_map(|side| (2..=side.te).map(move |i| (side.nodes[i - 1], side.nodes[i])))
            .map(|(a, b)| {
                // panel length projected on the free-stream direction
                let dx = (p.x[b.panel] - p.x[a.panel]) * ca + (p.y[b.panel] - p.y[a.panel]) * sa;
                0.5 * (a.outputs.wall_shear + b.outputs.wall_shear) * dx * 2.0 / (SPEED * SPEED)
            })
            .sum();
    }

    fn operating_point(&self, settings: &Settings, converged: bool, iterations: usize, rms: f64) -> OperatingPoint {
        let viscous = settings.is_viscous();
        let surface = self.surface_for(viscous);
        let cp_min = surface.cp_viscous[..surface.airfoil_nodes].iter().copied().fold(f64::INFINITY, f64::min);
        let (xtr_top, xtr_bottom) = match (&self.viscous.bl, viscous) {
            (Some(bl), true) => (bl.sides[UPPER].transition.x_over_c, bl.sides[LOWER].transition.x_over_c),
            _ => (1.0, 1.0),
        };
        OperatingPoint {
            alpha: self.alpha.to_degrees(),
            cl: self.loads.cl,
            cm: self.loads.cm,
            cd: self.cd,
            cdf: self.cd_friction,
            cdp: self.cd - self.cd_friction,
            cp_min,
            xtr_top,
            xtr_bottom,
            mach: self.numbers.mach,
            reynolds: self.numbers.reynolds,
            viscous,
            converged,
            iterations,
            rms_change: rms,
        }
    }

    pub fn surface(&self) -> SurfaceDistribution {
        self.surface_for(self.last.viscous)
    }

    fn surface_for(&self, viscous: bool) -> SurfaceDistribution {
        let compressibility = self.compressibility();
        let n = self.panels.n;
        let viscous = viscous && self.viscous.bl.is_some();
        let count = if viscous { self.panels.len() } else { n };
        let cp = |q: &[f64]| q[..count].iter().map(|&v| compressibility.cp_from_speed(v, SPEED)).collect::<Vec<_>>();
        let cp_inviscid = cp(&self.q_inviscid);
        SurfaceDistribution {
            x: self.panels.x[..count].to_vec(),
            y: self.panels.y[..count].to_vec(),
            cp_viscous: if viscous { cp(&self.q_viscous) } else { cp_inviscid.clone() },
            cp_inviscid,
            airfoil_nodes: n,
        }
    }

    pub fn boundary_layer(&self) -> Option<BoundaryLayerResult> {
        let bl = self.viscous.bl.as_ref()?;
        if !self.viscous.initialized {
            return None;
        }
        let p = &self.panels;
        let side = |index: usize| {
            let s = &bl.sides[index];
            let mut out = BoundaryLayerSide { te: s.te - 1, ..Default::default() };
            for (i, node) in s.stations() {
                let laminar = s.is_laminar(i);
                out.x.push(p.x[node.panel]);
                out.y.push(p.y[node.panel]);
                out.xi.push(node.xi);
                out.ue.push(node.ue);
                out.dstar.push(node.dstar);
                out.theta.push(node.theta);
                out.h.push((node.dstar - node.wake_gap) / node.theta);
                out.cf.push(node.outputs.wall_shear / (0.5 * SPEED * SPEED));
                out.n.push(if laminar { node.ampl_or_shear } else { 0.0 });
                out.ctau.push(if laminar { 0.0 } else { node.ampl_or_shear });
            }
            out
        };
        Some(BoundaryLayerResult { top: side(UPPER), bottom: side(LOWER) })
    }
}
