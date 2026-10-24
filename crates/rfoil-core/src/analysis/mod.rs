//! High-level airfoil analysis: the public entry point of the library.
//!
//! ```no_run
//! use rfoil_core::{Airfoil, Analysis, Settings};
//!
//! let mut analysis = Analysis::new(Airfoil::naca("2412")?)?;
//! analysis.set_settings(Settings { reynolds: 1.0e6, ..Settings::default() });
//! let point = analysis.solve_alpha(4.0);
//! println!("CL = {:.4}  CD = {:.5}", point.cl, point.cd);
//! # Ok::<(), rfoil_core::Error>(())
//! ```

mod flow;
mod results;
mod settings;

pub use results::{BoundaryLayerResult, BoundaryLayerSide, IterationInfo, OperatingPoint, SurfaceDistribution};
pub use settings::{FlowNumbers, Scaling, Settings};

use flow::{Flow, Target};

use crate::error::Result;
use crate::geometry::{panel, Airfoil, PanelingOptions, SectionShape};
use crate::inviscid::Panels;

/// Consecutive non-converged points after which a sweep stops.
const SWEEP_FAILURE_LIMIT: usize = 4;

/// An airfoil together with its paneling, settings and current solution.
pub struct Analysis {
    airfoil: Airfoil,
    paneling: PanelingOptions,
    settings: Settings,
    flow: Flow,
}

impl Analysis {
    /// Panels `airfoil` with default options and sets up the inviscid system.
    pub fn new(airfoil: Airfoil) -> Result<Self> {
        Self::with_paneling(airfoil, PanelingOptions::default())
    }

    pub fn with_paneling(airfoil: Airfoil, paneling: PanelingOptions) -> Result<Self> {
        let flow = Flow::new(paneled(&airfoil, &paneling)?);
        Ok(Self { airfoil, paneling, settings: Settings::default(), flow })
    }

    pub fn airfoil(&self) -> &Airfoil {
        &self.airfoil
    }

    pub fn paneling(&self) -> &PanelingOptions {
        &self.paneling
    }

    pub fn settings(&self) -> &Settings {
        &self.settings
    }

    /// The paneled airfoil (and wake, once a viscous solution exists).
    pub fn panels(&self) -> &Panels {
        &self.flow.panels
    }

    /// Thickness and camber of the paneled airfoil.
    pub fn shape(&self) -> SectionShape {
        self.flow.panels.curve().thickness_and_camber()
    }

    /// Replaces the airfoil; all solution state is discarded.
    pub fn set_airfoil(&mut self, airfoil: Airfoil) -> Result<()> {
        self.flow = Flow::new(paneled(&airfoil, &self.paneling)?);
        self.airfoil = airfoil;
        Ok(())
    }

    /// Re-panels the airfoil; all solution state is discarded.
    pub fn set_paneling(&mut self, paneling: PanelingOptions) -> Result<()> {
        self.flow = Flow::new(paneled(&self.airfoil, &paneling)?);
        self.paneling = paneling;
        Ok(())
    }

    /// Changes the settings.  An existing boundary layer is kept as the
    /// initial guess.
    pub fn set_settings(&mut self, settings: Settings) {
        if self.settings.changes_viscous_problem(&settings) {
            self.flow.invalidate_convergence();
        }
        if self.settings.wake_length != settings.wake_length {
            self.flow.invalidate_wake();
        }
        self.settings = settings;
    }

    /// Discards the boundary layer so that the next point starts afresh.
    pub fn reset_boundary_layer(&mut self) {
        self.flow.reset_boundary_layer();
    }

    /// Solves for a prescribed angle of attack in degrees.
    pub fn solve_alpha(&mut self, alpha_deg: f64) -> OperatingPoint {
        self.solve_alpha_with(alpha_deg, &mut |_| {})
    }

    /// Like [`Analysis::solve_alpha`], reporting each Newton iteration.
    pub fn solve_alpha_with(&mut self, alpha_deg: f64, progress: &mut dyn FnMut(&IterationInfo)) -> OperatingPoint {
        self.flow.solve(Target::Alpha(alpha_deg), &self.settings, self.settings.max_iterations, progress)
    }

    /// Solves for a prescribed lift coefficient (viscous CL when viscous).
    pub fn solve_cl(&mut self, cl: f64) -> OperatingPoint {
        self.solve_cl_with(cl, &mut |_| {})
    }

    pub fn solve_cl_with(&mut self, cl: f64, progress: &mut dyn FnMut(&IterationInfo)) -> OperatingPoint {
        self.flow.solve(Target::Lift(cl), &self.settings, self.settings.max_iterations, progress)
    }

    /// Sweeps α from `start` to `end` (degrees, inclusive) in steps of `step`.
    /// `on_point` sees every point and may return `false` to stop early.
    pub fn sweep_alpha(
        &mut self,
        start: f64,
        end: f64,
        step: f64,
        on_point: impl FnMut(&OperatingPoint) -> bool,
    ) -> Vec<OperatingPoint> {
        let targets = sequence(start, end, step).map(Target::Alpha);
        self.sweep(targets, on_point)
    }

    /// Sweeps CL from `start` to `end` (inclusive) in steps of `step`.
    pub fn sweep_cl(
        &mut self,
        start: f64,
        end: f64,
        step: f64,
        on_point: impl FnMut(&OperatingPoint) -> bool,
    ) -> Vec<OperatingPoint> {
        self.sweep(sequence(start, end, step).map(Target::Lift), on_point)
    }

    fn sweep(
        &mut self,
        targets: impl Iterator<Item = Target>,
        mut on_point: impl FnMut(&OperatingPoint) -> bool,
    ) -> Vec<OperatingPoint> {
        // sequences allow a few more iterations per point
        let iterations = self.settings.max_iterations + 5;
        let mut points = Vec::new();
        let mut failures = 0;
        for target in targets {
            let point = self.flow.solve(target, &self.settings, iterations, &mut |_| {});
            failures = if point.converged { 0 } else { failures + 1 };
            points.push(point);
            if !on_point(&point) || failures >= SWEEP_FAILURE_LIMIT {
                break;
            }
        }
        points
    }

    /// Pressure distribution of the current solution.
    pub fn surface(&self) -> SurfaceDistribution {
        self.flow.surface()
    }

    /// Boundary-layer distributions, if a viscous solution exists.
    pub fn boundary_layer(&self) -> Option<BoundaryLayerResult> {
        self.flow.boundary_layer()
    }
}

fn paneled(airfoil: &Airfoil, opts: &PanelingOptions) -> Result<Panels> {
    if opts.nodes < 20 {
        return Err(crate::Error::InvalidParameter(format!("at least 20 panel nodes needed, got {}", opts.nodes)));
    }
    Ok(Panels::from_curve(&panel(&airfoil.curve(), opts)))
}

/// `start, start + step, ...` up to and including `end` (direction follows `end`).
fn sequence(start: f64, end: f64, step: f64) -> impl Iterator<Item = f64> {
    let step = step.abs().copysign(end - start);
    let count = if step == 0.0 { 1 } else { ((end - start) / step + 0.5) as usize + 1 };
    (0..count).map(move |i| start + step * i as f64)
}
