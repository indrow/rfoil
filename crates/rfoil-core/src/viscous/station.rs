//! State of the boundary layer at one streamwise station.

use super::closure::{self, Fit};
use super::dual::{as_upstream, floor_value, slot, var, D};
use super::params::{Calibration, Freestream, GAMMA_MINUS_1, SUTHERLAND_RATIO};

/// Flow regime governing which closure relations apply.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Regime {
    Laminar,
    Turbulent,
    Wake,
}

/// Primary variables of a station as stored between Newton iterations.
#[derive(Clone, Copy, Debug, Default)]
pub struct Primary {
    /// Arc length from the stagnation point.
    pub xi: f64,
    /// Amplification exponent `N` (meaningful while laminar).
    pub ampl: f64,
    /// Shear-stress coefficient `√Cτ` (meaningful while turbulent).
    pub shear: f64,
    pub theta: f64,
    /// Displacement thickness including any wake gap.
    pub dstar: f64,
    /// Wake gap (dead-air thickness behind a blunt TE).
    pub wake_gap: f64,
    /// Incompressible edge speed.
    pub ue: f64,
}

/// A station with all derived quantities as duals.
#[derive(Clone, Copy, Debug, Default)]
pub struct Station {
    // primary variables
    pub xi: D,
    /// Compressible edge speed (partials taken w.r.t. the incompressible speed).
    pub ue: D,
    pub theta: D,
    /// Displacement thickness excluding the wake gap.
    pub dstar: D,
    /// Shear-stress coefficient `√Cτ`.
    pub shear: D,
    /// Amplification exponent `N`.
    pub ampl: D,
    pub wake_gap: f64,

    // kinematics
    /// Edge Mach number squared.
    pub msq: D,
    /// Edge density.
    pub rho: D,
    /// Shape parameter `δ*/θ`.
    pub h: D,
    /// Kinematic shape parameter.
    pub hk: D,
    /// Momentum-thickness Reynolds number.
    pub re_theta: D,

    // closure relations
    /// Density-thickness shape parameter `H**`.
    pub h_density: D,
    /// Kinetic-energy shape parameter `H*`.
    pub h_energy: D,
    /// Normalised wall slip velocity `Us`.
    pub slip: D,
    /// Equilibrium shear coefficient `√Cτ_eq`.
    pub shear_eq: D,
    pub cf: D,
    /// Dissipation coefficient `2·CD/H*`.
    pub dissipation: D,
    /// Boundary-layer thickness `δ`.
    pub delta: D,
}

impl Station {
    /// Seeds `p` as the downstream station of an interval and evaluates its
    /// kinematics.  Closures still need [`Station::close`].
    pub fn new(p: &Primary, fs: &Freestream) -> Self {
        let downstream = |v| slot::downstream(v);
        // Kármán–Tsien conversion of the incompressible edge speed
        let ue_incompressible = D::var(p.ue, downstream(var::UE));
        let lambda = fs.karman_tsien;
        let q = ue_incompressible / fs.speed;
        let ue = ue_incompressible * (1.0 - lambda) / (1.0 - lambda * q * q);
        Self {
            xi: D::var(p.xi, downstream(var::XI)),
            ue,
            theta: D::var(p.theta, downstream(var::THETA)),
            dstar: D::var(p.dstar - p.wake_gap, downstream(var::DSTAR)),
            shear: D::var(p.shear, downstream(var::AMPL_OR_SHEAR)),
            ampl: D::var(p.ampl, downstream(var::AMPL_OR_SHEAR)),
            wake_gap: p.wake_gap,
            ..Self::default()
        }
        .with_kinematics(fs)
    }

    /// Evaluates the kinematic quantities from the primary variables.
    pub fn with_kinematics(mut self, fs: &Freestream) -> Self {
        let ue_sq_h = self.ue * self.ue * fs.inv_stagnation_enthalpy;
        self.msq = ue_sq_h / (GAMMA_MINUS_1 * (1.0 - 0.5 * ue_sq_h));
        let temperature_ratio = 1.0 + 0.5 * GAMMA_MINUS_1 * self.msq;
        self.rho = fs.stagnation_density * temperature_ratio.powf(-1.0 / GAMMA_MINUS_1);
        self.h = self.dstar / self.theta;

        // kinematic viscosity from Sutherland's law in terms of the edge enthalpy
        let enthalpy_ratio = 1.0 - 0.5 * ue_sq_h;
        let viscosity =
            enthalpy_ratio.powf(1.5) * (1.0 + SUTHERLAND_RATIO) / (enthalpy_ratio + SUTHERLAND_RATIO) / fs.reynolds;

        let (hk, hk_h, hk_msq) = closure::kinematic_shape(self.h.value, self.msq.value);
        self.hk = D::chain(hk, &[(hk_h, self.h), (hk_msq, self.msq)]);
        self.re_theta = self.rho * self.ue * self.theta / viscosity;
        self
    }

    /// Re-expresses the station as the upstream station of the next interval.
    pub fn as_upstream(self) -> Self {
        let f = as_upstream;
        Self {
            xi: f(self.xi),
            ue: f(self.ue),
            theta: f(self.theta),
            dstar: f(self.dstar),
            shear: f(self.shear),
            ampl: f(self.ampl),
            msq: f(self.msq),
            rho: f(self.rho),
            h: f(self.h),
            hk: f(self.hk),
            re_theta: f(self.re_theta),
            h_density: f(self.h_density),
            h_energy: f(self.h_energy),
            slip: f(self.slip),
            shear_eq: f(self.shear_eq),
            cf: f(self.cf),
            dissipation: f(self.dissipation),
            delta: f(self.delta),
            ..self
        }
    }

    fn lift(&self, fit: Fit) -> D {
        fit.lift(self.hk, self.re_theta, self.msq)
    }

    /// Evaluates the closure relations for the given regime.
    pub fn close(&mut self, regime: Regime, cal: &Calibration) {
        let wake = regime == Regime::Wake;
        self.hk = floor_value(self.hk, if wake { 1.00005 } else { 1.05 });
        let (hk, re_theta, msq) = (self.hk.value, self.re_theta.value, self.msq.value);

        self.h_density = self.lift(closure::density_shape(hk, msq));
        self.h_energy = self.lift(match regime {
            Regime::Laminar => closure::energy_shape_laminar(hk),
            _ => closure::energy_shape_turbulent(hk, re_theta, msq),
        });

        // normalised slip velocity, clamped for numerical safety
        let slip = 0.5 * self.h_energy * (1.0 - (self.hk - 1.0) / (cal.g_beta_b * self.h));
        let slip_max = if wake { 0.99995 } else { 0.95 };
        self.slip = if slip.value > slip_max { D::constant(if wake { 0.99995 } else { 0.98 }) } else { slip };

        self.shear_eq = self.equilibrium_shear(regime, cal);
        self.cf = self.skin_friction(regime, cal);
        self.dissipation = self.dissipation_coefficient(regime, cal);

        // BL thickness from Green's correlation, limited to 12θ
        let delta = (3.15 + 1.72 / (self.hk - 1.0)) * self.theta + self.dstar;
        self.delta = if delta.value > 12.0 * self.theta.value { 12.0 * self.theta } else { delta };
    }

    /// Equilibrium shear coefficient `√Cτ_eq`.
    fn equilibrium_shear(&self, regime: Regime, cal: &Calibration) -> D {
        let (hk, h_energy, slip, h) = (self.hk, self.h_energy, self.slip, self.h);
        let hk_excess = match regime {
            Regime::Turbulent => low_reynolds_hk_excess(hk, self.re_theta, cal),
            _ => hk - 1.0,
        };
        (cal.ctau_eq_scale() * h_energy * (hk - 1.0) * hk_excess * hk_excess / ((1.0 - slip) * h * hk * hk)).sqrt()
    }

    fn skin_friction(&self, regime: Regime, cal: &Calibration) -> D {
        let (hk, re_theta, msq) = (self.hk.value, self.re_theta.value, self.msq.value);
        match regime {
            Regime::Wake => D::constant(0.0),
            Regime::Laminar => self.lift(closure::skin_friction_laminar(hk, re_theta)),
            Regime::Turbulent => {
                // laminar Cf only wins at unreasonably small Rθ
                let turbulent = closure::skin_friction_turbulent(hk, re_theta, msq, cal.cf_scale);
                self.lift(turbulent.max(closure::skin_friction_laminar(hk, re_theta)))
            }
        }
    }

    /// Dissipation coefficient `2·CD/H*`.
    fn dissipation_coefficient(&self, regime: Regime, cal: &Calibration) -> D {
        let (hk, re_theta, msq) = (self.hk.value, self.re_theta.value, self.msq.value);
        if regime == Regime::Laminar {
            return self.lift(closure::dissipation_laminar(hk, re_theta));
        }
        let (h_energy, slip) = (self.h_energy, self.slip);

        let mut di = if regime == Regime::Turbulent {
            // wall contribution, reduced for very low Hk where the wake layer vanishes
            let cf = self.lift(closure::skin_friction_turbulent(hk, re_theta, msq, cal.cf_scale));
            let hk_min = 1.0 + 2.1 / self.re_theta.ln();
            let fade = 0.5 + 0.5 * ((self.hk - 1.0) / (hk_min - 1.0)).tanh();
            cf * slip / h_energy * fade
        } else {
            D::constant(0.0)
        };

        // outer layer, plus the laminar stress contribution
        di += self.shear * self.shear * (0.995 - slip) * 2.0 / h_energy;
        di += 0.15 * (0.995 - slip) * (0.995 - slip) / self.re_theta * 2.0 / h_energy;

        match regime {
            Regime::Turbulent => {
                let laminar = self.lift(closure::dissipation_laminar(hk, re_theta));
                larger(laminar, di)
            }
            _ => {
                // two wake halves
                let laminar = self.lift(closure::dissipation_laminar_wake(hk, re_theta));
                2.0 * larger(laminar, di)
            }
        }
    }

    /// Quantities stored for output.
    pub fn outputs(&self) -> StationOutputs {
        let (rho, ue) = (self.rho.value, self.ue.value);
        StationOutputs {
            wall_shear: 0.5 * rho * ue * ue * self.cf.value,
            dissipation: 0.5 * rho * ue * ue * ue * self.dissipation.value * self.h_energy.value,
            shear_eq: self.shear_eq.value,
            delta: self.delta.value,
            energy_thickness: self.h_energy.value * self.theta.value,
        }
    }
}

/// `Hk − 1` reduced at low Reynolds numbers (the `C/(H·Rθ)` term of the
/// G–β locus), floored at 0.01.
pub(super) fn low_reynolds_hk_excess(hk: D, re_theta: D, cal: &Calibration) -> D {
    let excess = hk - 1.0 - cal.g_beta_c / re_theta;
    if excess.value < 0.01 {
        D::constant(0.01)
    } else {
        excess
    }
}

/// `a` if its value exceeds that of `b`, else `b`.
fn larger(a: D, b: D) -> D {
    if a.value > b.value {
        a
    } else {
        b
    }
}

/// Derived quantities kept for plotting and drag integration.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StationOutputs {
    /// Wall shear stress, normalised by `ρ∞·V∞²`.
    pub wall_shear: f64,
    pub dissipation: f64,
    pub shear_eq: f64,
    pub delta: f64,
    pub energy_thickness: f64,
}
