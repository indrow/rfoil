//! Kármán–Tsien compressibility correction and pressure integration.

/// Ratio of specific heats (air).
pub const GAMMA: f64 = 1.4;

/// Kármán–Tsien parameters for a free-stream Mach number.
#[derive(Clone, Copy, Debug)]
pub struct Compressibility {
    pub mach: f64,
    /// `β = √(1 − M²)` and `dβ/dM²`.
    beta: f64,
    beta_msq: f64,
    /// Denominator factor `M²/(2(1 + β))` and its M²-derivative.
    cp_factor: f64,
    cp_factor_msq: f64,
    /// Kármán–Tsien λ = `M²/(1 + β)²` and `dλ/dM²`.
    pub karman_tsien: f64,
    pub karman_tsien_msq: f64,
}

impl Compressibility {
    pub fn new(mach: f64) -> Self {
        let msq = mach * mach;
        let beta = (1.0 - msq).sqrt();
        let beta_msq = -0.5 / beta;
        let cp_factor = 0.5 * mach * mach / (1.0 + beta);
        let karman_tsien = msq / (1.0 + beta).powi(2);
        Self {
            mach,
            beta,
            beta_msq,
            cp_factor,
            cp_factor_msq: 0.5 / (1.0 + beta) - cp_factor / (1.0 + beta) * beta_msq,
            karman_tsien,
            karman_tsien_msq: 1.0 / (1.0 + beta).powi(2) - 2.0 * karman_tsien / (1.0 + beta) * beta_msq,
        }
    }

    /// Compressible Cp from the incompressible Cp.
    pub fn cp(&self, cp_incompressible: f64) -> f64 {
        cp_incompressible / (self.beta + self.cp_factor * cp_incompressible)
    }

    /// Cp from a local speed `q`.
    pub fn cp_from_speed(&self, q: f64, speed: f64) -> f64 {
        self.cp(1.0 - (q / speed).powi(2))
    }

    /// Cp at sonic conditions, or `None` in incompressible flow.
    pub fn sonic_cp(&self) -> Option<f64> {
        let msq = self.mach * self.mach;
        (self.mach > 0.0).then(|| {
            let g1 = GAMMA - 1.0;
            2.0 / (GAMMA * msq) * (((1.0 + 0.5 * g1 * msq) / (1.0 + 0.5 * g1)).powf(GAMMA / g1) - 1.0)
        })
    }

    /// Compressible Cp from a speed, with `dCp/dM²` and `dCp/dq`.
    fn cp_linearized(&self, q: f64, speed: f64) -> (f64, f64, f64) {
        let cp_incompressible = 1.0 - (q / speed).powi(2);
        let den = self.beta + self.cp_factor * cp_incompressible;
        let cp = cp_incompressible / den;
        let cp_msq = -cp / den * (self.beta_msq + self.cp_factor_msq * cp_incompressible);
        let cp_q = (1.0 - self.cp_factor * cp) / den * (-2.0 * q / (speed * speed));
        (cp, cp_msq, cp_q)
    }
}

/// Integrated pressure forces.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Loads {
    pub cl: f64,
    pub cm: f64,
    /// Pressure drag (meaningful only as a check in inviscid flow).
    pub cdp: f64,
    /// `dCL` per unit change of the speed perturbation `q_a` passed to [`integrate`].
    pub cl_speed: f64,
    /// `dCL/dα` from rotating the force axes at fixed surface pressures.
    pub cl_geometry: f64,
    pub cl_msq: f64,
}

impl Loads {
    /// Total `dCL/dα` when `q_a` is the α-derivative of the surface speed.
    pub fn cl_alpha(&self) -> f64 {
        self.cl_speed + self.cl_geometry
    }
}

/// Free-stream condition seen by the airfoil.
#[derive(Clone, Copy, Debug)]
pub struct Condition {
    pub alpha: f64,
    pub speed: f64,
    pub compressibility: Compressibility,
}

/// Pressure at one contour node with its sensitivities.
#[derive(Clone, Copy)]
struct NodePressure {
    cp: f64,
    cp_msq: f64,
    /// Change of Cp per unit speed perturbation.
    cp_perturbation: f64,
}

/// Integrates surface pressures over the closed contour `(x, y)` with surface
/// speeds `q` and a speed perturbation `q_a`.  Moments are taken about `moment_ref`.
pub fn integrate(x: &[f64], y: &[f64], q: &[f64], q_a: &[f64], cond: &Condition, moment_ref: (f64, f64)) -> Loads {
    let n = x.len();
    let Condition { alpha, speed, ref compressibility } = *cond;
    let (x_ref, y_ref) = moment_ref;
    let (sa, ca) = alpha.sin_cos();
    let mut loads = Loads::default();

    let pressure = |i: usize| {
        let (cp, cp_msq, cp_q) = compressibility.cp_linearized(q[i], speed);
        NodePressure { cp, cp_msq, cp_perturbation: cp_q * q_a[i] }
    };
    let mut p1 = pressure(0);
    for i in 0..n {
        let next = (i + 1) % n;
        let p2 = pressure(next);

        // panel extent in wind axes, and its α-derivative
        let (panel_dx, panel_dy) = (x[next] - x[i], y[next] - y[i]);
        let dx = panel_dx * ca + panel_dy * sa;
        let dy = panel_dy * ca - panel_dx * sa;
        let dx_alpha = -panel_dx * sa + panel_dy * ca;

        // panel midpoint relative to the moment reference, in wind axes
        let (xm, ym) = (0.5 * (x[next] + x[i]) - x_ref, 0.5 * (y[next] + y[i]) - y_ref);
        let arm_x = xm * ca + ym * sa;
        let arm_y = ym * ca - xm * sa;

        let cp_mean = 0.5 * (p2.cp + p1.cp);
        let cp_jump = p2.cp - p1.cp;
        loads.cl += dx * cp_mean;
        loads.cdp -= dy * cp_mean;
        loads.cm -= dx * (cp_mean * arm_x + cp_jump * dx / 12.0) + dy * (cp_mean * arm_y + cp_jump * dy / 12.0);
        loads.cl_speed += dx * 0.5 * (p2.cp_perturbation + p1.cp_perturbation);
        loads.cl_geometry += cp_mean * dx_alpha;
        loads.cl_msq += dx * 0.5 * (p2.cp_msq + p1.cp_msq);

        p1 = p2;
    }
    loads
}
