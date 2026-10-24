//! Empirical closure relations of the integral boundary-layer model.
//!
//! Each correlation returns its value together with the partial derivatives
//! with respect to its arguments, so that callers can lift them into dual
//! numbers with [`Fit::lift`].  Arguments are the kinematic shape parameter
//! `hk`, the momentum-thickness Reynolds number `re_theta` and the edge Mach
//! number squared `msq`.

use std::f64::consts::LN_10;

use super::dual::D;
use super::params::GAMMA_MINUS_1;

/// A correlation value with its partials `∂/∂Hk`, `∂/∂Rθ`, `∂/∂M²`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Fit {
    pub value: f64,
    pub d_hk: f64,
    pub d_re_theta: f64,
    pub d_msq: f64,
}

impl Fit {
    /// Applies the chain rule through dual-valued arguments.
    pub fn lift(self, hk: D, re_theta: D, msq: D) -> D {
        D::chain(self.value, &[(self.d_hk, hk), (self.d_re_theta, re_theta), (self.d_msq, msq)])
    }

    /// The larger of two fits by value.
    pub fn max(self, other: Fit) -> Fit {
        if other.value > self.value {
            other
        } else {
            self
        }
    }
}

/// Kinematic shape parameter from `H` and `M²` (Whitfield): `(Hk, ∂Hk/∂H, ∂Hk/∂M²)`.
pub fn kinematic_shape(h: f64, msq: f64) -> (f64, f64, f64) {
    let den = 1.0 + 0.113 * msq;
    let hk = (h - 0.29 * msq) / den;
    (hk, 1.0 / den, (-0.29 - 0.113 * hk) / den)
}

/// Density-thickness shape parameter `H**` (Whitfield).
pub fn density_shape(hk: f64, msq: f64) -> Fit {
    let base = 0.064 / (hk - 0.8) + 0.251;
    Fit { value: msq * base, d_hk: msq * (-0.064 / (hk - 0.8).powi(2)), d_re_theta: 0.0, d_msq: base }
}

/// Laminar kinetic-energy shape parameter `H*` (Falkner–Skan).
pub fn energy_shape_laminar(hk: f64) -> Fit {
    let t = hk - 4.35;
    let (value, d_hk) = if hk < 4.35 {
        let hp1 = hk + 1.0;
        (
            0.0111 * t * t / hp1 - 0.0278 * t.powi(3) / hp1 + 1.528 - 0.0002 * (t * hk).powi(2),
            0.0111 * (2.0 * t - t * t / hp1) / hp1
                - 0.0278 * (3.0 * t * t - t.powi(3) / hp1) / hp1
                - 0.0002 * 2.0 * t * hk * (t + hk),
        )
    } else {
        (0.015 * t * t / hk + 1.528, 0.015 * 2.0 * t / hk - 0.015 * t * t / (hk * hk))
    };
    Fit { value, d_hk, ..Fit::default() }
}

/// Turbulent kinetic-energy shape parameter `H*` (arctan(y⁺) + Schlichting profiles).
pub fn energy_shape_turbulent(hk: f64, re_theta: f64, msq: f64) -> Fit {
    const HS_MIN: f64 = 1.5;
    const DHS_INF: f64 = 0.015;

    // Hk of minimum H*, and Rθ floored at 200, each with its Rθ-derivative
    let (hk_min, hk_min_rt) =
        if re_theta > 400.0 { (3.0 + 400.0 / re_theta, -400.0 / (re_theta * re_theta)) } else { (4.0, 0.0) };
    let (rt, rt_rt) = if re_theta > 200.0 { (re_theta, 1.0) } else { (200.0, 0.0) };

    let (hs, hs_hk, hs_rt) = if hk < hk_min {
        // attached branch
        let hr = (hk_min - hk) / (hk_min - 1.0);
        let hr_hk = -1.0 / (hk_min - 1.0);
        let hr_rt = (1.0 - hr) / (hk_min - 1.0) * hk_min_rt;
        let a = 2.0 - HS_MIN - 4.0 / rt;
        let f = 1.5 / (hk + 0.5);
        (
            a * hr * hr * f + HS_MIN + 4.0 / rt,
            -a * hr * hr * 1.5 / (hk + 0.5).powi(2) + a * hr * 2.0 * f * hr_hk,
            a * hr * 2.0 * f * hr_rt + (hr * hr * f - 1.0) * 4.0 / (rt * rt) * rt_rt,
        )
    } else {
        // separated branch
        let ln_rt = rt.ln();
        let excess = hk - hk_min;
        let shifted = hk - hk_min + 4.0 / ln_rt;
        let g = 0.007 * ln_rt / (shifted * shifted) + DHS_INF / hk;
        let g_hk = -0.014 * ln_rt / shifted.powi(3) - DHS_INF / (hk * hk);
        let g_rt = -0.014 * ln_rt / shifted.powi(3) * (-hk_min_rt - 4.0 / (ln_rt * ln_rt) / rt * rt_rt)
            + 0.007 / (shifted * shifted) / rt * rt_rt;
        (
            excess * excess * g + HS_MIN + 4.0 / rt,
            excess * 2.0 * g + excess * excess * g_hk,
            excess * excess * g_rt - 4.0 / (rt * rt) * rt_rt + excess * 2.0 * g * (-hk_min_rt),
        )
    };

    // Whitfield's compressibility correction
    let fm = 1.0 + 0.014 * msq;
    let value = (hs + 0.028 * msq) / fm;
    Fit { value, d_hk: hs_hk / fm, d_re_theta: hs_rt / fm, d_msq: 0.028 / fm - 0.014 * value / fm }
}

/// Laminar skin friction `Cf` (Falkner–Skan).
pub fn skin_friction_laminar(hk: f64, re_theta: f64) -> Fit {
    let (num, num_hk) = if hk < 5.5 {
        let tmp = (5.5 - hk).powi(3) / (hk + 1.0);
        (0.0727 * tmp - 0.07, -0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0))
    } else {
        let tmp = 1.0 - 1.0 / (hk - 4.5);
        (0.015 * tmp * tmp - 0.07, 0.015 * tmp * 2.0 / (hk - 4.5).powi(2))
    };
    Fit { value: num / re_theta, d_hk: num_hk / re_theta, d_re_theta: -num / (re_theta * re_theta), d_msq: 0.0 }
}

/// Turbulent skin friction `Cf` (Coles), scaled by `cf_scale`.
pub fn skin_friction_turbulent(hk: f64, re_theta: f64, msq: f64, cf_scale: f64) -> Fit {
    let fc = (1.0 + 0.5 * GAMMA_MINUS_1 * msq).sqrt();
    let ln_rt = (re_theta / fc).ln().max(3.0);
    let exponent = -1.74 - 0.31 * hk;
    let arg = (-1.33 * hk).max(-20.0);
    let thk = (4.0 - hk / 0.875).tanh();

    let cf0 = cf_scale * 0.3 * arg.exp() * (ln_rt / LN_10).powf(exponent);
    let value = (cf0 + 1.1e-4 * (thk - 1.0)) / fc;
    Fit {
        value,
        d_hk: (-1.33 * cf0 - 0.31 * (ln_rt / LN_10).ln() * cf0 - 1.1e-4 * (1.0 - thk * thk) / 0.875) / fc,
        d_re_theta: exponent * cf0 / (fc * ln_rt) / re_theta,
        d_msq: exponent * cf0 / (fc * ln_rt) * (-0.25 * GAMMA_MINUS_1 / (fc * fc))
            - 0.25 * GAMMA_MINUS_1 * value / (fc * fc),
    }
}

/// Laminar dissipation coefficient `2·CD/H*` (Falkner–Skan).
pub fn dissipation_laminar(hk: f64, re_theta: f64) -> Fit {
    let (num, num_hk) = if hk < 4.0 {
        (0.00205 * (4.0 - hk).powf(5.5) + 0.207, -0.00205 * 5.5 * (4.0 - hk).powf(4.5))
    } else {
        let hkb = hk - 4.0;
        let den = 1.0 + 0.02 * hkb * hkb;
        (-0.0016 * hkb * hkb / den + 0.207, -0.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb * hkb / (den * den)))
    };
    Fit { value: num / re_theta, d_hk: num_hk / re_theta, d_re_theta: -num / (re_theta * re_theta), d_msq: 0.0 }
}

/// Laminar wake dissipation coefficient `2·CD/H*`.
pub fn dissipation_laminar_wake(hk: f64, re_theta: f64) -> Fit {
    let hs = energy_shape_laminar(hk);
    let rcd = 1.10 * (1.0 - 1.0 / hk).powi(2) / hk;
    let rcd_hk = 1.10 * (1.0 - 1.0 / hk) * 2.0 / hk.powi(3) - rcd / hk;
    let value = 2.0 * rcd / (hs.value * re_theta);
    Fit {
        value,
        d_hk: 2.0 * rcd_hk / (hs.value * re_theta) - value / hs.value * hs.d_hk,
        d_re_theta: -value / re_theta,
        d_msq: 0.0,
    }
}

/// Which amplification-rate model to use for eᴺ transition prediction.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Amplification {
    /// Falkner–Skan envelope method.
    #[default]
    Envelope,
    /// Envelope with Orr–Sommerfeld maximum rates for separated profiles.
    Modified,
}

/// Spatial amplification rate `dN/dξ` with its partials.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AmplificationRate {
    pub value: f64,
    pub d_hk: f64,
    pub d_theta: f64,
    pub d_re_theta: f64,
}

impl AmplificationRate {
    fn to_array(self) -> [f64; 4] {
        [self.value, self.d_hk, self.d_theta, self.d_re_theta]
    }

    fn from_array([value, d_hk, d_theta, d_re_theta]: [f64; 4]) -> Self {
        Self { value, d_hk, d_theta, d_re_theta }
    }
}

/// Envelope amplification rate at kinematic shape `hk`, momentum thickness
/// `theta` and Reynolds number `re_theta`.
pub fn amplification_rate(model: Amplification, hk: f64, theta: f64, re_theta: f64) -> AmplificationRate {
    /// Half-width, in log10(Rθ), of the ramp that switches amplification on.
    const RAMP: f64 = 0.08;
    let hmi = 1.0 / (hk - 1.0);
    let hmi_hk = -hmi * hmi;

    // log10 of the critical Rθ against H for Falkner–Skan profiles
    let aa = 2.492 * hmi.powf(0.43);
    let aa_hk = (aa / hmi) * 0.43 * hmi_hk;
    let bb = (14.0 * hmi - 9.24).tanh();
    let bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;
    let log_rt_crit = aa + 0.7 * (bb + 1.0);
    let log_rt_crit_hk = aa_hk + 0.7 * bb_hk;

    let log_rt = re_theta.log10();
    let log_rt_rt = 1.0 / (LN_10 * re_theta);

    let mut rate = AmplificationRate::default();
    if log_rt >= log_rt_crit - RAMP {
        // smooth cubic ramp switching amplification on around Rθ_crit
        let rnorm = (log_rt - (log_rt_crit - RAMP)) / (2.0 * RAMP);
        let (ramp, ramp_hk, ramp_rt) = if rnorm >= 1.0 {
            (1.0, 0.0, 0.0)
        } else {
            let ramp_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;
            (
                3.0 * rnorm * rnorm - 2.0 * rnorm.powi(3),
                ramp_rn * -log_rt_crit_hk / (2.0 * RAMP),
                ramp_rn * log_rt_rt / (2.0 * RAMP),
            )
        };

        // dN/dRθ
        let arg = 3.87 * hmi - 2.52;
        let ex = (-arg * arg).exp();
        let ex_hk = ex * (-2.0 * arg * 3.87 * hmi_hk);
        let dn_drt = 0.028 * (hk - 1.0) - 0.0345 * ex;
        let dn_drt_hk = 0.028 - 0.0345 * ex_hk;

        // θ·dRθ/dξ
        let (mut af, mut af_hmi) =
            (-0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi.powi(3), 2.7 - 11.0 * hmi + 9.0 * hmi * hmi);
        if model == Amplification::Modified {
            let e = (-20.0 * hmi).exp();
            af += 0.1 * e;
            af_hmi -= 2.0 * e;
        }
        let af_hk = af_hmi * hmi_hk;

        let value = af * dn_drt / theta * ramp;
        rate = AmplificationRate {
            value,
            d_hk: (af_hk * dn_drt / theta + af * dn_drt_hk / theta) * ramp + af * dn_drt / theta * ramp_hk,
            d_theta: -value / theta,
            d_re_theta: af * dn_drt / theta * ramp_rt,
        };
    }

    if model == Amplification::Modified {
        rate = blend_separated(rate, hk, theta, log_rt, log_rt_rt);
    }
    rate
}

/// Blends in the Orr–Sommerfeld maximum amplification rate for `3.5 < Hk < 4`.
fn blend_separated(attached: AmplificationRate, hk: f64, theta: f64, log_rt: f64, log_rt_rt: f64) -> AmplificationRate {
    const HK1: f64 = 3.5;
    const HK2: f64 = 4.0;
    if hk < HK1 {
        return attached;
    }
    let hnorm = (hk - HK1) / (HK2 - HK1);
    let (blend, blend_hk) = if hnorm >= 1.0 {
        (1.0, 0.0)
    } else {
        (3.0 * hnorm * hnorm - 2.0 * hnorm.powi(3), (6.0 * hnorm - 6.0 * hnorm * hnorm) / (HK2 - HK1))
    };

    let e = (-0.15 * (hk - 5.0)).exp();
    let log_rt0 = 0.30 + 0.35 * e;
    let log_rt0_hk = -0.35 * e * 0.15;
    let tnr = (1.2 * (log_rt - log_rt0)).tanh();
    let tnr_rt = (1.0 - tnr * tnr) * 1.2 * log_rt_rt;
    let tnr_hk = -(1.0 - tnr * tnr) * 1.2 * log_rt0_hk;

    let mut separated = [
        (0.086 * tnr - 0.25 / (hk - 1.0).powf(1.5)) / theta,
        (0.086 * tnr_hk + 1.5 * 0.25 / (hk - 1.0).powf(2.5)) / theta,
        0.0,
        0.086 * tnr_rt / theta,
    ];
    separated[2] = -separated[0] / theta;
    if separated[0] < 0.0 {
        separated = [0.0; 4];
    }

    let attached = attached.to_array();
    let mut rate = [0.0; 4];
    for k in 0..4 {
        rate[k] = blend * separated[k] + (1.0 - blend) * attached[k];
    }
    rate[1] += blend_hk * (separated[0] - attached[0]);
    AmplificationRate::from_array(rate)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_partials(f: impl Fn(f64, f64, f64) -> Fit, hk: f64, rt: f64, msq: f64) {
        let h = 1e-6;
        let fit = f(hk, rt, msq);
        let num_hk = (f(hk + h, rt, msq).value - f(hk - h, rt, msq).value) / (2.0 * h);
        let num_rt = (f(hk, rt + h * rt, msq).value - f(hk, rt - h * rt, msq).value) / (2.0 * h * rt);
        assert!((fit.d_hk - num_hk).abs() < 1e-5 * (1.0 + num_hk.abs()), "d/dhk {} vs {}", fit.d_hk, num_hk);
        assert!(
            (fit.d_re_theta - num_rt).abs() < 1e-5 * (1.0 + num_rt.abs()),
            "d/drt {} vs {}",
            fit.d_re_theta,
            num_rt
        );
    }

    #[test]
    fn analytic_partials_match_finite_differences() {
        for &(hk, rt) in &[(2.2, 300.0), (3.0, 1500.0), (4.6, 8000.0), (6.0, 250.0)] {
            check_partials(energy_shape_turbulent, hk, rt, 0.1);
            check_partials(|h, r, m| skin_friction_turbulent(h, r, m, 1.0), hk, rt, 0.1);
            check_partials(|h, r, _| skin_friction_laminar(h, r), hk, rt, 0.0);
            check_partials(|h, r, _| dissipation_laminar(h, r), hk, rt, 0.0);
            check_partials(|h, r, _| dissipation_laminar_wake(h, r), hk, rt, 0.0);
            check_partials(|h, _, _| energy_shape_laminar(h), hk, rt, 0.0);
        }
    }

    #[test]
    fn amplification_switches_on_above_critical_reynolds() {
        let below = amplification_rate(Amplification::Envelope, 2.6, 1e-3, 50.0);
        let above = amplification_rate(Amplification::Envelope, 2.6, 1e-3, 2000.0);
        assert_eq!(below.value, 0.0);
        assert!(above.value > 0.0);
    }
}
