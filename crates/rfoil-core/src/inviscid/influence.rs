//! Streamfunction influence of the airfoil vortex sheet and of the
//! airfoil/wake source sheets at a field point.
//!
//! The airfoil carries a linearly varying vortex strength `γ` and, for the
//! viscous coupling, a source strength `σ` that varies quadratically within
//! each panel (represented by two linear half-panels).  A constant-strength
//! source/vortex panel closes a blunt trailing edge.
//!
//! "Normal" derivatives `d/dn` are taken along the probe direction; with the
//! probe direction set to the surface normal, `dψ/dn` is the tangential speed.

use std::f64::consts::PI;

use super::panels::Panels;

const ONE_OVER_4PI: f64 = 0.25 / PI;
const ONE_OVER_2PI: f64 = 0.5 / PI;

/// A point at which the influence is evaluated.
#[derive(Clone, Copy, Debug)]
pub struct Probe {
    pub x: f64,
    pub y: f64,
    /// Direction along which the derivative `dψ/dn` is taken.
    pub nx: f64,
    pub ny: f64,
    /// Node index if the probe sits on a panel node (avoids self-singularities).
    pub node: Option<usize>,
}

impl Probe {
    /// Probe at node `i` using that node's normal.
    pub fn at_node(p: &Panels, i: usize) -> Self {
        Self { x: p.x[i], y: p.y[i], nx: p.nx[i], ny: p.ny[i], node: Some(i) }
    }
}

/// Uniform onset flow.
#[derive(Clone, Copy, Debug)]
pub struct UniformFlow {
    pub speed: f64,
    pub alpha: f64,
}

/// Influence of the airfoil vortex sheet (plus onset flow) at a probe.
#[derive(Clone, Debug)]
pub struct VortexInfluence {
    /// `dψ/dγ_j`.
    pub dpsi_dgamma: Vec<f64>,
    /// `d(dψ/dn)/dγ_j`: speed along the probe per unit node vorticity.
    pub dq_dgamma: Vec<f64>,
    /// `dψ/dn` of the given vorticity and onset flow.
    pub q: f64,
    /// `dψ/dn` of the α = 0° and 90° unit solutions.
    pub q_unit: [f64; 2],
}

/// Influence of a source sheet at a probe.
#[derive(Clone, Debug)]
pub struct SourceInfluence {
    /// `dψ/dσ_j`.
    pub dpsi_dsigma: Vec<f64>,
    /// `d(dψ/dn)/dσ_j`.
    pub dq_dsigma: Vec<f64>,
}

impl SourceInfluence {
    fn zeros(len: usize) -> Self {
        Self { dpsi_dsigma: vec![0.0; len], dq_dsigma: vec![0.0; len] }
    }
}

/// Coefficients of a panel whose strength varies linearly from `a` to `b`:
/// `ψ = sum·(a + b) + diff·(b − a)`, and likewise for `dψ/dn`.
#[derive(Clone, Copy)]
struct LinearCoefficients {
    sum: f64,
    diff: f64,
    sum_dn: f64,
    diff_dn: f64,
}

/// Coefficients of the constant-strength TE panel per unit source and vortex strength.
#[derive(Clone, Copy)]
struct TeCoefficients {
    source: f64,
    vortex: f64,
    source_dn: f64,
    vortex_dn: f64,
}

/// One end of a panel (or half-panel) seen from the probe, in the panel frame.
#[derive(Clone, Copy)]
struct PanelEnd {
    /// Probe position along the panel, measured from this end.
    x: f64,
    /// Squared distance to the probe, its logarithm and the angle it subtends.
    r_sq: f64,
    ln_r_sq: f64,
    theta: f64,
}

/// Probe position expressed in the local frame of one panel.
struct PanelFrame {
    start: PanelEnd,
    end: PanelEnd,
    /// Probe height above the panel.
    y: f64,
    /// Derivatives of the along-panel and normal coordinates along the probe direction.
    dx_dn: f64,
    dy_dn: f64,
    /// Reflection sign for probes below the sheet, and the matching angle offset.
    reflect: f64,
    branch_offset: f64,
}

impl PanelFrame {
    /// Frame of the panel from node `a` to node `b`.  `on_sheet` is true when
    /// the probe lies on the sheet being integrated; `branch_sign` selects
    /// the arctangent branch used off the sheet.
    fn new(p: &Panels, a: usize, b: usize, probe: &Probe, on_sheet: bool, branch_sign: f64) -> Option<Self> {
        let length = p.distance(a, b);
        if length == 0.0 {
            return None;
        }
        let (sx, sy) = ((p.x[b] - p.x[a]) / length, (p.y[b] - p.y[a]) / length);
        let (rx1, ry1) = (probe.x - p.x[a], probe.y - p.y[a]);
        let (rx2, ry2) = (probe.x - p.x[b], probe.y - p.y[b]);
        let y = sx * ry1 - sy * rx1;

        let reflect = if on_sheet { 1.0 } else { 1.0_f64.copysign(y) };
        let branch_offset = branch_sign * (0.5 - 0.5 * reflect) * PI;
        let end = |x: f64, r_sq: f64, node: usize| {
            let (ln_r_sq, theta) = if probe.node != Some(node) && r_sq > 0.0 {
                (r_sq.ln(), (reflect * x).atan2(reflect * y) + branch_offset)
            } else {
                (0.0, 0.0)
            };
            PanelEnd { x, r_sq, ln_r_sq, theta }
        };

        Some(Self {
            start: end(sx * rx1 + sy * ry1, rx1 * rx1 + ry1 * ry1, a),
            end: end(sx * rx2 + sy * ry2, rx2 * rx2 + ry2 * ry2, b),
            y,
            dx_dn: sx * probe.nx + sy * probe.ny,
            dy_dn: sx * probe.ny - sy * probe.nx,
            reflect,
            branch_offset,
        })
    }

    /// Streamfunction coefficients of a linear-strength vortex panel (times 4π).
    fn vortex(&self) -> LinearCoefficients {
        let (x1, x2, y) = (self.start.x, self.end.x, self.y);
        let (g1, g2) = (self.start.ln_r_sq, self.end.ln_r_sq);
        let (t1, t2) = (self.start.theta, self.end.theta);
        let inv_dx = 1.0 / (x1 - x2);

        let sum = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + y * (t1 - t2);
        let diff = ((x1 + x2) * sum + 0.5 * (self.end.r_sq * g2 - self.start.r_sq * g1 + x1 * x1 - x2 * x2)) * inv_dx;

        let (sum_x1, sum_x2, sum_y) = (0.5 * g1, -0.5 * g2, t1 - t2);
        let diff_x1 = ((x1 + x2) * sum_x1 + sum - x1 * g1 - diff) * inv_dx;
        let diff_x2 = ((x1 + x2) * sum_x2 + sum + x2 * g2 + diff) * inv_dx;
        let diff_y = ((x1 + x2) * sum_y - y * (g1 - g2)) * inv_dx;

        // both panel ends move together with the probe in the panel frame
        LinearCoefficients {
            sum,
            diff,
            sum_dn: (sum_x1 + sum_x2) * self.dx_dn + sum_y * self.dy_dn,
            diff_dn: (diff_x1 + diff_x2) * self.dx_dn + diff_y * self.dy_dn,
        }
    }

    /// Source and vortex coefficients of the constant-strength TE panel (times 2π).
    fn te_panel(&self, panel_angle: f64) -> TeCoefficients {
        let (x1, x2, y) = (self.start.x, self.end.x, self.y);
        let (g1, g2) = (self.start.ln_r_sq, self.end.ln_r_sq);
        let (t1, t2) = (self.start.theta, self.end.theta);
        TeCoefficients {
            source: 0.5 * y * (g1 - g2) + x2 * (t2 - panel_angle) - x1 * (t1 - panel_angle),
            vortex: 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + y * (t1 - t2),
            source_dn: (t2 - t1) * self.dx_dn + 0.5 * (g1 - g2) * self.dy_dn,
            vortex_dn: (0.5 * g1 - 0.5 * g2) * self.dx_dn + (t1 - t2) * self.dy_dn,
        }
    }

    /// Adds the influence of a quadratic source distribution on this panel,
    /// built from two linear half-panels meeting at the midpoint.
    fn add_sources(&self, panel_angle: f64, stencil: SourceStencil, out: &mut SourceInfluence) {
        let x_mid = 0.5 * (self.start.x + self.end.x);
        let r_sq_mid = x_mid * x_mid + self.y * self.y;
        let mid = PanelEnd {
            x: x_mid,
            r_sq: r_sq_mid,
            ln_r_sq: r_sq_mid.ln(),
            theta: (self.reflect * x_mid).atan2(self.reflect * self.y) + self.branch_offset,
        };

        let SourceStencil { before, start, end, after, inv_span_at_start, inv_length, inv_span_at_end } = stencil;
        // each half-panel's end strengths are linear in the node strengths
        // through the source slopes: (node, ∂(sum)/∂σ, ∂(diff)/∂σ)
        let halves = [
            (
                self.start,
                mid,
                [
                    (before, -inv_span_at_start, inv_span_at_start),
                    (start, -inv_length, -inv_length),
                    (end, inv_length + inv_span_at_start, inv_length - inv_span_at_start),
                ],
            ),
            (
                mid,
                self.end,
                [
                    (start, -(inv_span_at_end + inv_length), -(inv_span_at_end - inv_length)),
                    (end, inv_length, -inv_length),
                    (after, inv_span_at_end, inv_span_at_end),
                ],
            ),
        ];
        for (a, b, weights) in halves {
            let c = source_half_panel(a, b, self.y, panel_angle, self.dx_dn, self.dy_dn);
            for (node, w_sum, w_diff) in weights {
                out.dpsi_dsigma[node] += ONE_OVER_4PI * (c.sum * w_sum + c.diff * w_diff);
                out.dq_dsigma[node] += ONE_OVER_4PI * (c.sum_dn * w_sum + c.diff_dn * w_diff);
            }
        }
    }
}

/// Nodes around a source panel that define its quadratic strength variation.
#[derive(Clone, Copy)]
struct SourceStencil {
    /// Node before the panel start, the panel ends, and the node after the end
    /// (repeated at the ends of the sheet).
    before: usize,
    start: usize,
    end: usize,
    after: usize,
    /// `1/|end − before|`, `1/|end − start|` and `1/|after − start|`.
    inv_span_at_start: f64,
    inv_length: f64,
    inv_span_at_end: f64,
}

impl SourceStencil {
    /// Stencil of the panel starting at node `start` on the sheet `first..=last`.
    fn new(p: &Panels, start: usize, first: usize, last: usize) -> Self {
        let end = start + 1;
        let before = if start == first { start } else { start - 1 };
        let after = if end == last { end } else { end + 1 };
        Self {
            before,
            start,
            end,
            after,
            inv_span_at_start: 1.0 / p.distance(end, before),
            inv_length: 1.0 / p.distance(start, end),
            inv_span_at_end: 1.0 / p.distance(after, start),
        }
    }
}

/// Coefficients of a linear-strength source half-panel from `a` to `b` (times 4π).
fn source_half_panel(a: PanelEnd, b: PanelEnd, y: f64, panel_angle: f64, dx_dn: f64, dy_dn: f64) -> LinearCoefficients {
    let (ta, tb) = (a.theta - panel_angle, b.theta - panel_angle);
    let inv_dx = 1.0 / (a.x - b.x);
    let sum = b.x * tb - a.x * ta + 0.5 * y * (a.ln_r_sq - b.ln_r_sq);
    let diff = ((a.x + b.x) * sum + a.r_sq * ta - b.r_sq * tb + (b.x - a.x) * y) * inv_dx;

    let (sum_xa, sum_xb, sum_y) = (-ta, tb, 0.5 * (a.ln_r_sq - b.ln_r_sq));
    let diff_xa = ((a.x + b.x) * sum_xa + sum + 2.0 * a.x * ta - diff) * inv_dx;
    let diff_xb = ((a.x + b.x) * sum_xb + sum - 2.0 * b.x * tb + diff) * inv_dx;
    let diff_y = ((a.x + b.x) * sum_y + 2.0 * (b.x - a.x + y * (ta - tb))) * inv_dx;

    LinearCoefficients {
        sum,
        diff,
        sum_dn: (sum_xa + sum_xb) * dx_dn + sum_y * dy_dn,
        diff_dn: (diff_xa + diff_xb) * dx_dn + diff_y * dy_dn,
    }
}

/// Influence of the airfoil vortex sheet, TE panel and onset flow at `probe`.
/// `gamma` is the current vorticity, `unit_gamma` the α = 0°, 90° unit solutions.
/// If `with_sources` is set, the airfoil source-sheet sensitivities are also returned.
pub fn airfoil(
    p: &Panels,
    gamma: &[f64],
    unit_gamma: &[Vec<f64>; 2],
    onset: UniformFlow,
    probe: &Probe,
    with_sources: bool,
) -> (VortexInfluence, Option<SourceInfluence>) {
    let n = p.n;
    let on_airfoil = matches!(probe.node, Some(i) if i < n);
    let mut v = VortexInfluence { dpsi_dgamma: vec![0.0; n], dq_dgamma: vec![0.0; n], q: 0.0, q_unit: [0.0; 2] };
    let mut sources = with_sources.then(|| SourceInfluence::zeros(n));

    for a in 0..n - 1 {
        let b = a + 1;
        let Some(frame) = PanelFrame::new(p, a, b, probe, on_airfoil, 1.0) else {
            continue;
        };
        if let Some(sources) = sources.as_mut() {
            frame.add_sources(p.angle[a], SourceStencil::new(p, a, 0, n - 1), sources);
        }
        let c = frame.vortex();
        let q_of = |g: &[f64]| ONE_OVER_4PI * ((g[b] + g[a]) * c.sum_dn + (g[b] - g[a]) * c.diff_dn);
        v.dpsi_dgamma[a] += ONE_OVER_4PI * (c.sum - c.diff);
        v.dpsi_dgamma[b] += ONE_OVER_4PI * (c.sum + c.diff);
        v.dq_dgamma[a] += ONE_OVER_4PI * (c.sum_dn - c.diff_dn);
        v.dq_dgamma[b] += ONE_OVER_4PI * (c.sum_dn + c.diff_dn);
        v.q += q_of(gamma);
        v.q_unit[0] += q_of(&unit_gamma[0]);
        v.q_unit[1] += q_of(&unit_gamma[1]);
    }

    add_te_panel(p, gamma, unit_gamma, probe, on_airfoil, &mut v);

    let (sa, ca) = onset.alpha.sin_cos();
    v.q += onset.speed * (ca * probe.ny - sa * probe.nx);
    v.q_unit[0] += onset.speed * probe.ny;
    v.q_unit[1] -= onset.speed * probe.nx;
    (v, sources)
}

/// Adds the constant-strength source/vortex panel bridging the trailing-edge gap,
/// whose strengths follow from the vorticity jump across it.
fn add_te_panel(
    p: &Panels,
    gamma: &[f64],
    unit_gamma: &[Vec<f64>; 2],
    probe: &Probe,
    on_airfoil: bool,
    v: &mut VortexInfluence,
) {
    let (a, b) = (p.n - 1, 0);
    let negligible = (p.s[a] - p.s[0]) * 1.0e-5;
    if p.distance(a, b) < negligible {
        return;
    }
    let Some(frame) = PanelFrame::new(p, a, b, probe, on_airfoil, 1.0) else {
        return;
    };
    let c = frame.te_panel(p.angle[a]);
    let (source_weight, vortex_weight) = p.te.source_vortex_weights();

    // TE source strength = ½·w_s·Δγ, TE vortex strength = −½·w_v·Δγ
    let q_of = |g: &[f64]| {
        let jump = g[b] - g[a];
        ONE_OVER_2PI * (c.source_dn * 0.5 * source_weight * jump - c.vortex_dn * 0.5 * vortex_weight * jump)
    };
    v.q += q_of(gamma);
    v.q_unit[0] += q_of(&unit_gamma[0]);
    v.q_unit[1] += q_of(&unit_gamma[1]);

    let dpsi = ONE_OVER_2PI * 0.5 * (c.source * source_weight - c.vortex * vortex_weight);
    let dq = ONE_OVER_2PI * 0.5 * (c.source_dn * source_weight - c.vortex_dn * vortex_weight);
    v.dpsi_dgamma[a] -= dpsi;
    v.dpsi_dgamma[b] += dpsi;
    v.dq_dgamma[a] -= dq;
    v.dq_dgamma[b] += dq;
}

/// Sensitivities to the wake source sheet at `probe`, indexed by node
/// (only the wake entries `n..` are non-zero).
pub fn wake_sources(p: &Panels, probe: &Probe) -> SourceInfluence {
    let (first, last) = (p.n, p.len() - 1);
    let on_wake = matches!(probe.node, Some(i) if i >= first);
    let mut out = SourceInfluence::zeros(p.len());
    for a in first..last {
        if let Some(frame) = PanelFrame::new(p, a, a + 1, probe, on_wake, -1.0) {
            frame.add_sources(p.angle[a], SourceStencil::new(p, a, first, last), &mut out);
        }
    }
    out
}
