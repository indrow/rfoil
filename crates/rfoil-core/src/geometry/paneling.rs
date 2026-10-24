//! Curvature-weighted panel node distribution.
//!
//! Nodes are placed so that `(1 + C·κ)·Δs` is equal on both sides of every
//! node, where `κ` is a smoothed, normalised surface curvature.  Artificial
//! curvature is added at the trailing edge and in optional refinement zones
//! to cluster panels there.

use super::curve::Curve;
use crate::numeric::spline::segmented_slopes;
use crate::numeric::{solve_tridiagonal, Spline};

/// User-facing paneling parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct PanelingOptions {
    /// Number of panel nodes.
    pub nodes: usize,
    /// Curvature attraction: 0 gives uniform spacing, ~1 is typical.
    pub bunching: f64,
    /// TE / LE panel density ratio.
    pub te_le_ratio: f64,
    /// Refined-area / LE panel density ratio.
    pub refine_ratio: f64,
    /// Top-side refinement zone, in x/c.
    pub top_refine: (f64, f64),
    /// Bottom-side refinement zone, in x/c.
    pub bottom_refine: (f64, f64),
}

impl Default for PanelingOptions {
    fn default() -> Self {
        Self {
            nodes: 160,
            bunching: 1.0,
            te_le_ratio: 0.15,
            refine_ratio: 0.2,
            top_refine: (1.0, 1.0),
            bottom_refine: (1.0, 1.0),
        }
    }
}

/// Temporary node over-sampling factor for robust convergence.
const OVERSAMPLE: usize = 5;
/// Ratio of the TE panel length to its neighbour.
const TE_PANEL_RATIO: f64 = 0.667;

/// Leading-edge location on the input contour.
#[derive(Clone, Copy)]
struct LeadingEdge {
    s: f64,
    /// Index of the first of two coincident nodes if the LE is a corner.
    corner: Option<usize>,
}

/// Redistributes the nodes of `buffer` according to `opts`.
pub fn panel(buffer: &Curve, opts: &PanelingOptions) -> Curve {
    let s = buffer.leading_edge().s;
    let corner = (0..buffer.len() - 1).find(|&i| buffer.s[i] == s && buffer.s[i + 1] == s);
    let le = LeadingEdge { s, corner };

    let field = CurvatureField::new(buffer, opts, le);
    let oversampled = distribute(buffer, &field, opts, le);

    let s_nodes: Vec<f64> = oversampled.iter().step_by(OVERSAMPLE).copied().collect();
    let (x, y) = insert_corners(buffer, &s_nodes);
    Curve::new(x, y)
}

/// The three diagonals of a tridiagonal matrix.
struct Tridiagonal {
    lower: Vec<f64>,
    diag: Vec<f64>,
    upper: Vec<f64>,
}

impl Tridiagonal {
    fn zeros(n: usize) -> Self {
        Self { lower: vec![0.0; n], diag: vec![0.0; n], upper: vec![0.0; n] }
    }

    fn identity(n: usize) -> Self {
        Self { diag: vec![1.0; n], ..Self::zeros(n) }
    }

    fn set_row(&mut self, i: usize, (lower, diag, upper): (f64, f64, f64)) {
        (self.lower[i], self.diag[i], self.upper[i]) = (lower, diag, upper);
    }

    /// Makes row `i` read `x_i = rhs_i`.
    fn pin(&mut self, i: usize) {
        self.set_row(i, (0.0, 1.0, 0.0));
    }

    /// Solves rows `range` as an independent system, in place.
    fn solve_range(&mut self, range: std::ops::Range<usize>, rhs: &mut [f64]) {
        let r = range.clone();
        solve_tridiagonal(&mut self.diag[r.clone()], &self.lower[r.clone()], &mut self.upper[r], &mut rhs[range]);
    }
}

/// Normalised, smoothed curvature `κ(s)` that drives the panel density.
struct CurvatureField {
    s: Vec<f64>,
    kappa: Vec<f64>,
    dkappa_ds: Vec<f64>,
}

/// Row of the implicit smoothing operator `(1 − L²·d²/ds²)` at a node with
/// spacings `ds_before`, `ds_after` and mean spacing `ds_mean`.
fn smoothing_row(l_sq: f64, ds_before: f64, ds_after: f64, ds_mean: f64) -> (f64, f64, f64) {
    (-l_sq / ds_before / ds_mean, l_sq * (1.0 / ds_after + 1.0 / ds_before) / ds_mean + 1.0, -l_sq / ds_after / ds_mean)
}

impl CurvatureField {
    fn new(c: &Curve, opts: &PanelingOptions, le: LeadingEdge) -> Self {
        let n = c.len();
        let s = &c.s;
        let half_length = 0.5 * c.length();
        let curvature = |at: f64| c.curvature_at(at).abs() * half_length;

        let mut kappa: Vec<f64> = s.iter().map(|&si| curvature(si)).collect();
        let kappa_le = curvature(le.s);

        // LE curvature averaged over a few points within one LE radius
        let kappa_le_mean = if le.corner.is_some() {
            10.0
        } else {
            let half_count = 3;
            let total: f64 = (-half_count..=half_count)
                .map(|j| curvature(le.s + j as f64 / half_count as f64 * half_length / kappa_le.max(20.0)))
                .sum();
            total / (2 * half_count + 1) as f64
        };

        // artificial TE curvature
        kappa[0] = kappa_le_mean * opts.te_le_ratio;
        kappa[n - 1] = kappa[0];

        let smoothing_length = (1.0 / kappa_le_mean.max(20.0)).max(0.25 / (opts.nodes / 2) as f64);
        let l_sq = (smoothing_length * half_length).powi(2);
        let mut system = Tridiagonal::identity(n);
        for i in 1..n - 1 {
            let ds_before = s[i] - s[i - 1];
            let ds_after = s[i + 1] - s[i];
            if ds_before == 0.0 || ds_after == 0.0 {
                continue; // keep the curvature at corners
            }
            system.set_row(i, smoothing_row(l_sq, ds_before, ds_after, 0.5 * (s[i + 1] - s[i - 1])));
        }

        Self::pin_leading_edge(&mut system, &mut kappa, s, le, kappa_le, l_sq);

        // artificial curvature inside the refinement zones
        let le_point = c.point_at(le.s);
        let (xte, yte) = c.trailing_edge();
        let chord_sq = (xte - le_point.x).powi(2) + (yte - le_point.y).powi(2);
        for i in 1..n - 1 {
            let x_over_c =
                ((c.x[i] - le_point.x) * (xte - le_point.x) + (c.y[i] - le_point.y) * (yte - le_point.y)) / chord_sq;
            let (from, to) = if s[i] < le.s { opts.top_refine } else { opts.bottom_refine };
            if x_over_c > from && x_over_c < to {
                system.pin(i);
                kappa[i] = kappa_le * opts.refine_ratio;
            }
        }

        // a LE corner splits the contour into two independent smoothing problems
        match le.corner {
            None => system.solve_range(0..n, &mut kappa),
            Some(j) => {
                system.solve_range(0..j + 1, &mut kappa);
                system.solve_range(j + 1..n, &mut kappa);
            }
        }

        let kappa_max = kappa.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        kappa.iter_mut().for_each(|v| *v /= kappa_max);
        let dkappa_ds = segmented_slopes(&kappa, s);
        Self { s: s.clone(), kappa, dkappa_ds }
    }

    /// Holds the curvature at the leading edge at its true value.  If the LE
    /// falls between two nodes, the smoothing rows of those nodes are
    /// rewritten so the LE acts as an extra node with fixed curvature.
    fn pin_leading_edge(
        system: &mut Tridiagonal,
        kappa: &mut [f64],
        s: &[f64],
        le: LeadingEdge,
        kappa_le: f64,
        l_sq: f64,
    ) {
        let at_corner = |i: usize| Some(i) == le.corner || Some(i) == le.corner.map(|j| j + 1);
        for i in 1..s.len() - 1 {
            if s[i] == le.s || at_corner(i) {
                system.pin(i);
                kappa[i] = kappa_le;
            } else if s[i - 1] < le.s && s[i] > le.s {
                let (ds_before, ds_after, ds_mean) = (s[i - 1] - s[i - 2], le.s - s[i - 1], 0.5 * (le.s - s[i - 2]));
                let (lower, diag, _) = smoothing_row(l_sq, ds_before, ds_after, ds_mean);
                system.set_row(i - 1, (lower, diag, 0.0));
                kappa[i - 1] += l_sq * kappa_le / (ds_after * ds_mean);

                let (ds_before, ds_after, ds_mean) = (s[i] - le.s, s[i + 1] - s[i], 0.5 * (s[i + 1] - le.s));
                let (_, diag, upper) = smoothing_row(l_sq, ds_before, ds_after, ds_mean);
                system.set_row(i, (0.0, diag, upper));
                kappa[i] += l_sq * kappa_le / (ds_before * ds_mean);
                break;
            }
        }
    }

    /// Curvature and its slope at `s`.
    fn at(&self, s: f64) -> (f64, f64) {
        let spline = Spline::new(&self.kappa, &self.dkappa_ds, &self.s);
        (spline.value(s), spline.derivative(s))
    }
}

/// RMS curvature over an interval with its partials w.r.t. the end positions.
struct IntervalCurvature {
    value: f64,
    d_start: f64,
    d_end: f64,
}

impl IntervalCurvature {
    fn new((k0, dk0): (f64, f64), (k1, dk1): (f64, f64)) -> Self {
        let value = k0.hypot(k1);
        if value == 0.0 {
            Self { value: 0.0, d_start: 0.0, d_end: 0.0 }
        } else {
            Self { value, d_start: dk0 * k0 / value, d_end: dk1 * k1 / value }
        }
    }
}

/// Solves for the over-sampled node positions by Newton iteration.
fn distribute(c: &Curve, field: &CurvatureField, opts: &PanelingOptions, le: LeadingEdge) -> Vec<f64> {
    let count = OVERSAMPLE * (opts.nodes - 1) + 1;
    // the TE panel length constraint, written as a ratio of interval differences
    let te_factor = (TE_PANEL_RATIO - 1.0) * 2.0 + 1.0;
    let attraction = 6.0 * opts.bunching;
    let (mut s, le_node) = initial_nodes(c, opts, le, count, te_factor);

    for _ in 0..20 {
        let mut system = Tridiagonal::zeros(count);
        let mut rhs = vec![0.0; count];

        // equidistribute (1 + C·κ)·Δs between neighbouring intervals
        let mut k_here = field.at(s[1]);
        let mut before = IntervalCurvature::new(field.at(s[0]), k_here);
        for i in 1..count - 1 {
            let ds_before = s[i] - s[i - 1];
            let ds_after = s[i] - s[i + 1];
            let k_next = field.at(s[i + 1]);
            let after = IntervalCurvature::new(k_here, k_next);
            let weight_before = attraction * before.value + 1.0;
            let weight_after = attraction * after.value + 1.0;

            system.set_row(
                i,
                (
                    -weight_before + attraction * ds_before * before.d_start,
                    weight_after + weight_before + attraction * (ds_after * after.d_start + ds_before * before.d_end),
                    -weight_after + attraction * ds_after * after.d_end,
                ),
            );
            rhs[i] = -(ds_after * weight_after + ds_before * weight_before);

            k_here = k_next;
            before = after;
        }

        // fixed end points
        system.pin(0);
        system.pin(count - 1);

        // TE panels a fixed fraction of their neighbours
        rhs[1] = -((s[1] - s[0]) + te_factor * (s[1] - s[2]));
        system.set_row(1, (-1.0, 1.0 + te_factor, -te_factor));
        let i = count - 2;
        rhs[i] = -((s[i] - s[i + 1]) + te_factor * (s[i] - s[i - 1]));
        system.set_row(i, (-te_factor, 1.0 + te_factor, -1.0));

        if let Some(i) = le_node {
            system.pin(i);
            rhs[i] = le.s - s[i];
        }

        system.solve_range(0..count, &mut rhs);

        let (relaxation, max_step) = limit_step(&s, &rhs);
        for i in 1..count - 1 {
            s[i] += relaxation * rhs[i];
        }
        if max_step < 1.0e-3 {
            break;
        }
    }
    s
}

/// Uniformly spaced nodes with shortened TE panels.  With a LE corner, each
/// side is spaced separately and the node landing on the corner is returned.
fn initial_nodes(
    c: &Curve,
    opts: &PanelingOptions,
    le: LeadingEdge,
    count: usize,
    te_factor: f64,
) -> (Vec<f64>, Option<usize>) {
    let (s_first, s_last) = (c.s[0], c.s[c.len() - 1]);
    let mut s = vec![0.0; count];
    let mut le_node = None;
    match le.corner {
        None => {
            let ds = (s_last - s_first) / ((count - 3) as f64 + 2.0 * te_factor);
            for (i, v) in s.iter_mut().enumerate().take(count - 1).skip(1) {
                *v = s_first + ds * ((i - 1) as f64 + te_factor);
            }
        }
        Some(j) => {
            let upper_nodes = opts.nodes * (j + 1) / c.len();
            let upper_count = OVERSAMPLE * (upper_nodes - 1) + 1;
            let ds = (le.s - s_first) / ((upper_count - 2) as f64 + te_factor);
            for i in 1..upper_count {
                s[i] = s_first + ds * ((i - 1) as f64 + te_factor);
            }
            let lower_count = count - upper_count + 1;
            let ds = (s_last - le.s) / ((lower_count - 2) as f64 + te_factor);
            for i in 1..lower_count - 1 {
                s[i - 1 + upper_count] = le.s + ds * ((i - 1) as f64 + te_factor);
            }
            le_node = Some(upper_count - 1);
        }
    }
    s[0] = s_first;
    s[count - 1] = s_last;
    (s, le_node)
}

/// Under-relaxation keeping every interval within 0.2×–4× of its length,
/// and the largest node displacement of the full step.
fn limit_step(s: &[f64], step: &[f64]) -> (f64, f64) {
    let mut relaxation = 1.0_f64;
    let mut max_step = 0.0_f64;
    for i in 0..s.len() - 1 {
        let ds = s[i + 1] - s[i];
        let dds = step[i + 1] - step[i];
        let ratio = 1.0 + relaxation * dds / ds;
        if ratio > 4.0 {
            relaxation = 3.0 * ds / dds;
        }
        if ratio < 0.2 {
            relaxation = -0.8 * ds / dds;
        }
        max_step = max_step.max(step[i].abs());
    }
    (relaxation, max_step)
}

/// Evaluates the node coordinates and inserts any buffer corners (doubled points).
fn insert_corners(c: &Curve, s_nodes: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let is_corner = |ib: usize| c.s[ib] == c.s[ib + 1];
    let mut s = s_nodes.to_vec();
    for ib in (0..c.len() - 1).filter(|&ib| is_corner(ib)) {
        let corner = c.s[ib];
        if let Some(i) = s.iter().position(|&v| v > corner) {
            s.insert(i, corner);
            // keep the panels adjacent to the corner comparable in size
            if i >= 2 {
                s[i - 1] = 0.5 * (s[i] + s[i - 2]);
            }
            if i + 2 < s.len() {
                s[i + 1] = 0.5 * (s[i] + s[i + 2]);
            }
        }
    }

    // corner nodes take the exact buffer corner coordinates
    let mut x: Vec<f64> = s.iter().map(|&v| c.x_at(v)).collect();
    let mut y: Vec<f64> = s.iter().map(|&v| c.y_at(v)).collect();
    for ib in (0..c.len() - 1).filter(|&ib| is_corner(ib)) {
        if let Some(i) = s.iter().position(|&v| v == c.s[ib]) {
            (x[i], y[i]) = (c.x[ib], c.y[ib]);
        }
    }
    (x, y)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Airfoil;

    #[test]
    fn panels_naca0012_with_requested_node_count() {
        let foil = Airfoil::naca("0012").unwrap();
        let opts = PanelingOptions::default();
        let paneled = panel(&foil.curve(), &opts);
        assert_eq!(paneled.len(), opts.nodes);
        // nodes cluster at the LE: smallest panel near the middle of the list
        let ds: Vec<f64> = paneled.s.windows(2).map(|w| w[1] - w[0]).collect();
        let smallest = ds.iter().enumerate().min_by(|a, b| a.1.total_cmp(b.1)).unwrap().0;
        assert!((smallest as i64 - 79).abs() < 5, "smallest panel at {smallest}");
        assert!((paneled.chord() - 1.0).abs() < 1e-3);
    }
}
