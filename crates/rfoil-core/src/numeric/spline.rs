//! Cubic Hermite splines parameterised by arc length.
//!
//! A spline of `x(s)` is stored as node values `x[i]` and node slopes
//! `dx/ds[i]`.  Repeated `s` values mark segment joints (corners), across
//! which the slope may be discontinuous.

use super::linalg::solve_tridiagonal;

/// End condition used when computing spline slopes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum End {
    /// Zero second derivative ("natural" spline).
    ZeroCurvature,
    /// Zero third derivative (end interval is a parabola).
    ZeroJerk,
    /// Prescribed slope `dx/ds`.
    Slope(f64),
}

/// Node slopes `dx/ds` of the spline through `(s, x)` with the given end conditions.
pub fn slopes(x: &[f64], s: &[f64], first: End, last: End) -> Vec<f64> {
    let n = x.len();
    let mut slope = vec![0.0; n];
    if n < 2 {
        return slope;
    }
    let (mut diag, mut lower, mut upper) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);

    // continuity of the second derivative at interior nodes
    for i in 1..n - 1 {
        let ds_before = s[i] - s[i - 1];
        let ds_after = s[i + 1] - s[i];
        lower[i] = ds_after;
        diag[i] = 2.0 * (ds_before + ds_after);
        upper[i] = ds_before;
        slope[i] = 3.0 * ((x[i + 1] - x[i]) * ds_before / ds_after + (x[i] - x[i - 1]) * ds_after / ds_before);
    }

    let secant = |i: usize, j: usize| (x[j] - x[i]) / (s[j] - s[i]);
    match first {
        End::ZeroCurvature => (diag[0], upper[0], slope[0]) = (2.0, 1.0, 3.0 * secant(0, 1)),
        End::ZeroJerk => (diag[0], upper[0], slope[0]) = (1.0, 1.0, 2.0 * secant(0, 1)),
        End::Slope(v) => (diag[0], upper[0], slope[0]) = (1.0, 0.0, v),
    }
    // two zero-jerk ends on a single interval would be singular
    let last = if n == 2 && first == End::ZeroJerk && last == End::ZeroJerk { End::ZeroCurvature } else { last };
    let end = n - 1;
    match last {
        End::ZeroCurvature => (lower[end], diag[end], slope[end]) = (1.0, 2.0, 3.0 * secant(end - 1, end)),
        End::ZeroJerk => (lower[end], diag[end], slope[end]) = (1.0, 1.0, 2.0 * secant(end - 1, end)),
        End::Slope(v) => (lower[end], diag[end], slope[end]) = (0.0, 1.0, v),
    }

    solve_tridiagonal(&mut diag, &lower, &mut upper, &mut slope);
    slope
}

/// Like [`slopes`] with zero-jerk ends, but splines each segment separately.
/// Segments are delimited by identical successive `s` values.
pub fn segmented_slopes(x: &[f64], s: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut slope = vec![0.0; n];
    let mut start = 0;
    for i in 1..n.saturating_sub(2) {
        if s[i] == s[i + 1] {
            let segment = start..=i;
            slope[segment.clone()].copy_from_slice(&slopes(
                &x[segment.clone()],
                &s[segment],
                End::ZeroJerk,
                End::ZeroJerk,
            ));
            start = i + 1;
        }
    }
    slope[start..].copy_from_slice(&slopes(&x[start..], &s[start..], End::ZeroJerk, End::ZeroJerk));
    slope
}

/// Cumulative chord length of the polyline `(x, y)`.
pub fn arc_length(x: &[f64], y: &[f64]) -> Vec<f64> {
    let mut s = Vec::with_capacity(x.len());
    let mut total = 0.0;
    for i in 0..x.len() {
        if i > 0 {
            total += (x[i] - x[i - 1]).hypot(y[i] - y[i - 1]);
        }
        s.push(total);
    }
    s
}

/// A spline `x(s)` borrowed as node values, node slopes and node parameters.
#[derive(Clone, Copy, Debug)]
pub struct Spline<'a> {
    pub x: &'a [f64],
    pub slope: &'a [f64],
    pub s: &'a [f64],
}

/// The cubic on one spline interval in Hermite form.
struct Interval {
    /// Local parameter in `[0, 1]`.
    t: f64,
    /// Interval length.
    ds: f64,
    /// End values.
    x0: f64,
    x1: f64,
    /// End-slope deviations from the secant, scaled by `ds`.
    c0: f64,
    c1: f64,
}

impl Interval {
    /// `ds·dx/ds` at the local parameter.
    fn scaled_derivative(&self) -> f64 {
        let t = self.t;
        self.x1 - self.x0 + (1.0 - 4.0 * t + 3.0 * t * t) * self.c0 + t * (3.0 * t - 2.0) * self.c1
    }

    /// `ds²·d²x/ds²` at the local parameter.
    fn scaled_second_derivative(&self) -> f64 {
        let t = self.t;
        (6.0 * t - 4.0) * self.c0 + (6.0 * t - 2.0) * self.c1
    }
}

impl<'a> Spline<'a> {
    pub fn new(x: &'a [f64], slope: &'a [f64], s: &'a [f64]) -> Self {
        Self { x, slope, s }
    }

    /// Index `i` with `s[i-1] <= at < s[i]`, clamped to the node range.
    fn interval_index(&self, at: f64) -> usize {
        let (mut lo, mut hi) = (0, self.s.len() - 1);
        while hi - lo > 1 {
            let mid = (hi + lo) / 2;
            if at < self.s[mid] {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        hi
    }

    fn interval(&self, at: f64) -> Interval {
        let i = self.interval_index(at);
        let ds = self.s[i] - self.s[i - 1];
        let dx = self.x[i] - self.x[i - 1];
        Interval {
            t: (at - self.s[i - 1]) / ds,
            ds,
            x0: self.x[i - 1],
            x1: self.x[i],
            c0: ds * self.slope[i - 1] - dx,
            c1: ds * self.slope[i] - dx,
        }
    }

    /// `x(at)`.
    pub fn value(&self, at: f64) -> f64 {
        let Interval { t, x0, x1, c0, c1, .. } = self.interval(at);
        t * x1 + (1.0 - t) * x0 + (t - t * t) * ((1.0 - t) * c0 - t * c1)
    }

    /// `dx/ds(at)`.
    pub fn derivative(&self, at: f64) -> f64 {
        let interval = self.interval(at);
        interval.scaled_derivative() / interval.ds
    }

    /// `d²x/ds²(at)`.
    pub fn second_derivative(&self, at: f64) -> f64 {
        let interval = self.interval(at);
        interval.scaled_second_derivative() / (interval.ds * interval.ds)
    }

    /// Solves `x(s) = target` by Newton iteration from `guess`.
    /// Returns `guess` unchanged if the iteration fails.
    pub fn invert(&self, target: f64, guess: f64) -> f64 {
        let span = self.s[self.s.len() - 1] - self.s[0];
        let mut s = guess;
        for _ in 0..10 {
            let step = -(self.value(s) - target) / self.derivative(s);
            s += step;
            if (step / span).abs() < 1.0e-5 {
                return s;
            }
        }
        guess
    }
}

/// Signed curvature of the planar curve `(x(s), y(s))` at `at`.
pub fn curvature(x: Spline, y: Spline, at: f64) -> f64 {
    let (ix, iy) = (x.interval(at), y.interval(at));
    let (xd, yd) = (ix.scaled_derivative(), iy.scaled_derivative());
    let speed = xd.hypot(yd).max(0.001 * ix.ds);
    (xd * iy.scaled_second_derivative() - yd * ix.scaled_second_derivative()) / speed.powi(3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproduces_cubic_with_exact_slopes() {
        let s: Vec<f64> = (0..11).map(|i| i as f64 * 0.1).collect();
        let x: Vec<f64> = s.iter().map(|&t| t * t * t).collect();
        let slope = slopes(&x, &s, End::Slope(0.0), End::Slope(3.0));
        let spline = Spline::new(&x, &slope, &s);
        for &t in &[0.05, 0.33, 0.71, 0.99] {
            assert!((spline.value(t) - t.powi(3)).abs() < 1e-4);
            assert!((spline.derivative(t) - 3.0 * t * t).abs() < 2e-3);
        }
    }

    #[test]
    fn inverts_monotone_spline() {
        let s: Vec<f64> = (0..21).map(|i| i as f64 * 0.05).collect();
        let x: Vec<f64> = s.iter().map(|&t| t.sin()).collect();
        let slope = slopes(&x, &s, End::ZeroCurvature, End::ZeroCurvature);
        let si = Spline::new(&x, &slope, &s).invert(0.6_f64.sin(), 0.5);
        assert!((si - 0.6).abs() < 1e-4);
    }
}
