//! A closed airfoil contour splined against its arc length.

use crate::numeric::spline::{self, arc_length, segmented_slopes, Spline};

/// Contour `x(s), y(s)` running
/// from the TE over the upper surface to the LE, then back along the lower surface to the TE.
#[derive(Clone, Debug)]
pub struct Curve {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Arc length at each node.
    pub s: Vec<f64>,
    /// Spline slopes at each node.
    pub dx_ds: Vec<f64>,
    pub dy_ds: Vec<f64>,
}

/// A point on the contour together with its arc-length parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContourPoint {
    pub s: f64,
    pub x: f64,
    pub y: f64,
}

/// Thickness and camber summary of a section.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SectionShape {
    pub thickness: f64,
    pub x_thickness: f64,
    pub camber: f64,
    pub x_camber: f64,
}

/// Unit vector along the chord line, from the leading edge to the TE midpoint.
struct ChordFrame {
    le: ContourPoint,
    length: f64,
    dir: (f64, f64),
}

impl ChordFrame {
    /// Distance along the chord line from the leading edge.
    fn along(&self, x: f64, y: f64) -> f64 {
        (x - self.le.x) * self.dir.0 + (y - self.le.y) * self.dir.1
    }

    /// Distance normal to the chord line (positive above it).
    fn normal(&self, x: f64, y: f64) -> f64 {
        (y - self.le.y) * self.dir.0 - (x - self.le.x) * self.dir.1
    }
}

impl Curve {
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Self {
        let s = arc_length(&x, &y);
        Self::with_arc_length(x, y, s)
    }

    /// Splines a contour whose arc-length parameter is already known.
    pub fn with_arc_length(x: Vec<f64>, y: Vec<f64>, s: Vec<f64>) -> Self {
        let dx_ds = segmented_slopes(&x, &s);
        let dy_ds = segmented_slopes(&y, &s);
        Self { x, y, s, dx_ds, dy_ds }
    }

    pub fn len(&self) -> usize {
        self.x.len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Total arc length.
    pub fn length(&self) -> f64 {
        self.s[self.len() - 1] - self.s[0]
    }

    pub fn x_spline(&self) -> Spline<'_> {
        Spline::new(&self.x, &self.dx_ds, &self.s)
    }

    pub fn y_spline(&self) -> Spline<'_> {
        Spline::new(&self.y, &self.dy_ds, &self.s)
    }

    pub fn x_at(&self, s: f64) -> f64 {
        self.x_spline().value(s)
    }

    pub fn y_at(&self, s: f64) -> f64 {
        self.y_spline().value(s)
    }

    pub fn point_at(&self, s: f64) -> ContourPoint {
        ContourPoint { s, x: self.x_at(s), y: self.y_at(s) }
    }

    pub fn curvature_at(&self, s: f64) -> f64 {
        spline::curvature(self.x_spline(), self.y_spline(), s)
    }

    /// Midpoint of the two trailing-edge nodes.
    pub fn trailing_edge(&self) -> (f64, f64) {
        let last = self.len() - 1;
        (0.5 * (self.x[0] + self.x[last]), 0.5 * (self.y[0] + self.y[last]))
    }

    /// Leading edge: the point whose tangent is normal to the line joining it to the TE.
    pub fn leading_edge(&self) -> ContourPoint {
        self.point_at(self.leading_edge_s())
    }

    fn leading_edge_s(&self) -> f64 {
        let n = self.len();
        let (xte, yte) = self.trailing_edge();
        let tolerance = self.length() * 1.0e-5;

        // first guess: the first node where the contour turns back toward the TE
        let i = (2..n - 2)
            .find(|&i| {
                let to_te = (self.x[i] - xte, self.y[i] - yte);
                let step = (self.x[i + 1] - self.x[i], self.y[i + 1] - self.y[i]);
                to_te.0 * step.0 + to_te.1 * step.1 < 0.0
            })
            .unwrap_or(n - 2);
        let guess = self.s[i];
        if self.s[i] == self.s[i - 1] {
            return guess; // sharp LE corner
        }

        // Newton iteration on (r − r_te)·dr/ds = 0
        let (xs, ys) = (self.x_spline(), self.y_spline());
        let mut s = guess;
        for _ in 0..50 {
            let (dx, dy) = (xs.derivative(s), ys.derivative(s));
            let (ddx, ddy) = (xs.second_derivative(s), ys.second_derivative(s));
            let (cx, cy) = (xs.value(s) - xte, ys.value(s) - yte);

            let residual = cx * dx + cy * dy;
            let jacobian = dx * dx + dy * dy + cx * ddx + cy * ddy;
            let limit = 0.02 * (cx + cy).abs();
            let step = (-residual / jacobian).clamp(-limit, limit);
            s += step;
            if step.abs() < tolerance {
                return s;
            }
        }
        guess
    }

    fn chord_frame(&self, le_s: f64) -> ChordFrame {
        let le = self.point_at(le_s);
        let (xte, yte) = self.trailing_edge();
        let length = (xte - le.x).hypot(yte - le.y);
        ChordFrame { le, length, dir: ((xte - le.x) / length, (yte - le.y) / length) }
    }

    /// Chord length measured from the leading edge to the TE midpoint.
    pub fn chord(&self) -> f64 {
        self.chord_frame(self.leading_edge_s()).length
    }

    /// Arc length of the point on the opposite surface at the same chordwise
    /// station as the point at `s`.  `le_s` is the leading-edge arc length.
    pub fn opposite_point(&self, s: f64, le_s: f64) -> f64 {
        let last = self.len() - 1;
        let span = self.length();
        let frame = self.chord_frame(le_s);

        let (near_end, far_end) = if s < le_s { (self.s[0], self.s[last]) } else { (self.s[last], self.s[0]) };
        let fraction = (s - le_s) / (near_end - le_s);
        let linear_guess = le_s + fraction * (far_end - le_s);
        if fraction.abs() <= 1.0e-5 {
            return le_s;
        }

        let (xs, ys) = (self.x_spline(), self.y_spline());
        let target = frame.along(xs.value(s), ys.value(s));
        let mut opposite = linear_guess;
        for _ in 0..12 {
            let residual = frame.along(xs.value(opposite), ys.value(opposite)) - target;
            let slope = xs.derivative(opposite) * frame.dir.0 + ys.derivative(opposite) * frame.dir.1;
            if (residual / span).abs() < 1.0e-5 {
                return opposite;
            }
            if slope == 0.0 {
                break;
            }
            let step = -residual / slope;
            opposite += step;
            if (step / span).abs() < 1.0e-5 {
                return opposite;
            }
        }
        linear_guess
    }

    /// Maximum thickness and camber, sampled at the contour nodes.
    pub fn thickness_and_camber(&self) -> SectionShape {
        let le_s = self.leading_edge_s();
        let frame = self.chord_frame(le_s);

        let mut shape = SectionShape::default();
        for i in 0..self.len() {
            let opposite = self.opposite_point(self.s[i], le_s);
            let x_opposite = self.x_at(opposite);
            let here = frame.normal(self.x[i], self.y[i]);
            let there = frame.normal(x_opposite, self.y_at(opposite));
            let camber = 0.5 * (here + there);
            let thickness = (here - there).abs();
            if camber.abs() > shape.camber.abs() {
                (shape.camber, shape.x_camber) = (camber, x_opposite);
            }
            if thickness > shape.thickness {
                (shape.thickness, shape.x_thickness) = (thickness, x_opposite);
            }
        }
        shape
    }

    /// Enclosed cross-sectional area (positive for counter-clockwise ordering).
    pub fn area(&self) -> f64 {
        let n = self.len();
        (0..n)
            .map(|i| {
                let j = (i + 1) % n;
                0.5 * (self.y[i] + self.y[j]) * (self.x[i] - self.x[j])
            })
            .sum()
    }
}
