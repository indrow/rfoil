//! The paneled body: airfoil nodes followed by wake nodes.

use std::f64::consts::PI;

use crate::geometry::{ContourPoint, Curve};
use crate::numeric::spline::segmented_slopes;

/// Geometry of the panel that closes a blunt trailing edge.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TrailingEdge {
    /// TE gap projected normal to the TE bisector (the base thickness).
    pub normal_gap: f64,
    /// TE gap projected along the TE bisector.
    pub stream_gap: f64,
    /// Total TE gap.
    pub gap: f64,
    pub sharp: bool,
}

impl TrailingEdge {
    /// Fractions of the TE panel strength carried as source and as vorticity.
    pub fn source_vortex_weights(&self) -> (f64, f64) {
        if self.sharp {
            (1.0, 0.0)
        } else {
            (self.normal_gap / self.gap, self.stream_gap / self.gap)
        }
    }
}

/// Panel nodes.  Indices `0..n` are the airfoil, running from the TE over the upper
/// surface to the LE, then back along the lower surface to the TE; indices `n..` are
/// the wake, starting just behind the TE.
#[derive(Clone, Debug)]
pub struct Panels {
    /// Number of airfoil nodes.
    pub n: usize,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Arc length; continues along the wake.
    pub s: Vec<f64>,
    /// Unit outward normals.
    pub nx: Vec<f64>,
    pub ny: Vec<f64>,
    /// Orientation of the panel starting at each node.
    pub angle: Vec<f64>,
    /// Spline slopes of the airfoil contour (airfoil nodes only).
    pub dx_ds: Vec<f64>,
    pub dy_ds: Vec<f64>,
    pub le: ContourPoint,
    /// Midpoint of the two TE nodes.
    pub te_midpoint: (f64, f64),
    pub chord: f64,
    pub te: TrailingEdge,
}

impl Panels {
    pub fn from_curve(c: &Curve) -> Self {
        let n = c.len();
        let le = c.leading_edge();
        let te_midpoint = c.trailing_edge();
        let chord = (te_midpoint.0 - le.x).hypot(te_midpoint.1 - le.y);
        let (nx, ny) = normals(&c.x, &c.y, &c.s);
        let mut panels = Self {
            n,
            x: c.x.clone(),
            y: c.y.clone(),
            s: c.s.clone(),
            nx,
            ny,
            angle: vec![0.0; n],
            dx_ds: c.dx_ds.clone(),
            dy_ds: c.dy_ds.clone(),
            le,
            te_midpoint,
            chord,
            te: TrailingEdge::default(),
        };
        panels.te = panels.trailing_edge_geometry();
        panels.angle = panels.panel_angles();
        panels
    }

    /// Total number of nodes (airfoil + wake).
    pub fn len(&self) -> usize {
        self.x.len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Airfoil contour of this paneling.
    pub fn curve(&self) -> Curve {
        let n = self.n;
        Curve::with_arc_length(self.x[..n].to_vec(), self.y[..n].to_vec(), self.s[..n].to_vec())
    }

    /// Distance between nodes `a` and `b`.
    pub fn distance(&self, a: usize, b: usize) -> f64 {
        (self.x[a] - self.x[b]).hypot(self.y[a] - self.y[b])
    }

    /// Drops any wake nodes.
    pub fn clear_wake(&mut self) {
        let n = self.n;
        for v in [&mut self.x, &mut self.y, &mut self.s, &mut self.nx, &mut self.ny, &mut self.angle] {
            v.truncate(n);
        }
    }

    /// Appends a wake node with its normal; its panel angle is set later.
    pub(crate) fn push_wake_node(&mut self, x: f64, y: f64, s: f64, nx: f64, ny: f64) {
        self.x.push(x);
        self.y.push(y);
        self.s.push(s);
        self.nx.push(nx);
        self.ny.push(ny);
        self.angle.push(0.0);
    }

    /// Chordwise coordinates `(x/c, y/c)` of a point, in the chord-line frame.
    pub fn chord_coordinates(&self, x: f64, y: f64) -> (f64, f64) {
        let (cx, cy) = (self.te_midpoint.0 - self.le.x, self.te_midpoint.1 - self.le.y);
        let chord_sq = cx * cx + cy * cy;
        let (dx, dy) = (x - self.le.x, y - self.le.y);
        ((dx * cx + dy * cy) / chord_sq, (dy * cx - dx * cy) / chord_sq)
    }

    fn trailing_edge_geometry(&self) -> TrailingEdge {
        let last = self.n - 1;
        let (gap_x, gap_y) = (self.x[0] - self.x[last], self.y[0] - self.y[last]);
        // mean of the two TE tangents (both pointing downstream)
        let bisector_x = 0.5 * (-self.dx_ds[0] + self.dx_ds[last]);
        let bisector_y = 0.5 * (-self.dy_ds[0] + self.dy_ds[last]);
        let gap = gap_x.hypot(gap_y);
        TrailingEdge {
            normal_gap: bisector_x * gap_y - bisector_y * gap_x,
            stream_gap: bisector_x * gap_x + bisector_y * gap_y,
            gap,
            sharp: gap < 1.0e-4 * self.chord,
        }
    }

    fn panel_angles(&self) -> Vec<f64> {
        let n = self.n;
        let mut angle: Vec<f64> = (0..n - 1)
            .map(|i| {
                let (sx, sy) = (self.x[i + 1] - self.x[i], self.y[i + 1] - self.y[i]);
                if sx == 0.0 && sy == 0.0 {
                    (-self.ny[i]).atan2(-self.nx[i])
                } else {
                    sx.atan2(-sy)
                }
            })
            .collect();
        // the TE panel runs from the last node back to the first
        angle.push(if self.te.sharp {
            PI
        } else {
            let (sx, sy) = (self.x[0] - self.x[n - 1], self.y[0] - self.y[n - 1]);
            (-sx).atan2(sy) + PI
        });
        angle
    }
}

/// Unit normals at the nodes of a splined contour, averaged at corners.
fn normals(x: &[f64], y: &[f64], s: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let dx_ds = segmented_slopes(x, s);
    let dy_ds = segmented_slopes(y, s);
    let (mut nx, mut ny): (Vec<f64>, Vec<f64>) = dx_ds
        .iter()
        .zip(&dy_ds)
        .map(|(&dx, &dy)| {
            let len = dx.hypot(dy);
            (dy / len, -dx / len)
        })
        .unzip();
    for i in 0..s.len() - 1 {
        if s[i] == s[i + 1] {
            let (ax, ay) = (0.5 * (nx[i] + nx[i + 1]), 0.5 * (ny[i] + ny[i + 1]));
            let len = ax.hypot(ay);
            (nx[i], ny[i], nx[i + 1], ny[i + 1]) = (ax / len, ay / len, ax / len, ay / len);
        }
    }
    (nx, ny)
}
