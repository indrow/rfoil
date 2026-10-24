//! Wake trajectory: a streamline traced downstream from the trailing edge.

use super::influence::{self, Probe, UniformFlow};
use super::panels::Panels;

/// Number of wake nodes for an airfoil with `n` nodes.
pub fn wake_node_count(n: usize) -> usize {
    n / 8 + 2
}

/// Replaces the wake of `p` with a streamline of the flow given by the
/// vorticity `gamma`, extending `length` chords behind the trailing edge.
pub fn trace(p: &mut Panels, gamma: &[f64], onset: UniformFlow, length: f64) {
    p.clear_wake();
    let n = p.n;
    let count = wake_node_count(n);
    let first_step = 0.5 * (p.s[1] - p.s[0] + p.s[n - 1] - p.s[n - 2]);
    let spacing = geometric_spacing(first_step, length * p.chord, count);
    let zero = vec![0.0; n];
    let no_unit = [zero.clone(), zero];

    // the local velocity direction gives the normal and angle of the next wake panel
    let flow_normal = |p: &Panels, i: usize| {
        let speed_along = |nx, ny| {
            let probe = Probe { x: p.x[i], y: p.y[i], nx, ny, node: Some(i) };
            influence::airfoil(p, gamma, &no_unit, onset, &probe, false).0.q
        };
        let (psi_x, psi_y) = (speed_along(1.0, 0.0), speed_along(0.0, 1.0));
        let len = psi_x.hypot(psi_y);
        ((-psi_x / len, -psi_y / len), psi_y.atan2(psi_x))
    };

    // first node a tiny distance behind the TE, normal to the TE bisector
    let last = n - 1;
    let (sx, sy) = (0.5 * (p.dy_ds[last] - p.dy_ds[0]), 0.5 * (p.dx_ds[0] - p.dx_ds[last]));
    let len = sx.hypot(sy);
    let (nx, ny) = (sx / len, sy / len);
    let (xte, yte) = p.te_midpoint;
    p.push_wake_node(xte - 1.0e-4 * ny, yte + 1.0e-4 * nx, p.s[last], nx, ny);

    let mut next_normal;
    (next_normal, p.angle[n]) = flow_normal(p, n);
    for k in 1..count {
        let i = n + k;
        let ds = spacing[k] - spacing[k - 1];
        let (nx, ny) = next_normal;
        p.push_wake_node(p.x[i - 1] - ds * ny, p.y[i - 1] + ds * nx, p.s[i - 1] + ds, nx, ny);
        if k < count - 1 {
            (next_normal, p.angle[i]) = flow_normal(p, i);
        }
    }
}

/// Inviscid speeds of the α = 0°, 90° unit solutions over airfoil and wake.
pub fn unit_speeds(p: &Panels, unit_gamma: &[Vec<f64>; 2], speed: f64) -> [Vec<f64>; 2] {
    let n = p.n;
    let mut q = [unit_gamma[0][..n].to_vec(), unit_gamma[1][..n].to_vec()];
    let onset = UniformFlow { speed, alpha: 0.0 };
    let zero = vec![0.0; n];
    for i in n..p.len() {
        let q_unit = if i == n {
            [q[0][n - 1], q[1][n - 1]] // the first wake node sits at the TE
        } else {
            influence::airfoil(p, &zero, unit_gamma, onset, &Probe::at_node(p, i), false).0.q_unit
        };
        q[0].push(q_unit[0]);
        q[1].push(q_unit[1]);
    }
    q
}

/// `count` points from 0 to `total` whose spacing grows geometrically from `first`.
pub fn geometric_spacing(first: f64, total: f64, count: usize) -> Vec<f64> {
    let intervals = (count - 1) as f64;
    let length_ratio = total / first;
    let inv_intervals = 1.0 / intervals;

    // initial guess for the stretching ratio from the quadratic truncation of
    // Σ rᵏ = length_ratio, then Newton on the intervals-th root of that sum
    let cubic = intervals * (intervals - 1.0) * (intervals - 2.0) / 6.0;
    let quadratic = intervals * (intervals - 1.0) / 2.0;
    let constant = intervals - length_ratio;
    let mut ratio = if count == 3 {
        -constant / quadratic + 1.0
    } else {
        (-quadratic + (quadratic * quadratic - 4.0 * cubic * constant).max(0.0).sqrt()) / (2.0 * cubic) + 1.0
    };
    if ratio != 1.0 {
        for _ in 0..100 {
            let power = ratio.powf(intervals);
            let sum = (power - 1.0) / (ratio - 1.0);
            let residual = sum.powf(inv_intervals) - length_ratio.powf(inv_intervals);
            let slope = inv_intervals * sum.powf(inv_intervals) * (intervals * ratio.powf(intervals - 1.0) - sum)
                / (power - 1.0);
            let step = -residual / slope;
            ratio += step;
            if step.abs() < 1.0e-5 {
                break;
            }
        }
    }

    let mut s = Vec::with_capacity(count);
    let (mut position, mut ds) = (0.0, first);
    for _ in 0..count {
        s.push(position);
        position += ds;
        ds *= ratio;
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geometric_spacing_hits_the_end() {
        let s = geometric_spacing(0.01, 1.0, 22);
        assert_eq!(s.len(), 22);
        assert!((s[1] - 0.01).abs() < 1e-14);
        assert!((s[21] - 1.0).abs() < 1e-4);
    }
}
