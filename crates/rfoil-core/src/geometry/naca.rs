//! NACA 4- and 5-digit section generators.

use crate::error::{Error, Result};

/// Exponent of the chordwise point distribution's bunching toward the LE and TE.
const BUNCHING: f64 = 1.5;

/// Coordinates of a NACA section ordered
/// from the TE over the upper surface to the LE, then back along the lower surface to the TE.
pub struct NacaSection {
    pub name: String,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

/// Generates a NACA 4- or 5-digit airfoil from its designation, e.g. `"2412"` or `"23012"`.
pub fn naca(designation: &str, points_per_side: usize) -> Result<NacaSection> {
    let digits = designation.trim().trim_start_matches(|c: char| c.is_alphabetic() || c.is_whitespace());
    let code: u32 = digits.parse().map_err(|_| Error::InvalidNaca(designation.into()))?;
    let camber = match digits.len() {
        4 => four_digit_camber(code),
        5 => five_digit_camber(code).ok_or_else(|| Error::InvalidNaca(designation.into()))?,
        _ => return Err(Error::InvalidNaca(designation.into())),
    };
    let thickness = (code % 100) as f64 / 100.0;
    Ok(assemble(format!("NACA {digits}"), thickness, camber, points_per_side.max(3)))
}

type CamberLine = Box<dyn Fn(f64) -> f64>;

/// NACA `MPxx`: maximum camber `M`% at `P`/10 of the chord.
fn four_digit_camber(code: u32) -> CamberLine {
    let max_camber = (code / 1000) as f64 / 100.0;
    let x_max = ((code / 100) % 10) as f64 / 10.0;
    Box::new(move |x| {
        if x < x_max {
            max_camber / (x_max * x_max) * (2.0 * x_max * x - x * x)
        } else {
            max_camber / (1.0 - x_max).powi(2) * ((1.0 - 2.0 * x_max) + 2.0 * x_max * x - x * x)
        }
    })
}

/// NACA `2Pxx` with a standard (non-reflexed) mean line; `m` is the
/// transition point of the cubic and `k` its scale.
fn five_digit_camber(code: u32) -> Option<CamberLine> {
    let (m, k) = match code / 100 {
        210 => (0.0580, 361.4),
        220 => (0.1260, 51.64),
        230 => (0.2025, 15.957),
        240 => (0.2900, 6.643),
        250 => (0.3910, 3.230),
        _ => return None,
    };
    Some(Box::new(move |x| {
        if x < m {
            k / 6.0 * (x.powi(3) - 3.0 * m * x * x + m * m * (3.0 - m) * x)
        } else {
            k / 6.0 * m.powi(3) * (1.0 - x)
        }
    }))
}

/// Half thickness of the symmetric NACA 00xx section with thickness ratio `t`
/// (with the closed-TE coefficient 0.1015).
fn half_thickness(x: f64, t: f64) -> f64 {
    (0.29690 * x.sqrt() - 0.12600 * x - 0.35160 * x.powi(2) + 0.28430 * x.powi(3) - 0.10150 * x.powi(4)) * t / 0.20
}

fn assemble(name: String, thickness: f64, camber: CamberLine, points_per_side: usize) -> NacaSection {
    let exponent = BUNCHING + 1.0;
    let stations: Vec<f64> = (0..points_per_side)
        .map(|i| {
            let frac = i as f64 / (points_per_side - 1) as f64;
            if i == points_per_side - 1 {
                1.0
            } else {
                1.0 - exponent * frac * (1.0 - frac).powf(BUNCHING) - (1.0 - frac).powf(exponent)
            }
        })
        .collect();

    let upper = stations.iter().rev().map(|&x| (x, camber(x) + half_thickness(x, thickness)));
    let lower = stations.iter().skip(1).map(|&x| (x, camber(x) - half_thickness(x, thickness)));
    let (x, y) = upper.chain(lower).unzip();
    NacaSection { name, x, y }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn naca0012_is_symmetric_and_closed_ccw() {
        let s = naca("0012", 61).unwrap();
        assert_eq!(s.x.len(), 121);
        assert_eq!(s.x[0], 1.0);
        assert_eq!(s.x[60], 0.0);
        assert!((s.y[30] + s.y[90]).abs() < 1e-14);
        assert!(s.y[30] > 0.0);
    }

    #[test]
    fn rejects_bad_designations() {
        assert!(naca("12", 10).is_err());
        assert!(naca("26012", 10).is_err());
        assert!(naca("NACA 23012", 10).is_ok());
    }
}
