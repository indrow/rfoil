//! Raw ("buffer") airfoil coordinates as read from a file or generated.

use super::curve::{Curve, SectionShape};
use super::naca;
use crate::error::{Error, Result};

/// Points per side used for generated NACA sections.
const NACA_POINTS_PER_SIDE: usize = 120;

/// An airfoil as a named list of contour points, ordered
/// from the TE over the upper surface to the LE, then back along the lower surface to the TE.
#[derive(Clone, Debug)]
pub struct Airfoil {
    pub name: String,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

impl Airfoil {
    /// Builds an airfoil, reversing the point order if it was given clockwise.
    pub fn new(name: impl Into<String>, x: Vec<f64>, y: Vec<f64>) -> Result<Self> {
        if x.len() != y.len() || x.len() < 4 {
            return Err(Error::TooFewPoints(x.len().min(y.len())));
        }
        let mut foil = Self { name: name.into(), x, y };
        if foil.curve().area() < 0.0 {
            foil.x.reverse();
            foil.y.reverse();
        }
        Ok(foil)
    }

    /// NACA 4- or 5-digit section, e.g. `Airfoil::naca("2412")`.
    pub fn naca(designation: &str) -> Result<Self> {
        let section = naca::naca(designation, NACA_POINTS_PER_SIDE)?;
        Self::new(section.name, section.x, section.y)
    }

    /// Parses the common coordinate-file layouts: plain `x y` pairs (optionally
    /// preceded by a name line), ISES/MSES files (only the first element is
    /// read) and Lednicer files (upper and lower surfaces each listed
    /// from LE to TE).
    pub fn parse(text: &str) -> Result<Self> {
        let mut name = String::new();
        let mut pairs = Vec::new();
        for line in text.lines().map(str::trim).filter(|l| !l.is_empty() && !l.starts_with('#')) {
            match parse_pair(line) {
                Some(p) => pairs.push(p),
                None if pairs.is_empty() && name.is_empty() => name = line.to_string(),
                None => {}
            }
        }

        // ISES files carry a line of domain parameters after the name
        if pairs.len() > 1 && line_has_four_numbers(text) {
            pairs.remove(0);
        }
        // MSES element separator
        if let Some(end) = pairs.iter().position(|&(x, y)| x == 999.0 && y == 999.0) {
            pairs.truncate(end);
        }

        let (x, y) = match lednicer_counts(&pairs) {
            Some((upper, lower)) => lednicer_to_contour(&pairs[1..], upper, lower)?,
            None => pairs.into_iter().unzip(),
        };
        let name = if name.is_empty() { "Unnamed".to_string() } else { name };
        Self::new(name, x, y)
    }

    /// Writes the coordinates as a name line followed by `x y` pairs.
    pub fn to_dat(&self) -> String {
        let mut out = format!("{}\n", self.name);
        for (x, y) in self.x.iter().zip(&self.y) {
            out.push_str(&format!("{x:12.7} {y:12.7}\n"));
        }
        out
    }

    pub fn curve(&self) -> Curve {
        Curve::new(self.x.clone(), self.y.clone())
    }

    pub fn shape(&self) -> SectionShape {
        self.curve().thickness_and_camber()
    }
}

fn parse_pair(line: &str) -> Option<(f64, f64)> {
    let mut tokens = line.split(|c: char| c.is_whitespace() || c == ',').filter(|t| !t.is_empty());
    let x = tokens.next()?.parse().ok()?;
    let y = tokens.next()?.parse().ok()?;
    Some((x, y))
}

/// ISES files have, as their second non-blank line, four numbers (the domain box).
fn line_has_four_numbers(text: &str) -> bool {
    text.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .nth(1)
        .map(|l| {
            let tokens: Vec<_> = l.split_whitespace().collect();
            tokens.len() == 4 && tokens.iter().all(|t| t.parse::<f64>().is_ok())
        })
        .unwrap_or(false)
}

/// Lednicer files begin with the point counts of each surface, e.g. `17. 17.`.
fn lednicer_counts(pairs: &[(f64, f64)]) -> Option<(usize, usize)> {
    let &(upper, lower) = pairs.first()?;
    let is_count = |v: f64| v >= 2.0 && v.fract() == 0.0;
    (is_count(upper) && is_count(lower) && (upper + lower) as usize == pairs.len() - 1)
        .then_some((upper as usize, lower as usize))
}

/// Joins Lednicer surfaces (both listed from LE to TE) into one contour that
/// starts and ends at the TE.
fn lednicer_to_contour(points: &[(f64, f64)], upper_count: usize, lower_count: usize) -> Result<(Vec<f64>, Vec<f64>)> {
    let (upper, lower) = points.split_at(upper_count);
    if lower.len() != lower_count {
        return Err(Error::Parse("Lednicer surface counts do not match the data".into()));
    }
    // the upper surface reversed to run from TE to LE, then the lower surface without its duplicated LE point
    let skip = usize::from(lower.first() == upper.first());
    Ok(upper.iter().rev().chain(&lower[skip..]).copied().unzip())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_labelled_file_and_fixes_orientation() {
        let text = "My foil\n0.0 0.0\n0.5 -0.05\n1.0 0.0\n0.5 0.06\n0.0 0.0\n";
        let foil = Airfoil::parse(text).unwrap();
        assert_eq!(foil.name, "My foil");
        assert_eq!(foil.x.len(), 5);
        assert!(foil.curve().area() > 0.0);
    }

    #[test]
    fn parses_lednicer_file() {
        let text = "LED\n3. 3.\n\n0.0 0.0\n0.5 0.06\n1.0 0.0\n\n0.0 0.0\n0.5 -0.05\n1.0 0.0\n";
        let foil = Airfoil::parse(text).unwrap();
        assert_eq!(foil.x, vec![1.0, 0.5, 0.0, 0.5, 1.0]);
        assert_eq!(foil.y[1], 0.06);
    }
}
