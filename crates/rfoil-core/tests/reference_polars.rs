//! Regression against the reference polars described in the README
//! (double precision, 160 panels, Ncrit = 9, free transition).

use rfoil_core::{Airfoil, Analysis, OperatingPoint, Settings};

/// `(alpha, CL, CD, CM)` as printed in the reference polar file.
type Row = (f64, f64, f64, f64);

fn sweep(naca: &str, settings: Settings, start: f64, end: f64, step: f64) -> Vec<OperatingPoint> {
    let mut analysis = Analysis::new(Airfoil::naca(naca).unwrap()).unwrap();
    analysis.set_settings(Settings { max_iterations: 100, ..settings });
    analysis.sweep_alpha(start, end, step, |_| true)
}

fn assert_matches(points: &[OperatingPoint], reference: &[Row]) {
    for &(alpha, cl, cd, cm) in reference {
        let p = points
            .iter()
            .find(|p| (p.alpha - alpha).abs() < 1e-9)
            .unwrap_or_else(|| panic!("no point at alpha = {alpha}"));
        assert!(p.converged, "alpha = {alpha} did not converge");
        assert!((p.cl - cl).abs() < 6e-4, "alpha = {alpha}: CL {} vs reference {cl}", p.cl);
        assert!((p.cd - cd).abs() < 2e-5 + 2e-3 * cd, "alpha = {alpha}: CD {} vs reference {cd}", p.cd);
        assert!((p.cm - cm).abs() < 6e-4, "alpha = {alpha}: CM {} vs reference {cm}", p.cm);
    }
}

#[test]
fn naca0012_inviscid_lift() {
    let mut analysis = Analysis::new(Airfoil::naca("0012").unwrap()).unwrap();
    let p = analysis.solve_alpha(5.0);
    assert!((p.cl - 0.6033).abs() < 1e-4);
    assert_eq!(p.cd, 0.0);
}

#[test]
fn naca0012_re1e6() {
    let points = sweep("0012", Settings { reynolds: 1e6, ..Settings::default() }, 0.0, 15.0, 1.0);
    assert_matches(
        &points,
        &[
            (0.0, 0.0000, 0.00540, 0.0000),
            (4.0, 0.4278, 0.00728, 0.0060),
            (8.0, 0.9099, 0.01211, -0.0039),
            (12.0, 1.2454, 0.01936, 0.0133),
            (15.0, 1.3807, 0.03169, 0.0310),
        ],
    );
    let p4 = points.iter().find(|p| p.alpha == 4.0).unwrap();
    assert!((p4.xtr_top - 0.2537).abs() < 2e-3 && (p4.xtr_bottom - 0.9685).abs() < 2e-3);
}

#[test]
fn naca0012_re1e5_laminar_bubbles() {
    let points = sweep("0012", Settings { reynolds: 1e5, ..Settings::default() }, 0.0, 6.0, 1.0);
    assert_matches(
        &points,
        &[(2.0, 0.3730, 0.01445, -0.0294), (4.0, 0.5362, 0.01519, -0.0145), (6.0, 0.6861, 0.01958, -0.0008)],
    );
}

#[test]
fn naca23012_re3e6() {
    let points = sweep("23012", Settings { reynolds: 3e6, ..Settings::default() }, 0.0, 12.0, 2.0);
    assert_matches(
        &points,
        &[(0.0, 0.1290, 0.00616, -0.0102), (6.0, 0.7981, 0.00690, -0.0091), (12.0, 1.4342, 0.01298, -0.0082)],
    );
}

#[test]
fn naca2412_compressible() {
    let points = sweep("2412", Settings { reynolds: 6e6, mach: 0.6, ..Settings::default() }, 0.0, 6.0, 2.0);
    assert_matches(
        &points,
        &[(0.0, 0.3131, 0.00573, -0.0656), (4.0, 0.8926, 0.00691, -0.0564), (6.0, 1.1774, 0.00953, -0.0429)],
    );
}

#[test]
fn prescribed_lift_recovers_alpha() {
    let mut analysis = Analysis::new(Airfoil::naca("0012").unwrap()).unwrap();
    analysis.set_settings(Settings { reynolds: 1e6, ..Settings::default() });
    let p = analysis.solve_cl(0.4278);
    assert!(p.converged);
    assert!((p.alpha - 4.0).abs() < 0.01, "alpha = {}", p.alpha);
    assert!((p.cd - 0.00728).abs() < 3e-5);
}
