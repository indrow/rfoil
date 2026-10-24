//! Prints a polar: `polar <naca|file.dat> <Re> <alpha0> <alpha1> <dalpha> [mach] [ncrit]`
use rfoil_core::{Airfoil, Analysis, Settings};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let arg = |i: usize, default: f64| args.get(i).and_then(|s| s.parse().ok()).unwrap_or(default);
    let foil = match args.first().map(String::as_str) {
        Some(path) if path.ends_with(".dat") => Airfoil::parse(&std::fs::read_to_string(path).unwrap()).unwrap(),
        Some(code) => Airfoil::naca(code).unwrap(),
        None => Airfoil::naca("0012").unwrap(),
    };
    let mut analysis = Analysis::new(foil).unwrap();
    analysis.set_settings(Settings {
        reynolds: arg(1, 1e6),
        mach: arg(5, 0.0),
        ncrit: arg(6, 9.0),
        max_iterations: 100,
        ..Settings::default()
    });
    println!("  alpha      CL        CD       CDp       CM     Top_Xtr  Bot_Xtr  conv");
    analysis.sweep_alpha(arg(2, 0.0), arg(3, 10.0), arg(4, 1.0), |p| {
        println!(
            "{:7.3} {:8.4} {:9.5} {:9.5} {:8.4} {:8.4} {:8.4}  {}",
            p.alpha, p.cl, p.cd, p.cdp, p.cm, p.xtr_top, p.xtr_bottom, p.converged
        );
        true
    });
}
