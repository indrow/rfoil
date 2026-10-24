# rfoil

Viscous/inviscid airfoil analysis with a Rust solver and a Python/Qt desktop
interface. rfoil computes pressure distributions, boundary-layer development,
transition, lift, drag and moment for 2D airfoils at subsonic speeds.

The numerical method is the one introduced by Mark Drela in XFOIL (see
[References](#references)). The solver is a fresh implementation of it, and
its results agree with XFOIL 6.97 to every printed digit on the validation
cases below.

## Features

- **Inviscid flow:** linear-vorticity panel method with a Karman-Tsien
  compressibility correction.
- **Viscous flow:** two-equation integral boundary layer with lagged
  dissipation, a wake, and free or forced transition by the e^N method.
  The boundary layer and the inviscid flow are solved together by a global
  Newton method.
- **Operating points:** prescribed angle of attack or prescribed lift
  coefficient. Mach and Reynolds numbers can be fixed or scale with CL.
- **Geometry:** NACA 4- and 5-digit generator, and coordinate files in plain,
  labelled, ISES and Lednicer layouts, with curvature-based repaneling.
- **Desktop interface:** pressure, boundary-layer, polar and geometry views.
  Single points and sequences run in a background thread, and polars can be
  exported as text or CSV.
- **Libraries:** the solver can be used directly from Python or from Rust.

## Installation

Requirements:

- Rust, stable toolchain (edition 2021)
- Python 3.10 or newer
- [uv](https://docs.astral.sh/uv/)

From a clone of this repository:

```sh
uv sync           # builds the Rust extension and installs the Python dependencies
uv run rfoil      # starts the desktop application
```

`uv run python -m rfoil` also starts the application, and
`uv run rfoil path/to/airfoil.dat` opens a coordinate file on start.

To build a wheel instead, run `uv build`; the wheel is written to `dist/`.

## Using the application

The left sidebar holds three panels:

- **Airfoil:** generate a NACA section, open or save a coordinate file, and
  repanel.
- **Conditions:** Reynolds number (or inviscid), Mach number, Mach/Reynolds
  scaling with CL, critical amplification Ncrit, trip locations and the
  iteration limit.
- **Analysis:** solve a single point or a sequence at prescribed angle of
  attack or lift coefficient, and read out the latest result.

The tabs on the right show:

| Tab | Content |
| --- | --- |
| Pressure | Cp distribution, the airfoil with its displacement surface, wake and transition points |
| Boundary layer | Ue, displacement and momentum thickness, H, Cf, N and shear stress along both sides and the wake |
| Polars | Drag polar, lift curve, moment curve and transition locations of every accumulated polar, plus a table |
| Geometry | Panel nodes, mean line and wake |
| Log | Newton iteration history and messages |

Every converged point is added to the polar that matches the current airfoil
and conditions. Keyboard shortcuts: `Ctrl+Return` solves a point,
`Ctrl+Shift+Return` runs a sequence, `Esc` stops it, and `Ctrl+1` to `Ctrl+5`
switch tabs.

## Using the library

### Python

```python
import rfoil

analysis = rfoil.Analysis(
    rfoil.Airfoil.naca("4412"),
    settings=rfoil.Settings(reynolds=1e6, mach=0.0, ncrit=9.0),
)

for point in analysis.sweep_alpha(-4, 14, 1):
    print(f"{point.alpha:5.1f}  CL {point.cl:.4f}  CD {point.cd:.5f}  CM {point.cm:.4f}")

point = analysis.solve_cl(0.8)              # prescribed lift
top, bottom = analysis.boundary_layer()    # x, xi, ue, dstar, theta, h, cf, n, ctau
surface = analysis.surface()               # x, y, cp_viscous, cp_inviscid
```

Consecutive operating points start from the previous boundary layer, so
sweeps in small steps converge quickly. `analysis.reset_boundary_layer()`
starts the next point from scratch. Long solves release the GIL, so they can
run in a worker thread. Each call accepts an optional callback that sees
every Newton iteration (`solve_alpha`, `solve_cl`) or every point
(`sweep_alpha`, `sweep_cl`, where returning `False` stops the sweep).

`rfoil.polar.Polar` collects points and writes them as a text table or as
CSV.

### Rust

The solver is the `rfoil-core` crate in `crates/rfoil-core`. It is not on
crates.io yet, so depend on it by path or git:

```rust
use rfoil_core::{Airfoil, Analysis, Settings};

let mut analysis = Analysis::new(Airfoil::naca("2412")?)?;
analysis.set_settings(Settings { reynolds: 1.0e6, ..Settings::default() });

let point = analysis.solve_alpha(4.0);
println!("CL = {:.4}  CD = {:.5}", point.cl, point.cd);
```

A small command-line example prints a polar:

```sh
cargo run --release --example polar -- 0012 1e6 -4 16 1
# arguments: <naca|file.dat> <Re> <alpha start> <alpha end> <step> [mach] [ncrit]
```

## Project layout

```
crates/rfoil-core/          solver library (pure Rust)
  src/numeric/              splines, dense linear algebra, dual numbers
  src/geometry/             coordinate files, NACA sections, contour, paneling
  src/inviscid/             panels, influence kernels, panel system, wake, forces
  src/viscous/
    closure.rs              empirical correlations (H*, Cf, dissipation, amplification)
    station.rs              one boundary-layer station as dual numbers
    interval.rs             the three interval equations and their linearisation
    transition.rs           e^N amplification and transition location
    layer.rs                stations, sides and the stagnation point
    march.rs                station-by-station marching
    newton/                 global Newton system: assembly, block solve, update
  src/analysis/             public API: Analysis, Settings, results
  tests/reference_polars.rs regression against the reference polars
  examples/polar.rs         command-line polar
crates/rfoil-py/            Python bindings (module rfoil._core)
python/rfoil/               Python package
  polar.py                  polar accumulation and export
  gui/                      Qt application: session, worker thread, panels, plots
python/tests/               Python tests
```

## Design notes

- **Exact Jacobians from dual numbers.** Every boundary-layer quantity is a
  forward-mode dual number over the 14 variables of one interval: both
  stations' primary variables, the free-stream Mach number and Reynolds
  number, the trip location and a trial amplification. The residuals are
  written once as ordinary formulas and the Newton Jacobian follows
  automatically, where XFOIL derives each partial derivative by hand.
- **Typed state.** Stations, sides, the stagnation point and the Newton
  system are explicit types instead of shared global arrays, and cached
  solution state is held in `Option`s.
- **Threading.** The application thread only reads an immutable snapshot of
  the session; the worker thread owns the analysis while a job runs.

## Validation

The reference values come from XFOIL 6.97 built in double precision, run with
160 panel nodes, Ncrit = 9 and free transition. Across the NACA 0012, 2412,
4412 and 23012 sections and the Eppler E387, at Reynolds numbers from 1e5 to
6e6 and Mach numbers from 0 to 0.6, rfoil reproduces CL, CD, CM and the
transition locations to every printed digit. This holds through laminar
separation bubbles and past maximum lift.

A sample, NACA 0012 at Re = 1e6:

| alpha | CL | CD | CM |
| ---: | ---: | ---: | ---: |
| 0 | 0.0000 | 0.00540 | 0.0000 |
| 4 | 0.4278 | 0.00728 | 0.0060 |
| 8 | 0.9099 | 0.01211 | -0.0039 |
| 12 | 1.2454 | 0.01936 | 0.0133 |
| 15 | 1.3807 | 0.03169 | 0.0310 |

These and further cases are checked by `crates/rfoil-core/tests/reference_polars.rs`.
A 21-point viscous polar takes about half a second with a release build.

### Differences from XFOIL 6.97

- Coordinate files are always repaneled. XFOIL uses the raw points until the
  `PANE` command is given; with `PANE` the two agree.
- The derivative of the laminar wake dissipation with respect to Hk has a sign
  error in XFOIL (the first term of `RCD_HK` in `DILW`). It is corrected here.
  It affects only the Newton Jacobian, not converged results.
- Near maximum lift the viscous solution can be non-unique. A large jump in
  angle of attack can settle on a different branch than XFOIL does; in small
  steps the two agree.
- Inverse design (`QDES`, `MDES`), geometry editing (`GDES`, including flaps)
  and plot hardcopy are not implemented.

## References

1. M. Drela, "XFOIL: An Analysis and Design System for Low Reynolds Number
   Airfoils", in T. J. Mueller (ed.), *Low Reynolds Number Aerodynamics*,
   Lecture Notes in Engineering, vol. 54, Springer, 1989, pp. 1-12.
2. M. Drela and M. B. Giles, "Viscous-Inviscid Analysis of Transonic and Low
   Reynolds Number Airfoils", *AIAA Journal*, vol. 25, no. 10, 1987,
   pp. 1347-1355.
3. M. Drela and H. Youngren, XFOIL 6.97 source code and documentation,
   <https://web.mit.edu/drela/Public/web/xfoil/>.

### Correspondence with XFOIL 6.97

For readers comparing the two codes, the table lists the XFOIL routines that
each rfoil module corresponds to. The rfoil code is organised differently and
does not reuse XFOIL's names.

| rfoil | XFOIL 6.97 |
| --- | --- |
| `numeric/spline.rs` | `SPLIND`, `SEGSPL`, `SEVAL`, `DEVAL`, `D2VAL`, `CURV`, `SINVRT` |
| `numeric/linalg.rs` | `TRISOL`, `GAUSS`, `LUDCMP`, `BAKSUB` |
| `geometry/naca.rs` | `NACA4`, `NACA5` |
| `geometry/curve.rs` | `LEFIND`, `SOPPS` |
| `geometry/paneling.rs` | `PANGEN` |
| `inviscid/influence.rs` | `PSILIN`, `PSWLIN` |
| `inviscid/system.rs` | `GGCALC`, `QDCALC` |
| `inviscid/wake.rs` | `XYWAKE`, `QWCALC`, `SETEXP` |
| `inviscid/forces.rs` | `COMSET`, `CPCALC`, `CLCALC` |
| `viscous/closure.rs` | `HKIN`, `HCT`, `HSL`, `HST`, `CFL`, `CFT`, `DIL`, `DILW`, `DAMPL`, `DAMPL2` |
| `viscous/params.rs` | `BLPINI`, and the `TKBL`, `RSTBL`, `HSTINV`, `REYBL` setup |
| `viscous/station.rs` | `BLKIN`, `BLVAR` |
| `viscous/interval.rs` | `BLSYS`, `BLDIF`, `TESYS`, `TRDIF` |
| `viscous/transition.rs` | `AXSET`, `TRCHEK2` |
| `viscous/layer.rs` | `STFIND`, `IBLPAN`, `IBLSYS`, `XICALC`, `UESET`, `STMOVE`, `XIFSET`, `DSLIM` |
| `viscous/march.rs` | `MRCHUE`, `MRCHDU` |
| `viscous/newton/assemble.rs` | `SETBL` |
| `viscous/newton/solve.rs` | `BLSOLV` |
| `viscous/newton/update.rs` | `UPDATE` |
| `analysis/flow.rs` | `SPECAL`, `SPECCL`, `VISCAL`, `CDCALC` |
| `analysis/settings.rs` | `MRCL` |

Settings and their XFOIL counterparts:

| rfoil | XFOIL 6.97 |
| --- | --- |
| `Settings::scaling` (`Fixed`, `FixedLift`, `FixedLiftAndSpeed`) | polar type 1, 2, 3 |
| `Settings::ncrit`, `xtr_top`, `xtr_bottom` | `VPAR`: `N`, `XTR` |
| `Settings::vaccel` | `VPAR`: `VACC` |
| `Settings::max_iterations` | `ITER` |
| `Settings::wake_length` | `WAKLEN` |
| `Amplification::Envelope`, `Modified` | `DAMPL`, `DAMPL2` |
| `PanelingOptions` | `PPAR` menu |
| `Analysis::reset_boundary_layer` | `INIT` |

Boundary-layer calibration constants (`viscous::Calibration`):

| rfoil | XFOIL 6.97 | Default |
| --- | --- | ---: |
| `shear_lag` | `SCCON` | 5.6 |
| `g_beta_a` | `GACON` | 6.70 |
| `g_beta_b` | `GBCON` | 0.75 |
| `g_beta_c` | `GCCON` | 18.0 |
| `wake_dissipation_length` | `DLCON` | 0.9 |
| `transition_ctau_scale` | `CTRCON` | 1.8 |
| `transition_ctau_exponent` | `CTRCEX` | 3.3 |
| `equilibrium_dudx_weight` | `DUXCON` | 1.0 |
| `cf_scale` | `CFFAC` | 1.0 |
| `ctau_eq_scale()` | `CTCON` | derived |

File formats: coordinates are written as a name line followed by `x y` pairs,
the labelled format XFOIL reads and writes. Polars exported as text follow the
layout of XFOIL's polar save files: a header block with the airfoil name and
conditions, then the columns alpha, CL, CD, CDp, CM, Top_Xtr and Bot_Xtr.

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for the
development setup, coding conventions and how numerical changes are verified.

## License

rfoil is free software, licensed under the GNU General Public License,
version 2 or (at your option) any later version. See [LICENSE](LICENSE).

Copyright (C) 2026 Indra Wibisono.

The method and the empirical correlations derive from XFOIL by Mark Drela
and Harold Youngren, which is also distributed under the GNU General Public
License.
