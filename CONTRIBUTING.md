# Contributing to rfoil

Thank you for considering a contribution. Bug reports, validation cases,
documentation fixes and code are all welcome.

## Reporting problems

Please include:

- the airfoil (NACA designation or the coordinate file),
- the settings: Reynolds number, Mach number, Ncrit, trips, scaling,
- the angle of attack or lift coefficient, and whether it was a single point or
  part of a sequence (the boundary layer of one point is the starting guess for
  the next, so the order matters),
- what you expected, ideally with reference numbers from XFOIL or experiment.

## Development setup

You need a stable Rust toolchain, Python 3.10 or newer, and
[uv](https://docs.astral.sh/uv/).

```sh
uv sync                                   # build the extension, install dependencies
uv run rfoil                              # run the application
uv sync --reinstall-package rfoil         # rebuild after changing Rust code
```

## Checks

Run these before opening a pull request. All of them must pass.

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --release --workspace
uvx ruff check python
uvx ruff format --check python
uv run pytest
```

Use `--release` for the Rust tests: the regression tests solve full viscous
polars and are slow in debug builds.

## Code conventions

- **Small, focused files.** Each module does one job. Prefer adding a file
  over growing a long one.
- **Descriptive names.** Use names that say what a quantity is (`theta`,
  `dstar`, `re_theta`, `wall_shear`), not abbreviations from other codes.
  Standard symbols such as `hk`, `cf` or `msq` are fine when documented.
- **Write the method, not a translation.** Express the equations with the
  project's own types: dual numbers for derivatives, typed state instead of
  global arrays, functions with a single purpose.
- **Comments explain the physics or the reasoning.** Comments and doc comments
  should stand on their own. They do not refer to routine or variable names of
  other programs; the correspondence with XFOIL is documented in the
  README's References section, so update that table if you add or move a
  module.
- **Formatting** follows `rustfmt.toml` (120 columns) and the ruff settings
  in `pyproject.toml` (110 columns).
- **Python GUI.** Widgets talk to `gui/session.py` and listen to its signals.
  Only the worker thread touches the `Analysis` while a job runs; the GUI
  thread reads the session's snapshot.

## Changes that affect numerical results

rfoil's results match XFOIL 6.97 to every printed digit on the validation
cases, and changes should keep it that way unless the change is a deliberate,
documented improvement.

- **Refactoring** (renaming, restructuring) must leave the results
  bit-identical. Keep the arithmetic and its order of operations unchanged; for
  example, `0.5 * m * m` and `0.5 * (m * m)` round differently. A practical
  check is to print a few polars and boundary-layer sums with `{:?}` (full
  precision) before and after the change and compare the outputs with `diff`.
- **Deliberate numerical changes** need a justification in the pull request,
  updated reference values if the regression tests change, and an entry under
  "Differences from XFOIL 6.97" in the README.
- **New validation cases** go in `crates/rfoil-core/tests/reference_polars.rs`.
  State where the reference values come from.

### Building a reference XFOIL

To produce reference polars, build XFOIL 6.97 from its source in double
precision without the plot library:

1. Compile the Fortran sources with
   `gfortran -O2 -fdefault-real-8 -std=legacy -fno-range-check`.
   (`-fno-range-check` is needed for an integer constant in `xoper.f`.)
2. Instead of linking the X11 plot library, link a file of empty subroutines
   with the names of the plot routines the linker reports as missing.
3. Run it non-interactively with plotting disabled (`PLOP`, then `G`) and the
   commands piped on standard input.
4. When loading a coordinate file, issue `PANE` after `LOAD` so that XFOIL
   repanels it as rfoil does.

Use 160 panel nodes and the default settings unless the case needs others.

## Pull requests

- Keep each pull request to one topic, and describe what changed and why.
- Add or update tests for behaviour changes.
- Update the README when user-visible behaviour, settings or the module
  layout change.

## License

rfoil is licensed under the GNU General Public License, version 2 or later.
By contributing, you agree that your contributions are licensed under the
same terms.
