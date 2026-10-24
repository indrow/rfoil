"""Tests of the Python API (bindings and polar model)."""

import pytest
import rfoil
from rfoil.polar import Polar


@pytest.fixture
def viscous():
    return rfoil.Analysis(rfoil.Airfoil.naca("0012"), settings=rfoil.Settings(reynolds=1e6))


def test_naca_and_parse_round_trip():
    foil = rfoil.Airfoil.naca("2412")
    again = rfoil.Airfoil.parse(foil.to_dat())
    assert again.name == "NACA 2412"
    assert len(again) == len(foil)
    thickness, *_ = again.shape()
    assert thickness == pytest.approx(0.12, abs=2e-3)


def test_invalid_input_raises_value_error():
    with pytest.raises(ValueError):
        rfoil.Airfoil.naca("99")
    with pytest.raises(ValueError):
        rfoil.Analysis(rfoil.Airfoil.naca("0012"), settings=rfoil.Settings(mach=1.2))


def test_viscous_point_matches_reference(viscous):
    p = viscous.solve_alpha(4.0)
    assert p.converged
    assert p.cl == pytest.approx(0.4278, abs=5e-4)
    assert p.cd == pytest.approx(0.00728, abs=3e-5)
    _, bottom = viscous.boundary_layer()
    assert bottom.te < len(bottom.x)  # wake follows the trailing edge
    assert len(viscous.wake[0]) > 0


def test_progress_callback_sees_iterations(viscous):
    seen = []
    p = viscous.solve_alpha(2.0, progress=seen.append)
    assert len(seen) == p.iterations
    assert seen[-1].rms < 1e-4


def test_sweep_callback_can_stop(viscous):
    points = viscous.sweep_alpha(0, 10, 1, callback=lambda p: p.alpha < 3)
    assert [p.alpha for p in points] == [0, 1, 2, 3]


def test_callback_exception_propagates(viscous):
    def boom(_):
        raise KeyError("stop")

    with pytest.raises(KeyError):
        viscous.sweep_alpha(0, 4, 1, callback=boom)


def test_settings_round_trip(viscous):
    s = viscous.settings
    s.ncrit = 5.0
    s.scaling = "lift"
    viscous.settings = s
    assert viscous.settings.ncrit == 5.0
    assert viscous.settings.scaling == "lift"
    s.scaling = "nonsense"
    with pytest.raises(ValueError):
        viscous.settings = s


def test_polar_accumulates_and_exports(viscous):
    polar = Polar.for_settings("NACA 0012", viscous.settings)
    for p in viscous.sweep_alpha(4, 0, 2):
        polar.add(p)
    polar.add(viscous.solve_alpha(2.0))  # replaces the existing α = 2 point
    assert [p.alpha for p in polar.points] == [0, 2, 4]
    text = polar.to_text()
    assert "Calculated polar for: NACA 0012" in text
    assert len(polar.to_csv().splitlines()) == 4
