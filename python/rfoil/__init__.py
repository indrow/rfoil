"""rfoil: viscous/inviscid airfoil analysis with a Rust solver."""

from ._core import (
    Airfoil,
    Analysis,
    BoundaryLayerSide,
    Iteration,
    OperatingPoint,
    Paneling,
    Settings,
    Surface,
    __version__,
)

__all__ = [
    "Airfoil",
    "Analysis",
    "BoundaryLayerSide",
    "Iteration",
    "OperatingPoint",
    "Paneling",
    "Settings",
    "Surface",
    "__version__",
]
