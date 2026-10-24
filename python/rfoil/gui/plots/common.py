"""Helpers shared by the plot views."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg

from ... import BoundaryLayerSide
from .. import theme


def legend(plot: pg.PlotItem, corner: str = "top-right") -> pg.LegendItem:
    """A compact legend in a corner of ``plot``."""
    t = theme.current()
    offset = {"top-right": (-10, 10), "bottom-right": (-10, -10)}[corner]
    item = plot.addLegend(offset=offset, labelTextColor=t.muted, labelTextSize="8pt")
    item.setBrush(pg.mkBrush(t.surface + "cc"))
    item.setPen(pg.mkPen(t.border))
    return item


def split_at_leading_edge(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Index arrays of the upper (TE to LE) and lower (LE to TE) surface nodes."""
    ile = int(np.argmin(x))
    return np.arange(0, ile + 1), np.arange(ile, len(x))


def displacement_surface(side: BoundaryLayerSide, outward: int) -> tuple[np.ndarray, np.ndarray]:
    """Surface offset by δ* along the outward normal, up to the trailing edge.

    ``outward`` is +1 for the upper side and −1 for the lower side (both
    sides run from the stagnation point towards the trailing edge).
    """
    end = side.te + 1
    x, y, d = (np.asarray(v[:end]) for v in (side.x, side.y, side.dstar))
    if len(x) < 2:
        return x, y
    tx, ty = np.gradient(x), np.gradient(y)
    norm = np.hypot(tx, ty)
    norm[norm == 0] = 1.0
    nx, ny = -outward * ty / norm, outward * tx / norm
    return x + d * nx, y + d * ny


def without_zeros(values: list[float]) -> np.ndarray:
    """Replaces exact zeros (quantity undefined there) by NaN so lines break."""
    a = np.asarray(values, dtype=float)
    a[a == 0.0] = np.nan
    return a
