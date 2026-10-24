"""The paneled airfoil, with nodes, at true aspect ratio."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg

from .. import theme
from ..session import Session
from ..worker import Solution


class GeometryView(pg.PlotWidget):
    def __init__(self, session: Session) -> None:
        super().__init__()
        t = theme.current()
        self.session = session
        self.plot = self.getPlotItem()
        theme.style_plot(self.plot, "x/c", "y/c")
        self.plot.setAspectLocked(True)
        self.outline = self.plot.plot(pen=theme.pen(t.text, 1.3))
        self.nodes = self.plot.plot(
            pen=None, symbol="o", symbolSize=3.5, symbolBrush=t.accent, symbolPen=None
        )
        self.wake = self.plot.plot(
            pen=theme.pen(t.inviscid, 1.0), symbol="o", symbolSize=3, symbolBrush=t.inviscid, symbolPen=None
        )
        self.camber = self.plot.plot(pen=theme.pen(t.lower, 1.0, dashed=True))
        session.airfoil_changed.connect(self.refresh)
        session.solved.connect(self._show_wake)
        self.refresh()

    def refresh(self) -> None:
        x, y = (np.asarray(v) for v in self.session.panels())
        self.outline.setData(np.r_[x, x[:1]], np.r_[y, y[:1]])
        self.nodes.setData(x, y)
        self.camber.setData(*mean_line(x, y))
        self.wake.clear()
        self.plot.autoRange()

    def _show_wake(self, solution: Solution) -> None:
        self.wake.setData(*solution.wake)


def mean_line(x: np.ndarray, y: np.ndarray, samples: int = 120) -> tuple[np.ndarray, np.ndarray]:
    """Approximate camber line: midpoint of upper and lower surfaces at equal x."""
    ile = int(np.argmin(x))
    xu, yu = x[: ile + 1][::-1], y[: ile + 1][::-1]
    xl, yl = x[ile:], y[ile:]
    xs = 0.5 * (1 - np.cos(np.linspace(0, np.pi, samples))) * (min(xu[-1], xl[-1]) - x[ile]) + x[ile]
    return xs, 0.5 * (np.interp(xs, xu, yu) + np.interp(xs, xl, yl))
