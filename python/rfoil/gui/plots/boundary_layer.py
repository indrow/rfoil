"""Boundary-layer quantities along both surfaces and the wake."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pyqtgraph as pg
from PySide6.QtWidgets import QComboBox, QLabel, QVBoxLayout, QWidget

from .. import theme
from ..session import Session
from ..widgets import hbox
from ..worker import Solution
from .common import legend, without_zeros


@dataclass(frozen=True)
class Quantity:
    label: str
    attr: str
    axis: str
    #: whether the quantity continues into the wake
    in_wake: bool = True
    #: zeros mean "not defined here" (e.g. N on turbulent stations)
    sparse: bool = False


QUANTITIES = (
    Quantity("Edge velocity", "ue", "Ue / V∞"),
    Quantity("Displacement thickness", "dstar", "δ* / c"),
    Quantity("Momentum thickness", "theta", "θ / c"),
    Quantity("Shape factor", "h", "H = δ*/θ"),
    Quantity("Skin friction", "cf", "Cf", in_wake=False),
    Quantity("Amplification", "n", "N", in_wake=False, sparse=True),
    Quantity("Shear stress", "ctau", "√Cτ", sparse=True),
)


class BoundaryLayerView(QWidget):
    def __init__(self, session: Session) -> None:
        super().__init__()
        t = theme.current()
        self.solution: Solution | None = None

        self.selector = QComboBox()
        for q in QUANTITIES:
            self.selector.addItem(q.label, q)
        self.selector.currentIndexChanged.connect(self._redraw)
        self.hint = QLabel("Run a viscous analysis to see boundary-layer data")
        self.hint.setObjectName("Muted")

        self.canvas = pg.PlotWidget()
        self.plot = self.canvas.getPlotItem()
        legend(self.plot)
        self.upper = self.plot.plot(pen=theme.pen(t.upper, 2.0), name="upper")
        self.lower = self.plot.plot(pen=theme.pen(t.lower, 2.0), name="lower")
        self.wake = self.plot.plot(pen=theme.pen(t.lower, 1.4, dashed=True), name="wake")
        self.te = self.plot.addLine(x=1.0, pen=theme.pen(t.border, 1.0, dashed=True))

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(hbox(self.selector, self.hint))
        layout.addWidget(self.canvas, 1)

        session.solved.connect(self.show_solution)
        session.airfoil_changed.connect(self.clear)
        self._redraw()

    def clear(self) -> None:
        self.solution = None
        self._redraw()

    def show_solution(self, solution: Solution) -> None:
        self.solution = solution
        self._redraw()

    def _redraw(self) -> None:
        quantity: Quantity = self.selector.currentData()
        theme.style_plot(self.plot, "x/c", quantity.axis)
        layer = self.solution.boundary_layer if self.solution else None
        self.hint.setVisible(layer is None)
        if layer is None:
            for curve in (self.upper, self.lower, self.wake):
                curve.clear()
            return

        top, bottom = layer

        def values(side) -> np.ndarray:
            raw = getattr(side, quantity.attr)
            return without_zeros(raw) if quantity.sparse else np.asarray(raw, dtype=float)

        # both sides run from the stagnation point; the bottom continues into the wake
        end_top, end_bottom = top.te + 1, bottom.te + 1
        self.upper.setData(np.asarray(top.x[:end_top]), values(top)[:end_top])
        x_bottom, y_bottom = np.asarray(bottom.x), values(bottom)
        self.lower.setData(x_bottom[:end_bottom], y_bottom[:end_bottom])
        if quantity.in_wake:
            self.wake.setData(x_bottom[bottom.te :], y_bottom[bottom.te :])
        else:
            self.wake.clear()
        self.te.setValue(max(top.x[top.te], bottom.x[bottom.te]))
        self.plot.enableAutoRange()
