"""Pressure distribution with the airfoil and boundary-layer edge beneath it."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg

from .. import theme
from ..session import Session
from ..worker import Solution
from .common import displacement_surface, legend, split_at_leading_edge


class PressureView(pg.GraphicsLayoutWidget):
    def __init__(self, session: Session) -> None:
        super().__init__()
        self.session = session
        t = theme.current()
        self.ci.setSpacing(0)
        self.ci.layout.setContentsMargins(4, 8, 12, 4)

        self.cp = self.addPlot(row=0, col=0)
        theme.style_plot(self.cp, "", "Cp")
        self.cp.invertY(True)
        self.cp.getAxis("bottom").setStyle(showValues=False)
        self.cp.addLine(y=0, pen=pg.mkPen(t.border))
        self.sonic = self.cp.addLine(y=0, pen=theme.pen(t.bad, 1, dashed=True))
        self.sonic.hide()
        self.cp_legend = legend(self.cp, corner="bottom-right")
        self.inviscid = self.cp.plot(pen=theme.pen(t.inviscid, 1.3, dashed=True), name="inviscid")
        self.upper = self.cp.plot(pen=theme.pen(t.upper, 2.0), name="upper")
        self.lower = self.cp.plot(pen=theme.pen(t.lower, 2.0), name="lower")
        self.wake_cp = self.cp.plot(pen=theme.pen(t.lower, 1.2, dashed=True))
        self.summary = pg.TextItem(anchor=(1, 0), color=t.text)
        self.cp.addItem(self.summary, ignoreBounds=True)

        self.geom = self.addPlot(row=1, col=0)
        theme.style_plot(self.geom, "x/c", "")
        self.geom.getAxis("left").setStyle(showValues=False)
        self.geom.setXLink(self.cp)
        self.geom.setMouseEnabled(x=True, y=False)
        self.geom.showGrid(x=False, y=False)
        self.ci.layout.setRowStretchFactor(0, 5)
        self.ci.layout.setRowStretchFactor(1, 2)
        self._y_center = 0.0
        self.outline = self.geom.plot(pen=theme.pen(t.text, 1.4), fillLevel=None)
        self.dstar_upper = self.geom.plot(pen=theme.pen(t.upper, 1.0, dashed=True))
        self.dstar_lower = self.geom.plot(pen=theme.pen(t.lower, 1.0, dashed=True))
        self.wake_line = self.geom.plot(pen=theme.pen(t.inviscid, 1.0))
        self.transition = self.geom.plot(
            pen=None, symbol="o", symbolSize=7, symbolBrush=t.good, symbolPen=None
        )

        self.cp.sigRangeChanged.connect(self._place_summary)
        self.geom.getViewBox().sigXRangeChanged.connect(self._true_aspect)
        self.geom.getViewBox().sigResized.connect(self._true_aspect)
        session.airfoil_changed.connect(self.show_geometry)
        session.solved.connect(self.show_solution)
        self.show_geometry()

    def show_geometry(self) -> None:
        x, y = self.session.panels()
        self._y_center = 0.5 * (min(y) + max(y))
        self.outline.setData(np.r_[x, x[:1]], np.r_[y, y[:1]])
        for item in (
            self.inviscid,
            self.upper,
            self.lower,
            self.wake_cp,
            self.dstar_upper,
            self.dstar_lower,
            self.wake_line,
            self.transition,
        ):
            item.clear()
        self.summary.setText("")
        self.sonic.hide()
        self.cp.setXRange(-0.05, 1.15, padding=0)

    def show_solution(self, solution: Solution) -> None:
        surface, point = solution.surface, solution.point
        n = surface.airfoil_nodes
        x = np.asarray(surface.x)
        cp_viscous, cp_inviscid = np.asarray(surface.cp_viscous), np.asarray(surface.cp_inviscid)
        upper, lower = split_at_leading_edge(x[:n])

        if point.viscous:
            self.inviscid.setData(x[:n], cp_inviscid[:n])
        else:
            self.inviscid.clear()
        self.upper.setData(x[upper], cp_viscous[upper])
        self.lower.setData(x[lower], cp_viscous[lower])
        if point.viscous and len(x) > n:
            self.wake_cp.setData(x[n - 1 :], cp_viscous[n - 1 :])
        else:
            self.wake_cp.clear()

        self._show_boundary_layer(solution)
        self._show_sonic_line(point.mach)

        text = f"α = {point.alpha:.2f}°    CL = {point.cl:.4f}    CM = {point.cm:.4f}"
        if point.viscous:
            text += f"\nCD = {point.cd:.5f}    L/D = {point.ld:.1f}    Re = {point.reynolds:,.0f}"
        if not point.converged:
            text += "\n(not converged)"
        colors = theme.current()
        self.summary.setText(text, color=colors.text if point.converged else colors.bad)
        self._place_summary()

    def _show_boundary_layer(self, solution: Solution) -> None:
        if solution.boundary_layer is None:
            for item in (self.dstar_upper, self.dstar_lower, self.wake_line, self.transition):
                item.clear()
            return
        top, bottom = solution.boundary_layer
        self.dstar_upper.setData(*displacement_surface(top, outward=+1))
        self.dstar_lower.setData(*displacement_surface(bottom, outward=-1))
        self.wake_line.setData(*solution.wake)

        # transition dots on the contour (the upper surface runs from TE to LE)
        point = solution.point
        x, y = (np.asarray(v) for v in self.session.panels())
        upper, lower = split_at_leading_edge(x)
        dots = [
            (xtr, np.interp(xtr, x[nodes][::order], y[nodes][::order]))
            for xtr, nodes, order in ((point.xtr_top, upper, -1), (point.xtr_bottom, lower, 1))
            if xtr < 0.999
        ]
        self.transition.setData([d[0] for d in dots], [d[1] for d in dots])

    def _show_sonic_line(self, mach: float) -> None:
        if mach <= 0:
            self.sonic.hide()
            return
        self.sonic.setValue(sonic_cp(mach))
        self.sonic.show()

    def _true_aspect(self, *_) -> None:
        """Keeps the airfoil strip at 1:1 scale by deriving its y-range from x."""
        vb = self.geom.getViewBox()
        (x0, x1), rect = vb.viewRange()[0], vb.sceneBoundingRect()
        if rect.width() > 0:
            half = 0.5 * (x1 - x0) * rect.height() / rect.width()
            vb.setYRange(self._y_center - half, self._y_center + half, padding=0)

    def _place_summary(self) -> None:
        (_, x1), (y0, _) = self.cp.viewRange()
        # the y axis is inverted: the visual top is the minimum
        self.summary.setPos(x1, y0)


def sonic_cp(mach: float, gamma: float = 1.4) -> float:
    """Pressure coefficient at which the local flow becomes sonic."""
    stagnation = (1 + 0.5 * (gamma - 1) * mach**2) / (1 + 0.5 * (gamma - 1))
    return 2 / (gamma * mach**2) * (stagnation ** (gamma / (gamma - 1)) - 1)
