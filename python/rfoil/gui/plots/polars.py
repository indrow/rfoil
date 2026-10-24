"""Overlaid polars: drag polar, lift and moment curves, transition."""

from __future__ import annotations

from typing import NamedTuple

import pyqtgraph as pg

from ...polar import Polar
from .. import theme


class Pane(NamedTuple):
    """One polar plot: which operating-point attributes go on which axis."""

    row: int
    col: int
    x: str
    y: str
    x_label: str
    y_label: str


PANES = (
    Pane(0, 0, "cd", "cl", "CD", "CL"),
    Pane(0, 1, "alpha", "cl", "α (deg)", "CL"),
    Pane(1, 0, "alpha", "cm", "α (deg)", "CM"),
    Pane(1, 1, "xtr_top", "cl", "transition x/c", "CL"),
)


class PolarPlots(pg.GraphicsLayoutWidget):
    def __init__(self) -> None:
        super().__init__()
        self.ci.layout.setContentsMargins(4, 8, 12, 4)
        self.plots: list[pg.PlotItem] = []
        for pane in PANES:
            plot = self.addPlot(row=pane.row, col=pane.col)
            theme.style_plot(plot, pane.x_label, pane.y_label)
            self.plots.append(plot)

        # CL is the shared ordinate of the drag polar, lift curve and transition plot
        cl_plots = [plot for plot, pane in zip(self.plots, PANES) if pane.y == "cl"]
        for plot in cl_plots[1:]:
            plot.setYLink(cl_plots[0])
        lift_curve, moment_curve = self.plots[1], self.plots[2]
        moment_curve.setXLink(lift_curve)

        # ring around the most recent operating point
        self.markers = [
            plot.plot(pen=None, symbol="o", symbolSize=11, symbolBrush=None) for plot in self.plots
        ]

    def set_polars(self, polars: list[Polar], visible: list[bool], current=None) -> None:
        for plot, pane in zip(self.plots, PANES):
            for item in list(plot.listDataItems()):
                if item not in self.markers:
                    plot.removeItem(item)
            for index, (polar, shown) in enumerate(zip(polars, visible)):
                if shown and polar.points:
                    self._draw(plot, pane, polar, color_of(index))
        self._mark(current)

    @staticmethod
    def _draw(plot: pg.PlotItem, pane: Pane, polar: Polar, color: str) -> None:
        y = polar.column(pane.y)
        plot.plot(
            polar.column(pane.x),
            y,
            pen=theme.pen(color, 1.8),
            symbol="o",
            symbolSize=4,
            symbolBrush=color,
            symbolPen=None,
        )
        if pane.x == "xtr_top":
            # bottom-surface transition dashed in the same colour
            plot.plot(polar.column("xtr_bottom"), y, pen=theme.pen(color, 1.4, dashed=True))

    def _mark(self, point) -> None:
        show = point is not None and point.converged and point.viscous
        for marker, pane in zip(self.markers, PANES):
            if show:
                marker.setData([getattr(point, pane.x)], [getattr(point, pane.y)])
                marker.setSymbolPen(pg.mkPen(theme.current().text, width=1.5))
            else:
                marker.clear()


def color_of(index: int) -> str:
    series = theme.current().series
    return series[index % len(series)]
