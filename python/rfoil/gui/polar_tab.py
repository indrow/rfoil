"""Polar plots with a list of polars and a table of the selected one."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QIcon, QPixmap
from PySide6.QtWidgets import (
    QFileDialog,
    QHeaderView,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..polar import COLUMNS
from .plots import PolarPlots, color_of
from .session import Session
from .widgets import caption, hbox

HEADERS = ("α", "CL", "CD", "CDp", "CM", "Xtr top", "Xtr bot")
FORMATS = ("{:.3f}", "{:.4f}", "{:.5f}", "{:.5f}", "{:.4f}", "{:.4f}", "{:.4f}")


class PolarTab(QSplitter):
    def __init__(self, session: Session) -> None:
        super().__init__(Qt.Orientation.Horizontal)
        self.session = session
        self.plots = PolarPlots()
        self.addWidget(self.plots)

        side = QWidget()
        layout = QVBoxLayout(side)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(caption("POLARS  (check to show)"))
        self.list = QListWidget()
        self.list.itemChanged.connect(self._redraw)
        self.list.currentRowChanged.connect(self._fill_table)
        layout.addWidget(self.list, 1)

        export = QPushButton("Export…")
        export.clicked.connect(self._export)
        remove = QPushButton("Remove")
        remove.clicked.connect(self._remove)
        clear = QPushButton("Clear all")
        clear.clicked.connect(session.clear_polars)
        layout.addWidget(hbox(export, remove, clear))

        self.table = QTableWidget(0, len(HEADERS))
        self.table.setHorizontalHeaderLabels(HEADERS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.verticalHeader().setDefaultSectionSize(22)
        layout.addWidget(self.table, 2)
        self.addWidget(side)
        self.setStretchFactor(0, 3)
        self.setStretchFactor(1, 1)
        self.setSizes([820, 480])

        self._current = None
        session.polars_changed.connect(self.refresh)
        session.solved.connect(self._on_solved)

    def refresh(self) -> None:
        """Rebuilds the list from the session, keeping visibility and selection."""
        polars = self.session.polars
        hidden = {
            self.list.item(i).text()
            for i in range(self.list.count())
            if self.list.item(i).checkState() != Qt.CheckState.Checked
        }
        row = self.list.currentRow()
        self.list.blockSignals(True)
        self.list.clear()
        for i, polar in enumerate(polars):
            item = QListWidgetItem(_swatch(color_of(i)), polar.label)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Unchecked if polar.label in hidden else Qt.CheckState.Checked)
            item.setToolTip(f"{len(polar)} points, Ncrit {polar.ncrit:g}")
            self.list.addItem(item)
        self.list.blockSignals(False)
        if polars:
            self.list.setCurrentRow(min(max(row, 0), len(polars) - 1) if row >= 0 else len(polars) - 1)
        self._redraw()
        self._fill_table(self.list.currentRow())

    def _on_solved(self, solution) -> None:
        self._current = solution.point
        self._redraw()

    def _visible(self) -> list[bool]:
        return [self.list.item(i).checkState() == Qt.CheckState.Checked for i in range(self.list.count())]

    def _redraw(self) -> None:
        self.plots.set_polars(self.session.polars, self._visible(), self._current)

    def _fill_table(self, row: int) -> None:
        points = self.session.polars[row].points if 0 <= row < len(self.session.polars) else []
        self.table.setRowCount(len(points))
        for r, p in enumerate(points):
            for c, (attr, fmt) in enumerate(zip(COLUMNS, FORMATS)):
                item = QTableWidgetItem(fmt.format(getattr(p, attr)))
                item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.table.setItem(r, c, item)

    def _export(self) -> None:
        row = self.list.currentRow()
        if row < 0:
            return
        name = self.session.polars[row].label.replace(" ", "_") + ".txt"
        path, _ = QFileDialog.getSaveFileName(
            self, "Export polar", name, "Text polar (*.txt *.pol);;CSV (*.csv)"
        )
        if path:
            self.session.save_polar(row, path)

    def _remove(self) -> None:
        row = self.list.currentRow()
        if row >= 0:
            self.session.remove_polar(row)


def _swatch(color: str) -> QIcon:
    pixmap = QPixmap(12, 12)
    pixmap.fill(QColor(color))
    return QIcon(pixmap)
