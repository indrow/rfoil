"""The main window: sidebar of controls and tabs of plots."""

from __future__ import annotations

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction, QCloseEvent, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QScrollArea,
    QSplitter,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .. import __version__
from . import theme
from .panels import AirfoilPanel, AnalysisPanel, ConditionsPanel
from .plots import BoundaryLayerView, GeometryView, PressureView
from .polar_tab import PolarTab
from .session import Session

SIDEBAR_WIDTH = 350


class MainWindow(QMainWindow):
    def __init__(self, session: Session) -> None:
        super().__init__()
        self.session = session
        self.setWindowTitle("rfoil")
        self.resize(1440, 900)

        self.airfoil_panel = AirfoilPanel(session)
        self.conditions_panel = ConditionsPanel(session)
        self.analysis_panel = AnalysisPanel(session)

        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.addTab(PressureView(session), "Pressure")
        self.tabs.addTab(BoundaryLayerView(session), "Boundary layer")
        self.tabs.addTab(PolarTab(session), "Polars")
        self.tabs.addTab(GeometryView(session), "Geometry")
        self.log = QPlainTextEdit(readOnly=True)
        self.log.setMaximumBlockCount(5000)
        self.tabs.addTab(self.log, "Log")

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._sidebar())
        splitter.addWidget(self.tabs)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([SIDEBAR_WIDTH, 1100])
        splitter.setContentsMargins(8, 8, 8, 0)
        self.setCentralWidget(splitter)

        self.status = QLabel()
        self.statusBar().addWidget(self.status, 1)
        self._build_menus()

        session.message.connect(self._on_message)
        session.iteration.connect(self._on_iteration)
        session.solved.connect(self._on_solved)
        self._restore_geometry()
        self._on_message(f"rfoil {__version__}: {session.airfoil.name} loaded")

    def _sidebar(self) -> QWidget:
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(0, 0, 6, 0)
        layout.setSpacing(10)
        for panel in (self.airfoil_panel, self.conditions_panel, self.analysis_panel):
            layout.addWidget(panel)
        layout.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidget(content)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setMinimumWidth(max(SIDEBAR_WIDTH, content.minimumSizeHint().width() + 14))
        return scroll

    def _build_menus(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        self._action(file_menu, "&Open airfoil…", self.airfoil_panel.open_file, QKeySequence.StandardKey.Open)
        self._action(file_menu, "&Save airfoil…", self.airfoil_panel.save_file, QKeySequence.StandardKey.Save)
        file_menu.addSeparator()
        self._action(file_menu, "&Quit", self.close, QKeySequence.StandardKey.Quit)

        analysis = self.menuBar().addMenu("&Analysis")
        self._action(analysis, "Solve point", self.analysis_panel.solve.click, "Ctrl+Return")
        self._action(analysis, "Run sequence", self.analysis_panel.sweep.click, "Ctrl+Shift+Return")
        self._action(analysis, "Stop", self.session.stop, "Esc")
        analysis.addSeparator()
        self._action(analysis, "Re-initialise boundary layer", self.session.reset_boundary_layer)

        view = self.menuBar().addMenu("&View")
        for index in range(self.tabs.count()):
            self._action(
                view,
                self.tabs.tabText(index),
                lambda i=index: self.tabs.setCurrentIndex(i),
                f"Ctrl+{index + 1}",
            )

        help_menu = self.menuBar().addMenu("&Help")
        self._action(help_menu, "&About rfoil", self._about)

    def _action(self, menu, text: str, slot, shortcut=None) -> QAction:
        action = QAction(text, self)
        if shortcut is not None:
            action.setShortcut(QKeySequence(shortcut))
        action.triggered.connect(slot)
        menu.addAction(action)
        return action

    def _on_message(self, text: str) -> None:
        self.status.setText(text)
        self.log.appendPlainText(text)

    def _on_iteration(self, step) -> None:
        self.log.appendPlainText(
            f"  {step.iteration:3d}  rms {step.rms:9.3e}  max {step.max:10.3e}  rlx {step.relaxation:5.3f}"
            f"   α {step.alpha:7.3f}  CL {step.cl:7.4f}  CD {step.cd:8.5f}"
        )

    def _on_solved(self, solution) -> None:
        point = solution.point
        state = "converged" if point.converged else "NOT converged"
        drag = f"  CD {point.cd:.5f}" if point.viscous else ""
        self._on_message(f"α {point.alpha:.3f}°  CL {point.cl:.4f}{drag}  CM {point.cm:.4f}  ({state})")

    def _about(self) -> None:
        QMessageBox.about(
            self,
            "About rfoil",
            f"<h3>rfoil {__version__}</h3>"
            "<p>Viscous/inviscid airfoil analysis with a Rust solver and a Qt interface.</p>"
            "<p>The method follows Drela's viscous/inviscid formulation; see the README "
            "for references.</p>"
            "<p>Licensed under the GNU GPL v2 or later.</p>",
        )

    def _restore_geometry(self) -> None:
        state = QSettings("rfoil", "rfoil").value("window/geometry")
        if state is not None:
            self.restoreGeometry(state)

    def closeEvent(self, event: QCloseEvent) -> None:
        QSettings("rfoil", "rfoil").setValue("window/geometry", self.saveGeometry())
        self.session.shutdown()
        super().closeEvent(event)


def apply_initial_theme(app: QApplication) -> None:
    theme.apply(app, theme.DARK if theme.system_prefers_dark() else theme.LIGHT)
