"""Single points and sequences, plus a read-out of the latest result."""

from __future__ import annotations

from PySide6.QtWidgets import QButtonGroup, QProgressBar, QPushButton, QWidget

from .. import theme
from ..session import Session
from ..widgets import Card, ValueLabel, caption, hbox, spin
from ..worker import Solution


class AnalysisPanel(Card):
    def __init__(self, session: Session, parent: QWidget | None = None) -> None:
        super().__init__("Analysis", parent)
        self.session = session

        self.by_alpha = QPushButton("Angle α", checkable=True, checked=True)
        self.by_cl = QPushButton("Lift CL", checkable=True)
        group = QButtonGroup(self)
        group.addButton(self.by_alpha)
        group.addButton(self.by_cl)
        group.buttonToggled.connect(self._mode_changed)
        self.form.addRow(hbox(self.by_alpha, self.by_cl))

        self.value = self.row("Value", spin(4.0, -30.0, 30.0, 0.5, 2))
        self.solve = QPushButton("Solve point")
        self.solve.setObjectName("Primary")
        self.solve.clicked.connect(self._solve_point)
        self.form.addRow(self.solve)

        self.start = spin(-4.0, -30.0, 30.0, 0.5, 2)
        self.end = spin(14.0, -30.0, 30.0, 0.5, 2)
        self.step = spin(1.0, 0.01, 10.0, 0.25, 2)
        self.form.addRow(hbox(self.start, self.end, self.step))
        self.form.addRow(hbox(*(caption(t) for t in ("start", "end", "step"))))
        self.sweep = QPushButton("Run sequence")
        self.sweep.clicked.connect(self._solve_sequence)
        self.stop = QPushButton("Stop", enabled=False)
        self.stop.clicked.connect(session.stop)
        self.form.addRow(hbox(self.sweep, self.stop))
        self.progress = QProgressBar(textVisible=False)
        self.progress.setRange(0, 1)
        self.form.addRow(self.progress)

        self.readouts = {
            "alpha": self.row("α", ValueLabel("{:.3f}°")),
            "cl": self.row("CL", ValueLabel("{:.4f}")),
            "cd": self.row("CD", ValueLabel("{:.5f}")),
            "cm": self.row("CM", ValueLabel("{:.4f}")),
            "ld": self.row("L/D", ValueLabel("{:.1f}")),
            "xtr_top": self.row("Transition top", ValueLabel("{:.4f}")),
            "xtr_bottom": self.row("Transition bottom", ValueLabel("{:.4f}")),
            "iterations": self.row("Iterations", ValueLabel("{}")),
        }

        session.busy_changed.connect(self._set_busy)
        session.solved.connect(self.show_solution)
        session.progress.connect(self._progress)

    def _mode_changed(self) -> None:
        alpha = self.by_alpha.isChecked()
        lo, hi, step = (-30.0, 30.0, 0.5) if alpha else (-2.0, 3.0, 0.05)
        for box, value in (
            (self.value, 4.0 if alpha else 0.5),
            (self.start, -4.0 if alpha else 0.0),
            (self.end, 14.0 if alpha else 1.2),
        ):
            box.setRange(lo, hi)
            box.setSingleStep(step)
            box.setValue(value)
        self.step.setValue(1.0 if alpha else 0.1)

    def _solve_point(self) -> None:
        run = self.session.run_alpha if self.by_alpha.isChecked() else self.session.run_cl
        run(self.value.value())

    def _solve_sequence(self) -> None:
        run = self.session.run_alpha if self.by_alpha.isChecked() else self.session.run_cl
        run(self.start.value(), self.end.value(), self.step.value())

    def _set_busy(self, busy: bool) -> None:
        for w in (self.solve, self.sweep, self.by_alpha, self.by_cl):
            w.setEnabled(not busy)
        self.stop.setEnabled(busy)
        if busy:
            self.progress.setRange(0, 0)  # indeterminate until the first point lands

    def _progress(self, done: int, total: int) -> None:
        self.progress.setRange(0, total)
        self.progress.setValue(done)

    def show_solution(self, solution: Solution) -> None:
        point = solution.point
        viscous_only = ("cd", "ld", "xtr_top", "xtr_bottom", "iterations")
        for key, label in self.readouts.items():
            shown = point.viscous or key not in viscous_only
            label.show_value(getattr(point, key) if shown else None)
            label.setStyleSheet("" if point.converged else f"color: {theme.current().bad};")
