"""Flow conditions and viscous solver parameters."""

from __future__ import annotations

from PySide6.QtWidgets import QCheckBox, QComboBox, QPushButton, QWidget

from ... import Settings
from ..session import Session
from ..widgets import Card, NumberEdit, int_spin, spin

SCALINGS = {
    "fixed": "Fixed M, Re  (type 1)",
    "lift": "Fixed M√CL, Re√CL  (type 2)",
    "lift_speed": "Fixed M, Re·CL  (type 3)",
}


class ConditionsPanel(Card):
    def __init__(self, session: Session, parent: QWidget | None = None) -> None:
        super().__init__("Conditions", parent)
        self.session = session
        s = session.settings

        self.viscous = QCheckBox("Viscous")
        self.viscous.setChecked(s.viscous)
        self.form.addRow(self.viscous)
        self.reynolds = self.row(
            "Reynolds", NumberEdit(s.reynolds if s.viscous else 1e6, 1e3, 1e10), "Chord Reynolds number"
        )
        self.mach = self.row("Mach", spin(s.mach, 0.0, 0.95, 0.05, 3))
        self.scaling = self.row("Scaling", QComboBox(), "How Mach and Re vary with CL")
        for key, text in SCALINGS.items():
            self.scaling.addItem(text, key)
        self.scaling.setCurrentIndex(list(SCALINGS).index(s.scaling))
        self.ncrit = self.row(
            "Ncrit", spin(s.ncrit, 0.1, 20.0, 0.5, 2), "eᴺ critical amplification (9 = average wind tunnel)"
        )
        self.xtr_top = self.row(
            "Trip top x/c", spin(s.xtr_top, 0.0, 1.0, 0.05, 3), "Forced transition; 1 = free"
        )
        self.xtr_bottom = self.row(
            "Trip bottom x/c", spin(s.xtr_bottom, 0.0, 1.0, 0.05, 3), "Forced transition; 1 = free"
        )
        self.iterations = self.row("Max iterations", int_spin(s.max_iterations, 5, 500, 5))

        reset = QPushButton("Re-initialise boundary layer")
        reset.setToolTip("Start the next point from a fresh boundary layer")
        reset.clicked.connect(session.reset_boundary_layer)
        self.form.addRow(reset)

        self.viscous.toggled.connect(self._apply)
        self.reynolds.value_changed.connect(self._apply)
        self.scaling.currentIndexChanged.connect(self._apply)
        for box in (self.mach, self.ncrit, self.xtr_top, self.xtr_bottom, self.iterations):
            box.valueChanged.connect(self._apply)

        self._viscous_only = (
            self.reynolds,
            self.ncrit,
            self.xtr_top,
            self.xtr_bottom,
            self.iterations,
            reset,
        )
        self._editors = (self.viscous, self.mach, self.scaling, *self._viscous_only)
        session.busy_changed.connect(self._set_busy)
        self._update_enabled()

    def _apply(self) -> None:
        settings = Settings(
            reynolds=self.reynolds.value() if self.viscous.isChecked() else 0.0,
            mach=self.mach.value(),
            ncrit=self.ncrit.value(),
            xtr_top=self.xtr_top.value(),
            xtr_bottom=self.xtr_bottom.value(),
            max_iterations=self.iterations.value(),
            scaling=self.scaling.currentData(),
        )
        self.session.set_settings(settings)
        self._update_enabled()

    def _update_enabled(self) -> None:
        busy = self.session.busy
        for w in self._viscous_only:
            w.setEnabled(self.viscous.isChecked() and not busy)

    def _set_busy(self, busy: bool) -> None:
        for w in self._editors:
            w.setEnabled(not busy)
        self._update_enabled()
