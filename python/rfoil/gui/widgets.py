"""Small reusable widgets."""

from __future__ import annotations

from typing import TypeVar

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import (
    QAbstractSpinBox,
    QDoubleSpinBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

SpinBox = TypeVar("SpinBox", bound=QAbstractSpinBox)


class Card(QFrame):
    """A titled, rounded panel holding a form."""

    def __init__(self, title: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Card")
        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 10, 12, 12)
        outer.setSpacing(8)
        heading = QLabel(title.upper())
        heading.setObjectName("CardTitle")
        outer.addWidget(heading)
        self.form = QFormLayout()
        self.form.setHorizontalSpacing(10)
        self.form.setVerticalSpacing(6)
        self.form.setLabelAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        outer.addLayout(self.form)
        self.body = outer

    def row(self, label: str, widget: QWidget, tip: str = "") -> QWidget:
        if tip:
            widget.setToolTip(tip)
        self.form.addRow(muted(label), widget)
        return widget


def muted(text: str) -> QLabel:
    label = QLabel(text)
    label.setObjectName("Muted")
    return label


def caption(text: str) -> QLabel:
    """Small muted label used under a row of inputs."""
    label = muted(text)
    label.setStyleSheet("font-size: 8pt;")
    return label


def hbox(*widgets: QWidget, spacing: int = 6) -> QWidget:
    """Lays widgets out side by side in a margin-less container."""
    holder = QWidget()
    layout = QHBoxLayout(holder)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(spacing)
    for w in widgets:
        layout.addWidget(w)
    return holder


class NumberEdit(QLineEdit):
    """Line edit for a float that accepts scientific notation (``1e6``, ``3.5e5``)."""

    value_changed = Signal(float)

    def __init__(self, value: float, minimum: float = 0.0, maximum: float = 1e12) -> None:
        super().__init__()
        validator = QDoubleValidator(minimum, maximum, 6, self)
        validator.setNotation(QDoubleValidator.Notation.ScientificNotation)
        self.setValidator(validator)
        self.setAlignment(Qt.AlignmentFlag.AlignRight)
        self._value = value
        self.set_value(value)
        self.editingFinished.connect(self._commit)

    def value(self) -> float:
        return self._value

    def set_value(self, value: float) -> None:
        self._value = value
        self.setText(f"{value:g}")

    def _commit(self) -> None:
        try:
            value = float(self.text())
        except ValueError:
            self.set_value(self._value)
            return
        if value != self._value:
            self._value = value
            self.value_changed.emit(value)
        self.set_value(value)


def spin(
    value: float, lo: float, hi: float, step: float, decimals: int = 2, suffix: str = ""
) -> QDoubleSpinBox:
    box = QDoubleSpinBox()
    box.setRange(lo, hi)
    box.setDecimals(decimals)
    box.setSingleStep(step)
    box.setValue(value)
    if suffix:
        box.setSuffix(suffix)
    return _plain(box)


def int_spin(value: int, lo: int, hi: int, step: int = 1) -> QSpinBox:
    box = QSpinBox()
    box.setRange(lo, hi)
    box.setSingleStep(step)
    box.setValue(value)
    return _plain(box)


def _plain(box: SpinBox) -> SpinBox:
    """Right-aligned, commits on Enter or focus loss; arrows hidden by the theme.

    (``NoButtons`` would be simpler, but with a styled border Qt then
    collapses the text field to zero width.)
    """
    box.setAlignment(Qt.AlignmentFlag.AlignRight)
    box.setKeyboardTracking(False)
    box.lineEdit().setTextMargins(4, 0, 4, 0)
    return box


class ValueLabel(QLabel):
    """Monospaced read-out of a number."""

    def __init__(self, fmt: str = "{:.4f}") -> None:
        super().__init__("-")
        self.setObjectName("Value")
        self.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self._fmt = fmt

    def show_value(self, value: float | None) -> None:
        self.setText("-" if value is None else self._fmt.format(value))
