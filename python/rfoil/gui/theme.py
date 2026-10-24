"""Colours and styling shared by widgets and plots."""

from __future__ import annotations

from dataclasses import dataclass

import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QFont, QGuiApplication, QPalette
from PySide6.QtWidgets import QApplication


@dataclass(frozen=True)
class Theme:
    name: str
    background: str
    surface: str
    raised: str
    border: str
    text: str
    muted: str
    accent: str
    upper: str
    lower: str
    inviscid: str
    good: str
    bad: str
    series: tuple[str, ...]


DARK = Theme(
    name="dark",
    background="#15171c",
    surface="#1c1f26",
    raised="#252932",
    border="#323743",
    text="#e6e9ef",
    muted="#959eae",
    accent="#5aa9ff",
    upper="#5aa9ff",
    lower="#ffb454",
    inviscid="#7d8697",
    good="#4cc38a",
    bad="#ff6b6b",
    series=("#5aa9ff", "#ffb454", "#4cc38a", "#ff6b9a", "#b48cff", "#3cd0d0", "#e0c341", "#ff7f50"),
)

LIGHT = Theme(
    name="light",
    background="#f4f5f8",
    surface="#ffffff",
    raised="#eceef3",
    border="#d6dae2",
    text="#1b1f27",
    muted="#5b6474",
    accent="#2f6fde",
    upper="#2f6fde",
    lower="#d9822b",
    inviscid="#8a93a3",
    good="#1f9d61",
    bad="#d64545",
    series=("#2f6fde", "#d9822b", "#1f9d61", "#d6457a", "#7b52d9", "#1b9fa3", "#a88a12", "#e0603a"),
)

_current: Theme = DARK


def current() -> Theme:
    return _current


def system_prefers_dark() -> bool:
    hints = QGuiApplication.styleHints()
    return hints.colorScheme() == Qt.ColorScheme.Dark


def apply(app: QApplication, theme: Theme) -> None:
    """Applies ``theme`` to the application palette, style sheet and plots."""
    global _current
    _current = theme
    app.setStyle("Fusion")
    app.setPalette(_palette(theme))
    app.setStyleSheet(_style_sheet(theme))
    font = QFont(app.font())
    font.setPointSizeF(max(font.pointSizeF(), 9.5))
    app.setFont(font)
    pg.setConfigOptions(antialias=True, background=theme.surface, foreground=theme.muted)


def pen(color: str, width: float = 1.6, dashed: bool = False) -> pg.QtGui.QPen:
    p = pg.mkPen(color=color, width=width)
    if dashed:
        p.setStyle(Qt.PenStyle.DashLine)
    return p


def style_plot(plot: pg.PlotItem, x_label: str, y_label: str) -> None:
    """Uniform axis styling for a plot."""
    t = current()
    plot.showGrid(x=True, y=True, alpha=0.12)
    for side, label in (("bottom", x_label), ("left", y_label)):
        axis = plot.getAxis(side)
        axis.enableAutoSIPrefix(False)
        axis.setLabel(label, color=t.muted)
        axis.setPen(pg.mkPen(t.border))
        axis.setTextPen(pg.mkPen(t.muted))
    plot.getViewBox().setBackgroundColor(t.surface)


def _palette(t: Theme) -> QPalette:
    p = QPalette()
    roles = {
        QPalette.ColorRole.Window: t.background,
        QPalette.ColorRole.WindowText: t.text,
        QPalette.ColorRole.Base: t.surface,
        QPalette.ColorRole.AlternateBase: t.raised,
        QPalette.ColorRole.Text: t.text,
        QPalette.ColorRole.Button: t.raised,
        QPalette.ColorRole.ButtonText: t.text,
        QPalette.ColorRole.Highlight: t.accent,
        QPalette.ColorRole.HighlightedText: "#ffffff",
        QPalette.ColorRole.ToolTipBase: t.raised,
        QPalette.ColorRole.ToolTipText: t.text,
        QPalette.ColorRole.PlaceholderText: t.muted,
        QPalette.ColorRole.Mid: t.border,
        QPalette.ColorRole.Link: t.accent,
    }
    for role, color in roles.items():
        p.setColor(role, QColor(color))
    p.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, QColor(t.muted))
    p.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, QColor(t.muted))
    return p


def _style_sheet(t: Theme) -> str:
    return f"""
    QMainWindow, QDialog {{ background: {t.background}; }}
    QToolTip {{ background: {t.raised}; color: {t.text}; border: 1px solid {t.border}; padding: 4px; }}

    QFrame#Card {{
        background: {t.surface};
        border: 1px solid {t.border};
        border-radius: 8px;
    }}
    QLabel#CardTitle {{
        color: {t.muted};
        font-weight: 600;
        font-size: 8.5pt;
    }}
    QLabel#Muted {{ color: {t.muted}; }}
    QLabel#Value {{ font-family: "JetBrains Mono", "Cascadia Mono", "DejaVu Sans Mono", monospace; }}
    QLabel#Headline {{ font-size: 13pt; font-weight: 600; }}

    QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
        background: {t.background};
        border: 1px solid {t.border};
        border-radius: 5px;
        padding: 3px 6px;
        min-height: 20px;
        selection-background-color: {t.accent};
    }}
    QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{ border-color: {t.accent}; }}
    /* spin boxes are typed into; hide their arrow buttons */
    QSpinBox::up-button, QSpinBox::down-button,
    QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {{ width: 0; border: none; }}

    QPushButton {{
        background: {t.raised};
        border: 1px solid {t.border};
        border-radius: 5px;
        padding: 5px 12px;
    }}
    QPushButton:hover {{ border-color: {t.accent}; }}
    QPushButton:disabled {{ color: {t.muted}; }}
    QPushButton#Primary {{
        background: {t.accent};
        border-color: {t.accent};
        color: #ffffff;
        font-weight: 600;
    }}
    QPushButton#Primary:disabled {{ background: {t.raised}; border-color: {t.border}; color: {t.muted}; }}

    QTabWidget::pane {{ border: 1px solid {t.border}; border-radius: 8px; background: {t.surface}; top: -1px; }}
    QTabBar::tab {{
        background: transparent;
        color: {t.muted};
        padding: 7px 16px;
        border: none;
        border-bottom: 2px solid transparent;
    }}
    QTabBar::tab:selected {{ color: {t.text}; border-bottom-color: {t.accent}; }}
    QTabBar::tab:hover {{ color: {t.text}; }}

    QTableView {{
        background: {t.surface};
        gridline-color: {t.border};
        border: 1px solid {t.border};
        border-radius: 6px;
    }}
    QHeaderView::section {{
        background: {t.raised};
        color: {t.muted};
        border: none;
        border-right: 1px solid {t.border};
        padding: 4px 6px;
    }}
    QPlainTextEdit {{
        background: {t.surface};
        border: 1px solid {t.border};
        border-radius: 6px;
        font-family: "JetBrains Mono", "Cascadia Mono", "DejaVu Sans Mono", monospace;
        font-size: 9pt;
    }}
    QScrollArea {{ border: none; background: transparent; }}
    QScrollArea > QWidget > QWidget {{ background: transparent; }}
    QSplitter::handle {{ background: {t.background}; }}
    QScrollBar:vertical {{ background: transparent; width: 10px; margin: 0; }}
    QScrollBar:horizontal {{ background: transparent; height: 10px; margin: 0; }}
    QScrollBar::handle {{ background: {t.border}; border-radius: 4px; min-height: 24px; min-width: 24px; }}
    QScrollBar::handle:hover {{ background: {t.muted}; }}
    QScrollBar::add-line, QScrollBar::sub-line {{ width: 0; height: 0; }}
    QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; }}
    QCheckBox::indicator {{ width: 14px; height: 14px; }}
    QStatusBar {{ color: {t.muted}; }}
    QProgressBar {{
        border: 1px solid {t.border};
        border-radius: 4px;
        background: {t.background};
        max-height: 8px;
    }}
    QProgressBar::chunk {{ background: {t.accent}; border-radius: 3px; }}
    """
