from __future__ import annotations

from typing import TYPE_CHECKING
from qtpy import QtWidgets as QtW
from qtpy.QtCore import Qt
from superqt import QElidingLabel

from himena.qt.magicgui._toggle_switch import QLabeledToggleSwitch

if TYPE_CHECKING:
    from himena_bio.widgets.editor import QMultiSeqEdit


class QSeqControl(QtW.QWidget):
    def __init__(self, edit: QMultiSeqEdit):
        super().__init__()
        self._edit = edit
        layout = QtW.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._feature_label = QElidingLabel()
        self._feature_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self._feature_label, stretch=10)
        self._is_one_start = QLabeledToggleSwitch()
        self._is_one_start.setText("1-start")
        self._is_one_start.setToolTip(
            "Toggle between 0-based (Python style) and 1-based (ApE style) numbering for base positions."
        )
        self._is_one_start.toggled.connect(edit._selection_changed)
        layout.addWidget(self._is_one_start)

        self._hover_pos = QVParam("Pos", width=36, tooltip="Hovered position")
        layout.addWidget(self._hover_pos)

        self._sel = QVParam(
            "Selection", width=70, tooltip="Current text cursor selection"
        )
        layout.addWidget(self._sel)

        self._length = QVParam(
            "Length", width=45, tooltip="Length of the current text cursor selection"
        )
        layout.addWidget(self._length)

        self._tm = QVParam("Tm", width=40, tooltip="The calculated melting temperature")
        layout.addWidget(self._tm)
        self._percent_gc = QVParam(
            "%GC",
            width=45,
            tooltip="Percentage of GC content at the text cursor selection.",
        )
        layout.addWidget(self._percent_gc)

        self._topology = QtW.QComboBox()
        self._topology.setToolTip("Topology of the sequence.")
        self._topology.addItems(["linear", "circular"])
        self._topology.setFixedWidth(72)
        layout.addWidget(self._topology)


class QVParam(QtW.QWidget):
    def __init__(self, label: str, width: int = 96, tooltip: str = ""):
        super().__init__()
        layout = QtW.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(1)
        self._label = QtW.QLabel(label)
        self._value = QtW.QLabel()
        font = self._label.font()
        font.setPointSize(8)
        font.setBold(True)
        self._label.setFont(font)
        layout.addWidget(self._label)
        layout.addWidget(self._value)
        self.setFixedWidth(width)
        self.setToolTip(tooltip)

    def set_value(self, value: str):
        self._value.setText(value)

    def set_visible(self, visible: bool):
        self._label.setVisible(visible)
        self._value.setVisible(visible)
