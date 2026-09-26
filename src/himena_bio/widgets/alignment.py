from __future__ import annotations

from typing import Any
from himena import WidgetDataModel
from qtpy import QtWidgets as QtW, QtGui, QtCore
from himena.plugins import validate_protocol
from himena.consts import MonospaceFontFamily
from himena_bio.consts import Type
from Bio.Align import Alignment, MultipleSeqAlignment, PairwiseAlignments


class QAlignmentView(QtW.QWidget):
    def __init__(self):
        super().__init__()
        self._ith = QAlignmentSpinBox()
        self._score = QtW.QLabel()
        self._view = QtW.QPlainTextEdit()
        self._view.setReadOnly(True)
        self._view.setWordWrapMode(QtGui.QTextOption.WrapMode.NoWrap)
        self._view.setFont(QtGui.QFont(MonospaceFontFamily, 9))
        layout = QtW.QVBoxLayout(self)
        layout.addWidget(self._ith)
        layout.addWidget(self._score)
        layout.addWidget(self._view)
        self._ith.valueChanged.connect(self._on_index_changed)
        self._alignments: Any = []
        self._model_type = Type.ALIGNMENT

    @validate_protocol
    def update_model(self, model: WidgetDataModel):
        value = model.value
        if isinstance(value, (Alignment, MultipleSeqAlignment)):
            value = [value]
        elif not isinstance(value, (PairwiseAlignments, list, tuple)):
            raise ValueError(f"Invalid alignment type: {type(value)}")
        self._alignments = value
        self._model_type = model.type
        self._ith.setMaximum(_num_alignments(value) - 1)
        self._ith.setValue(0)
        self._on_index_changed(0)

    @validate_protocol
    def to_model(self) -> WidgetDataModel:
        return WidgetDataModel(value=self._alignments, type=self.model_type())

    @validate_protocol
    def model_type(self) -> str:
        return self._model_type

    @validate_protocol
    def size_hint(self) -> tuple[int, int]:
        return 420, 500

    def _on_index_changed(self, index: int):
        try:
            aln = self._alignments[index]
        except (IndexError, StopIteration):
            # PairwiseAlignments may not know the number of alignments in advance
            self._ith.setMaximum(index - 1)
            return
        self._score.setText(_alignment_summary(aln))
        if isinstance(aln, MultipleSeqAlignment):
            self._view.setPlainText(format(aln, "clustal"))
        else:
            self._view.setPlainText(str(aln))


def _num_alignments(alignments) -> float:
    try:
        return len(alignments)
    except (OverflowError, TypeError):
        return float("inf")


def _alignment_summary(aln) -> str:
    texts = []
    if (score := getattr(aln, "score", None)) is not None:
        texts.append(f"Score = {score:.2f}")
    if isinstance(aln, Alignment) and len(aln.sequences) == 2:
        try:
            counts = aln.counts()
        except Exception:
            pass
        else:
            length = aln.length
            if length > 0:
                texts.append(
                    f"Identity = {counts.identities}/{length} "
                    f"({counts.identities / length:.1%})"
                )
                texts.append(
                    f"Gaps = {counts.gaps}/{length} ({counts.gaps / length:.1%})"
                )
    return ", ".join(texts)


class QAlignmentSpinBox(QtW.QWidget):
    valueChanged = QtCore.Signal(int)

    def __init__(self):
        super().__init__()
        layout = QtW.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self._left = QtW.QPushButton("◀")
        self._left.clicked.connect(self._on_prev)
        self._left.setFixedWidth(30)
        self._label = QtW.QLabel("0")
        self._label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self._right = QtW.QPushButton("▶")
        self._right.clicked.connect(self._on_next)
        self._right.setFixedWidth(30)

        layout.addWidget(self._left)
        layout.addWidget(self._label)
        layout.addWidget(self._right)

        self._value = 0
        self._max_value: float = float("inf")
        self._update_buttons()

    def _on_next(self):
        if self._value < self._max_value:
            self.setValue(self._value + 1)
            self.valueChanged.emit(self._value)

    def _on_prev(self):
        if self._value > 0:
            self.setValue(self._value - 1)
            self.valueChanged.emit(self._value)

    def value(self) -> int:
        return self._value

    def setValue(self, value: int):
        self._value = int(max(min(value, self._max_value), 0))
        self._update_buttons()

    def setMaximum(self, value: float):
        self._max_value = max(value, 0)
        if self._value > self._max_value:
            self.setValue(int(self._max_value))
        self._update_buttons()

    def _update_buttons(self):
        if self._max_value == float("inf"):
            self._label.setText(str(self._value))
        else:
            self._label.setText(f"{self._value} / {int(self._max_value)}")
        self._left.setEnabled(self._value > 0)
        self._right.setEnabled(self._value < self._max_value)
