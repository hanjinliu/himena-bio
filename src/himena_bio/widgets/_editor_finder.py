from __future__ import annotations

from bisect import bisect_left
from qtpy import QtWidgets as QtW, QtGui
from qtpy.QtCore import Qt
from typing import TYPE_CHECKING
from superqt import QToggleSwitch
from Bio.Seq import Seq

if TYPE_CHECKING:
    from himena_bio.widgets.editor import QSeqEdit


class QSeqFinder(QtW.QWidget):
    def __init__(self, parent, textedit: QSeqEdit):
        super().__init__(parent)
        self._seqedit = textedit
        _layout = QtW.QHBoxLayout(self)
        _layout.setContentsMargins(0, 0, 0, 0)
        _layout.setSpacing(2)
        _line = QtW.QLineEdit()
        _line.setPlaceholderText("Find (case insensitive)")
        _also_find_reverse_complement = QToggleSwitch("RevComp")
        _also_find_reverse_complement.setToolTip(
            "Also find the reverse complement of the sequence."
        )
        _btn_prev = QtW.QPushButton("◀")
        _btn_next = QtW.QPushButton("▶")
        _btn_hide = QtW.QPushButton("✕")
        _btn_prev.setFixedSize(18, 18)
        _btn_next.setFixedSize(18, 18)
        _btn_hide.setFixedSize(18, 18)
        _btn_prev.setToolTip("Find previous (Shift+Enter)")
        _btn_next.setToolTip("Find next (Enter)")
        _layout.addWidget(_line)
        _layout.addWidget(_also_find_reverse_complement)
        _layout.addWidget(_btn_prev)
        _layout.addWidget(_btn_next)
        _layout.addWidget(_btn_hide)
        _btn_prev.clicked.connect(self._btn_prev_clicked)
        _btn_next.clicked.connect(self._btn_next_clicked)
        _btn_hide.clicked.connect(self.hide)
        _line.textChanged.connect(self._find_update)
        _also_find_reverse_complement.toggled.connect(self._find_update)
        self._line_edit = _line
        self._also_find_reverse_complement = _also_find_reverse_complement
        self._btn_prev = _btn_prev
        self._btn_next = _btn_next
        self._btn_hide = _btn_hide

    def show(self):
        super().show()
        self._line_edit.setFocus()
        self._line_edit.selectAll()

    def set_text(self, text: str):
        """Set the query text (this will trigger a search from the cursor)."""
        self._line_edit.setText(text)

    def _btn_prev_clicked(self):
        self._find_prev()
        self._line_edit.setFocus()

    def _btn_next_clicked(self):
        self._find_next()
        self._line_edit.setFocus()

    def _all_matches(self) -> tuple[list[int], int]:
        """Return all the (sorted) start positions of matches and the query length."""
        query = self._line_edit.text().upper()
        if query == "":
            return [], 0
        ref = self._seqedit.toPlainText().upper()
        queries = {query}
        if self._also_find_reverse_complement.isChecked():
            try:
                queries.add(str(Seq(query).reverse_complement()))
            except ValueError:  # e.g. protein sequence
                pass
        positions: set[int] = set()
        for q in queries:
            pos = ref.find(q)
            while pos >= 0:
                positions.add(pos)
                pos = ref.find(q, pos + 1)
        return sorted(positions), len(query)

    def _find_next(self, include_current: bool = False) -> bool:
        """Select the next match after the current cursor.

        If `include_current` is true, the match starting at the current selection is
        also a candidate. This is needed for incremental search, where the query is
        extended while the current selection should be kept.
        """
        positions, length = self._all_matches()
        if length == 0:
            self._set_found(True)
            return False
        cursor = self._seqedit.textCursor()
        start = cursor.selectionStart()
        if cursor.hasSelection() and not include_current:
            start += 1
        idx = bisect_left(positions, start)
        if idx >= len(positions):
            idx = 0  # wrap around
        return self._select_match(positions, idx, length)

    def _find_prev(self) -> bool:
        """Select the previous match before the current cursor."""
        positions, length = self._all_matches()
        if length == 0:
            self._set_found(True)
            return False
        cursor = self._seqedit.textCursor()
        idx = bisect_left(positions, cursor.selectionStart()) - 1
        # idx == -1 is the last one (wrap around)
        return self._select_match(positions, idx, length)

    def _select_match(self, positions: list[int], idx: int, length: int) -> bool:
        if len(positions) == 0:
            self._set_found(False)
            return False
        pos = positions[idx]
        cursor = self._seqedit.textCursor()
        cursor.setPosition(pos)
        cursor.setPosition(pos + length, QtGui.QTextCursor.MoveMode.KeepAnchor)
        self._seqedit.setTextCursor(cursor)
        self._set_found(True)
        return True

    def _set_found(self, found: bool):
        if found:
            self._line_edit.setStyleSheet("")
        else:
            self._line_edit.setStyleSheet("QLineEdit { color: red; }")

    def _find_update(self):
        self._find_next(include_current=True)

    def keyPressEvent(self, a0: QtGui.QKeyEvent) -> None:
        if a0.key() == Qt.Key.Key_Escape:
            self.hide()
            self._seqedit.setFocus()
        elif a0.key() in (Qt.Key.Key_Enter, Qt.Key.Key_Return):
            if a0.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                self._find_prev()
            else:
                self._find_next()
        return super().keyPressEvent(a0)
