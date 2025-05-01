from __future__ import annotations

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
        _layout.addWidget(_line)
        _layout.addWidget(_also_find_reverse_complement)
        _layout.addWidget(_btn_prev)
        _layout.addWidget(_btn_next)
        _layout.addWidget(_btn_hide)
        _btn_prev.clicked.connect(self._btn_prev_clicked)
        _btn_next.clicked.connect(self._btn_next_clicked)
        _btn_hide.clicked.connect(self.hide)
        _line.textChanged.connect(self._find_update)
        self._line_edit = _line
        self._also_find_reverse_complement = _also_find_reverse_complement
        self._btn_prev = _btn_prev
        self._btn_next = _btn_next
        self._btn_hide = _btn_hide
        self._parent_widget = parent

    def show(self):
        super().show()
        self._line_edit.setFocus()

    def _btn_prev_clicked(self):
        self._find_prev()
        self._line_edit.setFocus()

    def _btn_next_clicked(self):
        self._find_next()
        self._line_edit.setFocus()

    def _find_prev(self):
        text = self._line_edit.text()
        if text == "":
            return
        qtext = self._seqedit
        cursor = qtext.textCursor()
        pos = cursor.selectionEnd() - 1
        all_text = qtext.toPlainText()
        if self._rv_find(cursor, text, all_text[:pos], 0):
            qtext.setTextCursor(cursor)
        elif self._rv_find(cursor, text, all_text[pos:], pos):
            qtext.setTextCursor(cursor)

    def _rv_find(
        self,
        cursor: QtGui.QTextCursor,
        text: str,
        ref_text: str,
        offset: int,
    ) -> bool:
        ind_fw = ref_text.rfind(text)
        if self._also_find_reverse_complement.isChecked():
            ind_rc = ref_text.rfind(str(Seq(text).reverse_complement()))
            ind = max(ind_fw, ind_rc)
        else:
            ind = ind_fw
        if ind >= 0:
            cursor.setPosition(ind + offset)
            cursor.setPosition(
                ind + offset + len(text),
                QtGui.QTextCursor.MoveMode.KeepAnchor,
            )
            return True
        else:
            return False

    def _find_next(self):
        text = self._line_edit.text()
        if text == "":
            return
        qtext = self._seqedit
        cursor = qtext.textCursor()
        pos = cursor.selectionStart() + 1
        all_text = qtext.toPlainText()
        if self._fw_find(cursor, text, all_text[pos:], pos):
            qtext.setTextCursor(cursor)
        elif self._fw_find(cursor, text, all_text[:pos], 0):
            qtext.setTextCursor(cursor)

    def _fw_find(
        self,
        cursor: QtGui.QTextCursor,
        text: str,
        ref_text: str,
        offset: int,
    ) -> bool:
        ind_fw = ref_text.find(text)
        if self._also_find_reverse_complement.isChecked():
            ind_rc = ref_text.find(str(Seq(text).reverse_complement()))
            if ind_fw < 0 and ind_rc < 0:
                ind = -1
            elif ind_fw < 0:
                ind = ind_rc
            elif ind_rc < 0:
                ind = ind_fw
            else:
                ind = min(ind_fw, ind_rc)
        else:
            ind = ind_fw
        if ind >= 0:
            cursor.setPosition(ind + offset)
            cursor.setPosition(
                ind + offset + len(text),
                QtGui.QTextCursor.MoveMode.KeepAnchor,
            )
            return True
        else:
            return False

    _find_update = _find_next

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
