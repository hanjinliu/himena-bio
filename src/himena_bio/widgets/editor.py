from __future__ import annotations

from typing import Callable, Iterable
from qtpy import QtWidgets as QtW
from qtpy import QtCore, QtGui
from qtpy.QtCore import Qt
from superqt import QElidingLabel
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as _Tm

from magicgui.widgets import Dialog
from cmap import Color
from himena import WidgetDataModel
from himena.widgets import set_status_tip
from himena.qt.magicgui import get_type_map
from himena.qt.magicgui._toggle_switch import QLabeledToggleSwitch
from himena.consts import MonospaceFontFamily
from himena.plugins import validate_protocol
from himena_bio.consts import Keys, ApeAnnotation, SeqMeta, Type
from himena_bio._utils import parse_ape_color, get_feature_label
from himena_bio.widgets._feature_view import QFeatureView


def _char_to_qt_key(char: str) -> Qt.Key:
    if char.isalnum():
        return getattr(Qt.Key, f"Key_{char.upper()}")
    if char == " ":
        return Qt.Key.Key_Space
    if char == "*":
        return Qt.Key.Key_Asterisk
    if char == "-":
        return Qt.Key.Key_Minus
    raise NotImplementedError(f"Unsupported character: {char}")


_KEYS_MOVE = frozenset(
    [Qt.Key.Key_Left, Qt.Key.Key_Right, Qt.Key.Key_Up, Qt.Key.Key_Down,
    Qt.Key.Key_Home, Qt.Key.Key_End, Qt.Key.Key_PageUp, Qt.Key.Key_PageDown]
)  # fmt: skip


class QSeqEdit(QtW.QPlainTextEdit):
    hovered = QtCore.Signal(object)  # int or None
    edited = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFont(QtGui.QFont(MonospaceFontFamily, 9))
        self.setWordWrapMode(QtGui.QTextOption.WrapMode.WrapAnywhere)
        self.setLineWrapMode(QtW.QPlainTextEdit.LineWrapMode.WidgetWidth)
        self.setUndoRedoEnabled(True)
        self.setTabChangesFocus(True)
        self.setAcceptDrops(True)
        self.setMouseTracking(True)
        self._record = SeqRecord(Seq(""))
        self.set_keys_allowed(Keys.DNA)
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._show_context_menu)

    def set_record(self, record: SeqRecord):
        self._record = record
        self.setPlainText(str(record.seq))
        self.update_highlight()
        return None

    def update_highlight(self):
        cursor = self.textCursor()
        cursor.movePosition(QtGui.QTextCursor.MoveOperation.Start)
        cursor.movePosition(
            QtGui.QTextCursor.MoveOperation.End, QtGui.QTextCursor.MoveMode.KeepAnchor
        )
        cursor.setCharFormat(QtGui.QTextCharFormat())
        cursor = self.textCursor()
        format = QtGui.QTextCharFormat()
        for feature in self._record.features:
            if isinstance(loc := feature.location, SimpleLocation):
                parts = [(int(loc.start), int(loc.end))]
            elif isinstance(loc := feature.location, CompoundLocation):
                parts = [(int(sub.start), int(sub.end)) for sub in loc]
            else:
                parts = []
            for start, end in parts:
                cursor.setPosition(start)
                cursor.setPosition(end, QtGui.QTextCursor.MoveMode.KeepAnchor)
                if colors := feature.qualifiers.get(ApeAnnotation.FWCOLOR):
                    color = parse_ape_color(colors[0])
                    format.setBackground(color)
                    r, g, b, _ = color.getRgbF()
                    luminance = 0.299 * r + 0.587 * g + 0.114 * b
                    format.setForeground(
                        Qt.GlobalColor.white
                        if luminance < 0.5
                        else Qt.GlobalColor.black
                    )
                    cursor.setCharFormat(format)
        self.edited.emit()
        return None

    def to_record(self) -> SeqRecord:
        return SeqRecord(
            seq=Seq(self.toPlainText()),
            id=self._record.id,
            name=self._record.name,
            description=self._record.description,
            dbxrefs=self._record.dbxrefs,
            features=self._record.features,
            annotations=self._record.annotations,
            letter_annotations=self._record.letter_annotations,
        )

    def keyPressEvent(self, event: QtGui.QKeyEvent):
        _mod = event.modifiers()
        _key = event.key()
        if _key in _KEYS_MOVE:
            return super().keyPressEvent(event)
        if _mod & Qt.KeyboardModifier.ControlModifier:
            _shift_on = bool(_mod & Qt.KeyboardModifier.ShiftModifier)
            if _key == Qt.Key.Key_C:
                self._copy_selection(_shift_on)
            elif _key == Qt.Key.Key_X:
                self._cut_selection(_shift_on)
            elif _key == Qt.Key.Key_V:
                self._paste(_shift_on)
            return
        elif _mod & Qt.KeyboardModifier.AltModifier:
            return None
        else:
            if _key in self._keys_allowed:
                _key_char = event.text()
                self.insert_text(_key_char, self.textCursor())
            elif _key == Qt.Key.Key_Backspace:
                self._backspace_event()
            elif _key == Qt.Key.Key_Delete:
                self._delete_event()
            return None

    def mouseMoveEvent(self, e):
        self.hovered.emit(self.cursorForPosition(e.pos()).position())
        return super().mouseMoveEvent(e)

    def leaveEvent(self, a0):
        self.hovered.emit(None)
        return super().leaveEvent(a0)

    def _copy_selection(self, rev_comp: bool = False):
        if rev_comp:
            seq = str(Seq(self.textCursor().selectedText()).reverse_complement())
        else:
            seq = self.textCursor().selectedText()
        clipboard = QtW.QApplication.clipboard()
        clipboard.setText(seq)

    def _paste(self, rev_comp: bool = False):
        clipboard = QtW.QApplication.clipboard()
        seq = clipboard.text()
        if rev_comp:
            seq = str(Seq(seq).reverse_complement())
        self.textCursor().insertText(seq)

    def _cut_selection(self, rev_comp: bool = False):
        self._copy_selection(rev_comp)
        self.textCursor().removeSelectedText()

    def insert_text(self, text: str, cursor: QtGui.QTextCursor):
        start = cursor.selectionStart()
        if cursor.hasSelection():
            self.delete_text(cursor)
            cursor.clearSelection()
        nchars = len(text)
        for feature in self._record.features:
            _shift_feature(feature, start, nchars)

        cursor.insertText(text)
        self.update_highlight()

    def delete_text(self, cursor: QtGui.QTextCursor):
        if not cursor.hasSelection():
            raise ValueError("No text selected")
        start = cursor.selectionStart()
        nchars = cursor.selectionEnd() - start
        for feature in self._record.features:
            _shift_feature(feature, start, -nchars)
        cursor.removeSelectedText()
        self.update_highlight()

    def _backspace_event(self):
        cursor = self.textCursor()
        if not cursor.hasSelection():
            cursor.movePosition(
                QtGui.QTextCursor.MoveOperation.PreviousCharacter,
                QtGui.QTextCursor.MoveMode.KeepAnchor,
            )
        self.delete_text(cursor)

    def _delete_event(self):
        cursor = self.textCursor()
        if not cursor.hasSelection():
            cursor.movePosition(
                QtGui.QTextCursor.MoveOperation.NextCharacter,
                QtGui.QTextCursor.MoveMode.KeepAnchor,
            )
        self.delete_text(cursor)

    def _make_context_menu(self) -> QtW.QMenu:
        menu = QtW.QMenu()
        cursor = self.textCursor()
        pos = cursor.position()

        cut_action = menu.addAction("Cut", self._cut_selection)
        cut_rc_action = menu.addAction(
            "Cut (RevComp)", lambda: self._cut_selection(True)
        )
        copy_action = menu.addAction("Copy", self._copy_selection)
        copy_rc_action = menu.addAction(
            "Copy (RevComp)", lambda: self._copy_selection(True)
        )
        paste_action = menu.addAction("Paste", self._paste)
        menu.addSeparator()
        new_feature_action = menu.addAction("New Feature", self._new_feature)
        edit_feature_action = menu.addAction(
            "Edit Feature", lambda: self._edit_feature(self._get_front_feature())
        )
        del_feature_action = menu.addAction(
            "Delete Feature", lambda: self._delete_feature(self._get_front_feature())
        )
        move_feature_front_action = menu.addAction(
            "Move Feature Front",
            lambda: self._move_feature_front(self._get_front_feature()),
        )
        move_feature_back_action = menu.addAction(
            "Move Feature Back",
            lambda: self._move_feature_back(self._get_front_feature()),
        )
        if not cursor.hasSelection():
            cut_action.setEnabled(False)
            cut_rc_action.setEnabled(False)
            copy_action.setEnabled(False)
            copy_rc_action.setEnabled(False)
            new_feature_action.setEnabled(False)
            move_feature_front_action.setEnabled(False)
            move_feature_back_action.setEnabled(False)

        if len(self._features_under_pos(pos)) == 0:
            edit_feature_action.setEnabled(False)
            del_feature_action.setEnabled(False)

        if QtW.QApplication.clipboard().text() == "":
            paste_action.setEnabled(False)
        return menu

    def _show_context_menu(self, pos: QtCore.QPoint):
        cursor_pos = self.cursorForPosition(pos).position()
        cursor = self.textCursor()
        start, end = cursor.selectionStart(), cursor.selectionEnd()
        if not start <= cursor_pos < end:
            cursor.setPosition(cursor_pos)
            self.setTextCursor(cursor)
        menu = self._make_context_menu()
        menu.exec(self.mapToGlobal(pos))

    def _features_under_pos(self, pos: int) -> list[SeqFeature]:
        return [feat for feat in self._record.features if pos in feat]

    def _new_feature(self):
        cursor = self.textCursor()
        start, end = cursor.selectionStart(), cursor.selectionEnd()
        feature = SeqFeature(
            location=SimpleLocation(start, end),
            **self._feature_qualifiers_from_dialog(),
        )
        self._record.features.append(feature)
        self.update_highlight()

    def _get_front_feature(self) -> SeqFeature:
        pos = self.textCursor().position()
        features = self._features_under_pos(pos)
        return features[-1]

    def _edit_feature(self, feature):
        kwargs = self._feature_qualifiers_from_dialog(
            name=feature.type,
            fcolor=feature.qualifiers.get(ApeAnnotation.FWCOLOR, ["cyan"])[0],
            rcolor=feature.qualifiers.get(ApeAnnotation.RVCOLOR, ["cyan"])[0],
        )
        feature.type = feature.id = kwargs["type"]
        feature.qualifiers.update(kwargs["qualifiers"])
        self.update_highlight()

    def _delete_feature(self, feature):
        self._record.features.remove(feature)
        self.update_highlight()

    def _move_feature_front(self, feature):
        self._record.features.remove(feature)
        self._record.features.append(feature)
        self.update_highlight()

    def _move_feature_back(self, feature):
        self._record.features.remove(feature)
        self._record.features.insert(0, feature)
        self.update_highlight()

    def _feature_qualifiers_from_dialog(
        self,
        name: str = "Unnamed",
        fcolor: str = "cyan",
        rcolor: str = "cyan",
    ) -> dict:
        typemap = get_type_map()
        w_name = typemap.create_widget(value=name, label="Name")
        w_fcolor = typemap.create_widget(value=Color(fcolor), label="Forward Color")
        w_rcolor = typemap.create_widget(value=Color(rcolor), label="Reverse Color")
        dlg = Dialog(widgets=[w_name, w_fcolor, w_rcolor])
        dlg.exec()
        return {
            "type": w_name.value,
            "qualifiers": {
                ApeAnnotation.LABEL: [w_name.value],
                ApeAnnotation.FWCOLOR: [Color(w_fcolor.value).hex],
                ApeAnnotation.RVCOLOR: [Color(w_rcolor.value).hex],
            },
        }

    def set_keys_allowed(self, keys: Iterable[str]):
        _keys = frozenset(_char_to_qt_key(char) for char in keys)
        self._keys_allowed = _keys


class QMultiSeqEdit(QtW.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QtW.QVBoxLayout(self)

        self._seq_choices = QtW.QComboBox(self)
        layout.addWidget(self._seq_choices)

        self._feature_view = QFeatureView(self)
        self._feature_view.setFixedHeight(40)
        self._feature_view.clicked.connect(self._on_view_clicked)
        self._feature_view.hovered.connect(self._on_view_hovered)
        layout.addWidget(self._feature_view)

        self._seq_edit = QSeqEdit(self)
        self._seq_edit.hovered.connect(self._seq_edit_hovered)
        self._seq_edit.edited.connect(self._seq_edited)
        layout.addWidget(self._seq_edit)

        self._comment = QtW.QPlainTextEdit(self)
        self._comment.setFont(QtGui.QFont(MonospaceFontFamily, 9))
        self._comment.setWordWrapMode(QtGui.QTextOption.WrapMode.WordWrap)
        self._comment.setFixedHeight(80)
        layout.addWidget(self._comment)

        self._tm_method: Callable[[Seq], float] = _Tm.Tm_GC

        self._control = QSeqControl(self)
        self._seq_edit.selectionChanged.connect(self._selection_changed)
        self._model_type = Type.DNA

    @validate_protocol
    def update_model(self, model: WidgetDataModel):
        recs = model.value
        if len(recs) == 0:
            return
        choices: list[str] = []
        recs_normed = []
        for rec in recs:
            if not isinstance(rec, SeqRecord):
                rec = SeqRecord(Seq(rec))
            choices.append(rec.name)
            recs_normed.append(rec)
        self._seq_choices.addItems(choices)
        self._set_record(recs_normed[0])
        self._seq_choices.setVisible(len(recs_normed) > 1)
        self._model_type = model.type
        if model.type == Type.DNA:
            self._seq_edit.set_keys_allowed(Keys.DNA)
        elif model.type == Type.RNA:
            self._seq_edit.set_keys_allowed(Keys.RNA)
        elif model.type == Type.PROTEIN:
            self._seq_edit.set_keys_allowed(Keys.PROTEIN)
        else:
            raise NotImplementedError(f"Unsupported model type: {model.type}")

    @validate_protocol
    def to_model(self) -> WidgetDataModel:
        cursor = self._seq_edit.textCursor()
        return WidgetDataModel(
            value=[self._seq_edit.to_record()],
            type=self._model_type,
            metadata=SeqMeta(
                current_index=self._seq_choices.currentIndex(),
                selection=(cursor.selectionStart(), cursor.selectionEnd()),
            ),
        )

    @validate_protocol
    def model_type(self) -> str:
        return self._model_type

    @validate_protocol
    def control_widget(self) -> QSeqControl:
        return self._control

    @validate_protocol
    def size_hint(self) -> tuple[int, int]:
        return 400, 400

    @validate_protocol
    def widget_added_callback(self):
        self._feature_view.auto_range()

    def setFocus(self):
        self._seq_edit.setFocus()

    def _set_record(self, record: SeqRecord):
        self._seq_edit.set_record(record)
        self._feature_view.set_record(record)
        cursor = self._seq_edit.textCursor()
        cursor.setPosition(0)
        self._seq_edit.setTextCursor(cursor)
        comment = record.annotations.get(ApeAnnotation.COMMENT, "")
        self._comment.setPlainText(_remove_ape_meta(comment))

    def _selection_changed(self):
        cursor = self._seq_edit.textCursor()
        selection = cursor.selectedText()
        offset = 1 if self._control._is_one_start.isChecked() else 0
        if selection:
            self._control._sel.set_value(
                f"{cursor.selectionStart() + offset} - {cursor.selectionEnd() + offset}"
            )
            self._control._length.set_value(str(len(selection)))
            gc_count = sum(1 for nucleotide in selection if nucleotide in "GC")
            self._control._percent_gc.set_value(
                f"{(gc_count / len(selection)) * 100:.2f}%"
            )
            tm = self._tm_method(Seq(selection))
            if tm < 0:
                self._control._tm.set_value("-- °C")
            else:
                self._control._tm.set_value(f"{tm:.1f}°C")
            _visible = True
        else:
            pos = cursor.position() + offset
            self._control._sel.set_value(f"{pos - 1} | {pos}")
            _visible = False
        self._control._length.set_visible(_visible)
        self._control._tm.set_visible(_visible)
        self._control._percent_gc.set_visible(_visible)

        # translate the selection and show in status tip
        adjusted_length = len(selection) // 3 * 3
        amino_acids = str(Seq(selection[:adjusted_length]).translate())
        if len(amino_acids) > 20:
            amino_acids = amino_acids[:10] + "..." + amino_acids[-10:]
        set_status_tip(amino_acids, duration=5)

    def _seq_edit_hovered(self, pos: int | None):
        if pos is None:
            self._control._feature_label.setText("")
            self._control._hover_pos.set_value("")
            return
        offset = 1 if self._control._is_one_start.isChecked() else 0
        self._control._hover_pos.set_value(str(pos + offset))
        feature_labels = [
            get_feature_label(feat) for feat in self._seq_edit._features_under_pos(pos)
        ]
        tip = ", ".join(feature_labels)
        self._control._feature_label.setText(tip)

    def _seq_edited(self):
        self._feature_view.set_record(self._seq_edit.to_record())

    def _on_view_clicked(self, feature: SeqFeature | None, nth: int):
        if feature is None:
            cursor = self._seq_edit.textCursor()
            cursor.setPosition(nth)
            self._seq_edit.ensureCursorVisible()
            self._seq_edit.setTextCursor(cursor)
            return
        if isinstance(loc := feature.location, SimpleLocation):
            start, end = int(loc.start), int(loc.end)
        elif isinstance(loc := feature.location, CompoundLocation):
            start, end = int(loc[nth].start), int(loc[nth].end)
        else:
            return
        cursor = self._seq_edit.textCursor()
        cursor.setPosition(start)
        self._seq_edit.ensureCursorVisible()
        cursor.setPosition(end, QtGui.QTextCursor.MoveMode.KeepAnchor)
        self._seq_edit.setTextCursor(cursor)

    def _on_view_hovered(self, feature: SeqFeature | None, nth: int):
        if isinstance(feature, SeqFeature):
            label = get_feature_label(feature)
            self._control._feature_label.setText(label)
        else:
            self._control._feature_label.setText("")
        self._control._hover_pos.set_value(str(nth))


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


def _shift_feature(feature: SeqFeature, start: int, shift: int):
    if isinstance(loc := feature.location, SimpleLocation):
        if loc.end > start:
            start_new, end_new = loc.start, loc.end
            end_new += shift
            if loc.start >= start:
                start_new += shift
            loc_new = SimpleLocation(start_new, end_new)
            feature.location = loc_new
    elif isinstance(loc := feature.location, CompoundLocation):
        new_parts: list[SimpleLocation] = []
        for part in loc:
            part_new: SimpleLocation = part
            if part.end > start:
                start_new, end_new = part.start, part.end
                end_new += shift
                if part.start >= start:
                    start_new += shift
                part_new = SimpleLocation(start_new, end_new)
            new_parts.append(part_new)
        loc_new = CompoundLocation(new_parts)
        feature.location = loc_new


def _remove_ape_meta(comment: str) -> str:
    lines = comment.splitlines()
    if lines and lines[-1].startswith("ApEinfo:methylated:"):
        lines.pop(-1)
    return "\n".join(lines)
