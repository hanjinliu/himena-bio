from __future__ import annotations

from contextlib import contextmanager
from typing import Callable, Iterable
from qtpy import QtWidgets as QtW
from qtpy import QtCore, QtGui
from qtpy.QtCore import Qt
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from Bio.SeqUtils import MeltingTemp

from magicgui.widgets import Dialog
from cmap import Color
from himena import WidgetDataModel
from himena.widgets import set_status_tip
from himena.types import is_subtype
from himena.qt.magicgui import get_type_map
from himena.consts import MonospaceFontFamily
from himena.plugins import validate_protocol
from himena.utils.collections import UndoRedoStack

from himena_bio.consts import Keys, ApeAnnotation, SeqMeta, Type
from himena_bio._utils import parse_ape_color, get_feature_label
from himena_bio.widgets._feature_view import QFeatureView
from himena_bio.widgets._editor_finder import QSeqFinder
from himena_bio.widgets._base import char_to_qt_key, infer_seq_type
from himena_bio.widgets._editor_control import QSeqControl
from himena_bio.widgets import _editor_actions as _ea

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
        self._actions: dict[str, QtW.QAction] = {}
        self._setup_actions()
        self._undo_redo_stack = UndoRedoStack[_ea.EditorAction](size=100)
        self._data_modified = False
        self._pend_update_highlight = False

        @self.edited.connect
        def _():
            self._data_modified = True

    def _setup_actions(self):
        """Setup actions for keyboard shortcuts and context menu."""

        def add_action(name: str, shortcut: str, callback):
            action = QtW.QAction(name, self)
            action.setShortcut(shortcut)
            action.triggered.connect(callback)
            action.setShortcutContext(Qt.ShortcutContext.WidgetShortcut)
            self._actions[name] = action
            self.addAction(action)

        add_action("Copy", "Ctrl+C", lambda: self._copy_selection(False))
        add_action("Cut", "Ctrl+X", lambda: self._cut_selection(False))
        add_action("Paste", "Ctrl+V", lambda: self._paste(False))
        add_action("Undo", "Ctrl+Z", self.undo)
        add_action("Redo", "Ctrl+Y", self.redo)
        add_action("Find", "Ctrl+F", self._open_finder)
        add_action("Copy (RevComp)", "Ctrl+Alt+C", lambda: self._copy_selection(True))
        add_action("Cut (RevComp)", "Ctrl+Alt+X", lambda: self._cut_selection(True))
        add_action("Paste (RevComp)", "Ctrl+Alt+V", lambda: self._paste(True))
        add_action("New Feature", "Alt+N", self._new_feature)
        add_action(
            "Edit Feature",
            "Alt+E",
            lambda: self._edit_feature(self._get_front_feature()),
        )
        add_action(
            "Delete Feature",
            "Alt+D",
            lambda: self._delete_feature(self._get_front_feature()),
        )
        add_action(
            "Move Feature Front",
            "Alt+F",
            lambda: self._move_feature_front(self._get_front_feature()),
        )
        add_action(
            "Move Feature Back",
            "Alt+B",
            lambda: self._move_feature_back(self._get_front_feature()),
        )

    def set_record(self, record: SeqRecord):
        self._record = record
        self.setPlainText(str(record.seq))
        self.update_highlight()
        return None

    def update_highlight(self):
        if self._pend_update_highlight:
            return
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
            _alt_on = bool(_mod & Qt.KeyboardModifier.AltModifier)
            if _key == Qt.Key.Key_C:
                self._copy_selection(_alt_on)
            elif _key == Qt.Key.Key_X:
                self._cut_selection(_alt_on)
            elif _key == Qt.Key.Key_A:
                self.selectAll()
            elif _key == Qt.Key.Key_V:
                self._paste(_alt_on)
            elif _key == Qt.Key.Key_Z and not _alt_on:
                self.undo()
            elif _key == Qt.Key.Key_Y and not _alt_on:
                self.redo()
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
        self.insert_text(seq, self.textCursor())

    def _cut_selection(self, rev_comp: bool = False):
        self._copy_selection(rev_comp)
        self.delete_text(self.textCursor())

    @contextmanager
    def pend_update_highlight(self):
        """Pend calling update_highlight in this context."""
        was_pending = self._pend_update_highlight
        self._pend_update_highlight = True
        try:
            yield
        finally:
            self._pend_update_highlight = was_pending

    def insert_text(
        self,
        text: str,
        cursor: QtGui.QTextCursor,
        record_undo: bool = True,
    ):
        """Insert the text using the given text cursor."""
        start = cursor.selectionStart()
        actions = []
        if cursor.hasSelection():
            actions.append(
                _ea.DeleteSeqAction(
                    start=cursor.selectionStart(),
                    length=cursor.selectionEnd() - cursor.selectionStart(),
                    seq=cursor.selectedText(),
                )
            )
            self.delete_text(cursor, record_undo=False)
            cursor.clearSelection()
        nchars = len(text)
        for index, feature in enumerate(self._record.features):
            if (feature_shifted := _shift_feature(feature, start, nchars)) is not None:
                action = _ea.EditFeatureAction(
                    index=index, old=feature, new=feature_shifted
                )
                action.apply(self)
                actions.append(action)
        actions.append(_ea.InsertSeqAction(pos=start, seq=text))
        cursor.insertText(text)
        self.update_highlight()
        if record_undo:
            self._undo_redo_stack.push(_ea.CompositeAction(actions))

    def delete_text(self, cursor: QtGui.QTextCursor, record_undo: bool = True):
        """Delete the selected text using the given text cursor."""
        if not cursor.hasSelection():
            raise ValueError("No text selected")
        start = cursor.selectionStart()
        nchars = cursor.selectionEnd() - start
        actions = [
            _ea.DeleteSeqAction(
                start=start,
                length=nchars,
                seq=cursor.selectedText(),
            )
        ]
        for index, feature in enumerate(self._record.features):
            if (feature_shifted := _shift_feature(feature, start, -nchars)) is not None:
                action = _ea.EditFeatureAction(
                    index=index, old=feature, new=feature_shifted
                )
                action.apply(self)
                actions.append(action)
        cursor.removeSelectedText()
        self.update_highlight()
        if record_undo:
            self._undo_redo_stack.push(_ea.CompositeAction(actions))

    def undo(self):
        """Undo the last edit action."""
        if action := self._undo_redo_stack.undo():
            action.invert().apply(self)
            with self.pend_update_highlight():
                self.update_highlight()
        return None

    def redo(self):
        """Redo the last undone edit action."""
        if action := self._undo_redo_stack.redo():
            action.apply(self)
            with self.pend_update_highlight():
                self.update_highlight()
        return None

    def _open_finder(self):
        self.parentWidget()._finder_widget.show()

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
        menu.addAction(self._actions["Cut"])
        menu.addAction(self._actions["Cut (RevComp)"])
        menu.addAction(self._actions["Copy"])
        menu.addAction(self._actions["Copy (RevComp)"])
        menu.addAction(self._actions["Paste"])
        menu.addSeparator()
        menu.addAction(self._actions["New Feature"])
        menu.addAction(self._actions["Edit Feature"])
        menu.addAction(self._actions["Delete Feature"])
        menu.addAction(self._actions["Move Feature Front"])
        menu.addAction(self._actions["Move Feature Back"])

        if QtW.QApplication.clipboard().text() == "":
            self._actions["Paste"].setEnabled(False)
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
        index = len(self._record.features)
        self._record.features.append(feature)
        self._undo_redo_stack.push(
            _ea.EditFeatureAction(index=index, old=None, new=feature)
        )
        self.update_highlight()

    def _get_front_feature(self) -> SeqFeature:
        pos = self.textCursor().position()
        features = self._features_under_pos(pos)
        return features[-1]

    def _edit_feature(self, feature: SeqFeature):
        kwargs = self._feature_qualifiers_from_dialog(
            name=feature.type,
            fcolor=feature.qualifiers.get(ApeAnnotation.FWCOLOR, ["cyan"])[0],
            rcolor=feature.qualifiers.get(ApeAnnotation.RVCOLOR, ["cyan"])[0],
        )
        index = self._record.features.index(feature)
        feature_new = copy_feature(feature)
        feature_new.type = feature_new.id = kwargs["type"]
        feature_new.qualifiers = feature_new.qualifiers.copy()
        feature_new.qualifiers.update(kwargs["qualifiers"])
        self._record.features[index] = feature_new
        self.update_highlight()
        self._undo_redo_stack.push(
            _ea.EditFeatureAction(index=index, old=feature, new=feature_new)
        )

    def _delete_feature(self, feature: SeqFeature):
        """Delete the feature from the current record."""
        index = self._record.features.index(feature)
        self._record.features.pop(index)
        self.update_highlight()
        self._undo_redo_stack.push(
            _ea.EditFeatureAction(index=index, old=feature, new=None)
        )

    def _move_feature_front(self, feature: SeqFeature):
        """Make sure the feature is visible by moving it the last of the list."""
        old = self._record.features.index(feature)
        new = len(self._record.features)
        if old == new:
            return
        action = _ea.MoveFeatureAction(old=old, new=new)
        action.apply(self)
        self.update_highlight()
        self._undo_redo_stack.push(action)

    def _move_feature_back(self, feature: SeqFeature):
        old = self._record.features.index(feature)
        new = 0
        if old == new:
            return
        action = _ea.MoveFeatureAction(old=old, new=new)
        action.apply(self)
        self.update_highlight()
        self._undo_redo_stack.push(action)

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
        _keys = frozenset(char_to_qt_key(char) for char in keys)
        self._keys_allowed = _keys


class QMultiSeqEdit(QtW.QWidget):
    """Sequence editor widget."""

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QtW.QVBoxLayout(self)

        self._seq_choices = QtW.QComboBox(self)
        self._feature_view = QFeatureView(self)
        self._seq_edit = QSeqEdit(self)
        self._finder_widget = QSeqFinder(self, self._seq_edit)

        layout.addWidget(self._seq_choices)

        self._feature_view.setFixedHeight(40)
        self._feature_view.clicked.connect(self._on_view_clicked)
        self._feature_view.hovered.connect(self._on_view_hovered)
        layout.addWidget(self._feature_view)

        layout.addWidget(self._finder_widget)
        self._finder_widget.hide()

        self._seq_edit.hovered.connect(self._seq_edit_hovered)
        self._seq_edit.edited.connect(self._seq_edited)
        layout.addWidget(self._seq_edit)

        self._comment = QtW.QPlainTextEdit(self)
        self._comment.setFont(QtGui.QFont(MonospaceFontFamily, 9))
        self._comment.setWordWrapMode(QtGui.QTextOption.WrapMode.WordWrap)
        self._comment.setFixedHeight(80)
        layout.addWidget(self._comment)

        self._tm_method: Callable[[Seq], float] = MeltingTemp.Tm_GC

        self._control = QSeqControl(self)
        self._seq_edit.selectionChanged.connect(self._selection_changed)
        self._seq_edit.cursorPositionChanged.connect(self._selection_changed)
        self._model_type = Type.DNA
        self._extension_default = ".ape"

    @validate_protocol
    def update_model(self, model: WidgetDataModel):
        recs = model.value
        if len(recs) == 0:
            return
        was_empty = self._seq_choices.count() == 0
        choices: list[str] = []
        recs_normed: list[SeqRecord] = []
        for rec in recs:
            if not isinstance(rec, SeqRecord):
                rec = SeqRecord(Seq(rec))
            choices.append(rec.name)
            recs_normed.append(rec)
        self._seq_choices.addItems(choices)
        self._set_record(recs_normed[0])
        self._seq_choices.setVisible(len(recs_normed) > 1)
        if model.type == Type.SEQS:
            _mtype = infer_seq_type(str(recs_normed[0].seq))
        else:
            _mtype = model.type
        if is_subtype(_mtype, Type.DNA):
            self._seq_edit.set_keys_allowed(Keys.DNA)
        elif is_subtype(_mtype, Type.RNA):
            self._seq_edit.set_keys_allowed(Keys.RNA)
        elif is_subtype(_mtype, Type.PROTEIN):
            self._seq_edit.set_keys_allowed(Keys.PROTEIN)
        else:
            raise NotImplementedError(f"Unsupported model type: {_mtype}")
        self._model_type = _mtype
        if ext := model.extension_default:
            self._extension_default = ext
        if was_empty:
            self._seq_edit._data_modified = False

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
            extension_default=self._extension_default,
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
    def is_modified(self) -> bool:
        return self._seq_edit._data_modified

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
            # nucleotide specific info
            if self.model_type() in (Type.DNA, Type.RNA):
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

        if self.model_type() in (Type.DNA, Type.RNA):
            # translate the selection and show in status tip
            adjusted_length = len(selection) // 3 * 3
            amino_acids = str(Seq(selection[:adjusted_length]).translate())
            if len(amino_acids) > 20:
                amino_acids = amino_acids[:10] + "..." + amino_acids[-10:]
            set_status_tip(amino_acids, duration=5)

        has_selection = cursor.hasSelection()
        self._seq_edit._actions["Cut"].setEnabled(has_selection)
        self._seq_edit._actions["Cut (RevComp)"].setEnabled(has_selection)
        self._seq_edit._actions["Copy"].setEnabled(has_selection)
        self._seq_edit._actions["Copy (RevComp)"].setEnabled(has_selection)
        self._seq_edit._actions["New Feature"].setEnabled(has_selection)

        has_feature = len(self._seq_edit._features_under_pos(cursor.position())) > 0
        self._seq_edit._actions["Edit Feature"].setEnabled(has_feature)
        self._seq_edit._actions["Delete Feature"].setEnabled(has_feature)

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


def _shift_feature(feature: SeqFeature, start: int, shift: int) -> SeqFeature | None:
    """Shift the feature by the given amount, return None if not shifted."""
    if isinstance(loc := feature.location, SimpleLocation):
        if loc.end > start:
            start_new, end_new = loc.start, loc.end
            end_new += shift
            if loc.start >= start:
                start_new += shift
            loc_new = SimpleLocation(start_new, end_new)
            return copy_feature(feature, loc_new)
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
        return copy_feature(feature, loc_new)
    return None


def _remove_ape_meta(comment: str) -> str:
    lines = comment.splitlines()
    if lines and lines[-1].startswith("ApEinfo:methylated:"):
        lines.pop(-1)
    return "\n".join(lines)


def copy_feature(
    feature: SeqFeature, loc: SimpleLocation | CompoundLocation | None = None
) -> SeqFeature:
    """Shallow copy of a SeqFeature, optionally with a new location."""
    if loc is None:
        loc = feature.location
    return SeqFeature(
        location=loc,
        type=feature.type,
        id=feature.id,
        qualifiers={k: v for k, v in feature.qualifiers.items()},
    )
