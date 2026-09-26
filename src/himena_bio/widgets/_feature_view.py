from __future__ import annotations

from typing import TYPE_CHECKING
from qtpy import QtWidgets as QtW
from qtpy import QtCore, QtGui
from qtpy.QtCore import Qt
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

from himena.widgets import set_clipboard
from himena.qt import qimage_to_ndarray
from himena_bio._utils import feature_to_slice, feature_color, get_feature_label
from himena_bio.widgets._base import QBaseGraphicsView

if TYPE_CHECKING:
    from himena_bio.widgets.editor import QMultiSeqEdit


class QFeatureRectitem(QtW.QGraphicsRectItem):
    def __init__(self, loc: SimpleLocation, feature: SeqFeature, nth: int = 0):
        super().__init__(float(loc.start), -0.5, float(loc.end - loc.start), 1)
        self._feature = feature
        self._nth = nth
        self.setAcceptedMouseButtons(Qt.MouseButton.LeftButton)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        pen = QtGui.QPen(QtGui.QColor(Qt.GlobalColor.gray), 1)
        pen.setCosmetic(True)
        self.setPen(pen)


class QFeatureItem(QtW.QGraphicsItemGroup):
    def __init__(self, feature: SeqFeature):
        super().__init__()
        self._feature = feature
        self._rects: list[QFeatureRectitem] = []
        self.setAcceptedMouseButtons(Qt.MouseButton.LeftButton)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        if (color := feature_color(feature)) is None:
            color = QtGui.QColor(Qt.GlobalColor.gray)
        label = get_feature_label(feature)
        for ith, part in enumerate(feature.location.parts):
            rect_item = QFeatureRectitem(part, feature, ith)
            rect_item.setBrush(QtGui.QBrush(color))
            rect_item.setToolTip(label)
            self._rects.append(rect_item)
            self.addToGroup(rect_item)


class QFeatureView(QBaseGraphicsView):
    """The interactive viewer for the features.

    This viewer renders the features of a sequence in a human-readable format like:
    ---[    ]-[ ]---
    """

    clicked = QtCore.Signal(object, int)  # feature or None, nth part or position
    hovered = QtCore.Signal(object, object)  # feature or None, position or None

    def __init__(self, parent: QMultiSeqEdit):
        super().__init__()
        self._mseq_edit = parent
        self.setMouseTracking(True)
        self.setStyleSheet("QFeatureView { border: none; }")
        self._center_line = QtW.QGraphicsLineItem(0, 0, 1, 0)
        pen = QtGui.QPen(QtGui.QColor(Qt.GlobalColor.gray), 2)
        pen.setCosmetic(True)
        self._record = SeqRecord(Seq(""))
        self._center_line.setPen(pen)
        self._feature_items: list[QFeatureItem] = []
        self.scene().addItem(self._center_line)

        self._drag_start: QtCore.QPoint | None = None
        self._drag_prev = QtCore.QPoint()
        self._last_btn = Qt.MouseButton.NoButton
        self._is_auto_range = True

    def set_record(self, record: SeqRecord):
        for item in self._feature_items:
            self.scene().removeItem(item)
        self._feature_items.clear()
        for feature in record.features:
            if feature.location is None:
                continue
            item = QFeatureItem(feature)
            self._feature_items.append(item)
            self.scene().addItem(item)
        self._center_line.setLine(0, 0, len(record.seq), 0)
        self._record = record
        if self._is_auto_range:
            self.auto_range()

    def wheelEvent(self, event: QtGui.QWheelEvent):
        if event.angleDelta().y() < 0:
            self.scale(0.9, 1)
        else:
            self.scale(1.1, 1)
        self._is_auto_range = False

    def auto_range(self):
        _len = max(self._center_line.line().x2(), 1)
        self.fitInView(QtCore.QRectF(0, -1, _len, 2))
        self._is_auto_range = True

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self._is_auto_range:
            self.auto_range()

    def _seq_pos(self, pos: QtCore.QPoint) -> int:
        """Convert the widget position to the sequence position."""
        x = self.mapToScene(pos).x()
        return int(min(max(round(x), 0), len(self._record.seq)))

    def leaveEvent(self, a0):
        self.hovered.emit(None, None)

    def mousePressEvent(self, event):
        self._drag_start = self._drag_prev = event.pos()
        self._last_btn = event.button()
        return super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QtGui.QMouseEvent):
        if self._drag_start is None:
            # is hovering
            if isinstance(item := self.itemAt(event.pos()), QFeatureRectitem):
                self.hovered.emit(item._feature, self._seq_pos(event.pos()))
            else:
                self.hovered.emit(None, self._seq_pos(event.pos()))
        else:
            pos = event.pos()
            dpos = pos - self._drag_prev
            self._drag_prev = pos
            if dpos.x() != 0:
                self._is_auto_range = False
            self.horizontalScrollBar().setValue(
                self.horizontalScrollBar().value() - dpos.x()
            )
        return super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self._drag_start is None:
            return super().mouseReleaseEvent(event)
        is_click = (event.pos() - self._drag_start).manhattanLength() < 5
        self._drag_start = None
        if is_click:
            item = self.itemAt(event.pos())
            if isinstance(item, QFeatureRectitem):
                self.clicked.emit(item._feature, item._nth)
            else:
                self.clicked.emit(None, self._seq_pos(event.pos()))
            if self._last_btn == Qt.MouseButton.RightButton:
                if isinstance(item, QFeatureRectitem):
                    menu = self._make_menu_for_feature(item._feature, item._nth)
                else:
                    menu = self._make_menu_for_blank()
                menu.exec(event.globalPos())
        self._last_btn = Qt.MouseButton.NoButton
        return super().mouseReleaseEvent(event)

    def _make_menu_for_feature(self, feature: SeqFeature, nth: int) -> QtW.QMenu:
        seq_edit = self._mseq_edit._seq_edit
        menu = QtW.QMenu()
        menu.addAction("Copy", lambda: self._copy_feature(feature, nth))
        menu.addAction("Edit", lambda: seq_edit._edit_feature(feature))
        menu.addAction("Delete", lambda: seq_edit._delete_feature(feature))
        menu.addAction("Move Front", lambda: seq_edit._move_feature_front(feature))
        menu.addAction("Move Back", lambda: seq_edit._move_feature_back(feature))
        return menu

    def _make_menu_for_blank(self) -> QtW.QMenu:
        menu = QtW.QMenu()
        menu.addAction("Reset View", self.auto_range)
        menu.addAction("Copy as image", self._copy_as_image)
        return menu

    def _copy_feature(self, feature: SeqFeature, nth: int):
        x0, x1 = feature_to_slice(feature, nth)
        set_clipboard(text=str(self._record.seq[x0:x1]), internal_data=feature)
        return

    def _copy_as_image(self):
        arr = qimage_to_ndarray(self.grab().toImage()).copy()
        set_clipboard(image=arr)
