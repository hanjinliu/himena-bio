from __future__ import annotations

from typing import TYPE_CHECKING
from logging import getLogger
from qtpy import QtWidgets as QtW
from qtpy import QtCore, QtGui
from qtpy.QtCore import Qt
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation

from himena.widgets import set_clipboard
from himena.qt import qimage_to_ndarray
from himena_bio.consts import ApeAnnotation
from himena_bio._utils import (
    feature_to_slice,
    parse_ape_color,
    get_feature_label,
)
from himena_bio.widgets._base import QBaseGraphicsView

if TYPE_CHECKING:
    from himena_bio.widgets.editor import QMultiSeqEdit

_LOGGER = getLogger(__name__)


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
        if colors := feature.qualifiers.get(ApeAnnotation.FWCOLOR):
            color = parse_ape_color(colors[0])
        else:
            color = QtGui.QColor(Qt.GlobalColor.gray)
        for loc in feature.location.parts:
            if isinstance(loc, SimpleLocation):
                rect_item = QFeatureRectitem(loc, feature)
                rect_item.setBrush(QtGui.QBrush(color))
                rect_item.setToolTip(get_feature_label(feature))
                self._rects.append(rect_item)
                self.addToGroup(rect_item)
            elif isinstance(loc, CompoundLocation):
                for ith, part in enumerate(loc.parts):
                    rect_item = QFeatureRectitem(part, feature, ith)
                    rect_item.setBrush(QtGui.QBrush(color))
                    rect_item.setToolTip(get_feature_label(feature))
                    self._rects.append(rect_item)
                    self.addToGroup(rect_item)


class QFeatureView(QBaseGraphicsView):
    """The interactive viewer for the features.

    This viewer renders the features of a sequence in a human-readable format like:
    ---[    ]-[ ]---
    """

    clicked = QtCore.Signal(object, int)
    hovered = QtCore.Signal(object, int)

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

        self._drag_start = QtCore.QPoint()
        self._drag_prev = QtCore.QPoint()
        self._last_btn = Qt.MouseButton.NoButton

    def set_record(self, record: SeqRecord):
        for item in self._feature_items:
            self.scene().removeItem(item)
        self._feature_items.clear()
        for feature in record.features:
            item = QFeatureItem(feature)
            self._feature_items.append(item)
            self.scene().addItem(item)
        _len = len(record.seq) - 1
        self._center_line.setLine(0, 0, _len, 0)
        self._record = record
        self.auto_range()

    def wheelEvent(self, event: QtGui.QWheelEvent):
        if event.angleDelta().y() < 0:
            self.scale(0.9, 1)
        else:
            self.scale(1.1, 1)

    def auto_range(self):
        _len = self._center_line.line().x2()
        self.fitInView(QtCore.QRectF(0, -1, _len, 2))

    def leaveEvent(self, a0):
        self.hovered.emit(None, 0)

    def mousePressEvent(self, event):
        self._drag_start = self._drag_prev = event.pos()
        self._last_btn = event.button()
        return super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QtGui.QMouseEvent):
        if self._drag_start.isNull():
            # is hovering
            if isinstance(item := self.itemAt(event.pos()), QFeatureRectitem):
                self.hovered.emit(item._feature, event.pos().x())
            else:
                self.hovered.emit(None, event.pos().x())
        else:
            pos = event.pos()
            dpos = pos - self._drag_prev
            self._drag_prev = pos
            self.horizontalScrollBar().setValue(
                self.horizontalScrollBar().value() - dpos.x()
            )
        return super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        self._drag_start = QtCore.QPoint()
        ds = self._drag_start - self._drag_prev
        is_click = ds.x() + ds.y() < 5
        if is_click:
            item = self.itemAt(event.pos())
            if isinstance(item, QFeatureRectitem):
                self.clicked.emit(item._feature, item._nth)
            else:
                pos = event.pos().x()
                self.clicked.emit(None, pos)
            if self._last_btn == Qt.MouseButton.RightButton:
                if isinstance(item, QFeatureRectitem):
                    menu = self._make_menu_for_feature(item._feature, item._nth)
                else:
                    menu = self._make_menu_for_blank()
                menu.exec(event.globalPos())
        self._last_btn = Qt.MouseButton.NoButton
        return super().mouseReleaseEvent(event)

    def _make_menu_for_feature(self, feature: SeqFeature, nth: int) -> QtW.QMenu:
        menu = QtW.QMenu()
        menu.addAction("Copy", lambda: self._copy_feature(feature, nth))
        menu.addAction("Edit", lambda: self._mseq_edit._seq_edit._edit_feature(feature))
        menu.addAction(
            "Delete", lambda: self._mseq_edit._seq_edit._delete_feature(feature)
        )
        menu.addAction(
            "Move Front", lambda: self._mseq_edit._seq_edit._move_feature_front(feature)
        )
        menu.addAction(
            "Move Back", lambda: self._mseq_edit._seq_edit._move_feature_back(feature)
        )
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
