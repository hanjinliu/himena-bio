from __future__ import annotations

from logging import getLogger
from qtpy import QtWidgets as QtW
from qtpy import QtCore, QtGui
from qtpy.QtCore import Qt
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from himena_plasmid_editor.consts import ApeAnnotation
from himena_plasmid_editor._utils import parse_ape_color, get_feature_label

_LOGGER = getLogger(__name__)


class QBaseGraphicsScene(QtW.QGraphicsScene):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._grab_source: QtW.QGraphicsItem | None = None

    def grabSource(self) -> QtW.QGraphicsItem | None:
        return self._grab_source

    def setGrabSource(self, item: QtW.QGraphicsItem | None):
        self._grab_source = item


class QBaseGraphicsView(QtW.QGraphicsView):
    def __init__(self):
        scene = QBaseGraphicsScene()
        super().__init__(scene)
        self.setAlignment(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignHCenter)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setMouseTracking(True)

    def addItem(self, item: QtW.QGraphicsItem):
        self.scene().addItem(item)

    def scene(self) -> QBaseGraphicsScene:
        return super().scene()


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
    clicked = QtCore.Signal(object, int)
    hovered = QtCore.Signal(object, int)

    def __init__(self):
        super().__init__()
        self.setMouseTracking(True)
        self.setStyleSheet("QFeatureView { border: none; }")
        self._center_line = QtW.QGraphicsLineItem(0, 0, 1, 0)
        pen = QtGui.QPen(QtGui.QColor(Qt.GlobalColor.gray), 2)
        pen.setCosmetic(True)
        self._center_line.setPen(pen)
        self._feature_items: list[QFeatureItem] = []
        self.scene().addItem(self._center_line)

        self._drag_start = QtCore.QPoint()
        self._drag_prev = QtCore.QPoint()

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
        self.fitInView(QtCore.QRectF(0, -1, _len, 2))

    def wheelEvent(self, event: QtGui.QWheelEvent):
        if event.angleDelta().y() < 0:
            self.scale(0.9, 1)
        else:
            self.scale(1.1, 1)

    def leaveEvent(self, a0):
        self.hovered.emit(None, 0)

    def mousePressEvent(self, event):
        self._drag_start = self._drag_prev = event.pos()
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
            if isinstance(item := self.itemAt(event.pos()), QFeatureRectitem):
                self.clicked.emit(item._feature, item._nth)
            else:
                pos = event.pos().x()
                self.clicked.emit(None, pos)
        return super().mouseReleaseEvent(event)
