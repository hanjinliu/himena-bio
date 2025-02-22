from __future__ import annotations

import numpy as np
from qtpy import QtWidgets as QtW
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from himena import WidgetDataModel, StandardType
from himena.standards import plotting as hplt
from himena.plugins import validate_protocol
from himena.consts import MonospaceFontFamily
from himena_builtins.qt.plot import model_matplotlib_canvas
from himena_bio.consts import Type


class QAB1View(QtW.QWidget):
    def __init__(self):
        super().__init__()
        self._canvas = model_matplotlib_canvas()
        self._model_type = Type.DNA
        self._extension_default = ".ab1"
        self._record = SeqRecord(Seq(""))
        layout = QtW.QVBoxLayout(self)
        layout.addWidget(self._canvas)

    @validate_protocol
    def update_model(self, model: WidgetDataModel):
        self._model_type = model.type
        if ext := model.extension_default:
            self._extension_default = ext
        rec = model.value
        if isinstance(rec, list) and len(rec) == 1 and isinstance(rec[0], SeqRecord):
            rec = rec[0]
        else:
            assert isinstance(rec, SeqRecord)
        self._record = rec

        fig = hplt.figure()
        xloc = np.array(rec.annotations["abif_raw"]["PLOC1"])
        qscore = np.array(rec.letter_annotations["phred_quality"])
        sig_a = np.array(rec.annotations["abif_raw"]["DATA9"])
        sig_t = np.array(rec.annotations["abif_raw"]["DATA10"])
        sig_c = np.array(rec.annotations["abif_raw"]["DATA11"])
        sig_g = np.array(rec.annotations["abif_raw"]["DATA12"])
        sig_max = np.concatenate([sig_a, sig_t, sig_c, sig_g]).max()
        qscore_normed = qscore * sig_max * 0.88 / qscore.max()
        zeros = np.zeros_like(xloc)
        fig.band(xloc, zeros, qscore_normed, color="#80808054", width=0, style="-")
        fig.plot(np.arange(sig_a.size), sig_a, width=1, color="red")
        fig.plot(np.arange(sig_t.size), sig_t, width=1, color="green")
        fig.plot(np.arange(sig_c.size), sig_c, width=1, color="blue")
        fig.plot(np.arange(sig_g.size), sig_g, width=1, color="orange")
        fig.text(
            xloc,
            np.full(xloc.size, sig_max * 1.1),
            rec.seq,
            size=10,
            family=MonospaceFontFamily,
        )
        if xloc.size <= 50:
            fig.axes.x.lim = (0, xloc[-1])
        else:
            fig.axes.x.lim = (0, xloc[50])
        fig.axes.y.lim = (0, sig_max * 1.2)
        self._canvas.update_model(WidgetDataModel(value=fig, type=StandardType.PLOT))

    @validate_protocol
    def to_model(self) -> WidgetDataModel:
        return WidgetDataModel(value=[self._record], type=self.model_type())

    @validate_protocol
    def model_type(self) -> str:
        return self._model_type

    @validate_protocol
    def size_hint(self) -> tuple[int, int]:
        return 600, 220

    @validate_protocol
    def control_widget(self) -> QtW.QWidget:
        return self._canvas.control_widget()
