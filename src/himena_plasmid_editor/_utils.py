from __future__ import annotations

from typing import TYPE_CHECKING
import re
from qtpy import QtGui
from cmap import Color

if TYPE_CHECKING:
    from Bio.SeqFeature import SeqFeature

_GRAY_PATTERN = re.compile(r"gray(\d+)")


def parse_ape_color(color: str) -> QtGui.QColor:
    if match := _GRAY_PATTERN.match(color):
        val = round(255 * int(match.group(1)) / 100)
        return QtGui.QColor(val, val, val)
    return QtGui.QColor(Color(color).hex)


def get_feature_label(feature: SeqFeature) -> str:
    d = feature.qualifiers
    out = d.get("label", d.get("locus_tag", d.get("ApEinfo_label", None)))
    if isinstance(out, list):
        out = out[0]
    if not isinstance(out, str):
        out = feature.type
    return out
