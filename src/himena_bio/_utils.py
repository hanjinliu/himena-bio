from __future__ import annotations

from typing import TYPE_CHECKING
import re
from qtpy import QtGui
from cmap import Color

if TYPE_CHECKING:
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqRecord import SeqRecord

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


def feature_to_slice(feature: SeqFeature, nth: int) -> tuple[int, int]:
    from Bio.SeqFeature import SimpleLocation, CompoundLocation

    if isinstance(loc := feature.location, SimpleLocation):
        start, end = int(loc.start), int(loc.end)
    elif isinstance(loc := feature.location, CompoundLocation):
        start, end = int(loc[nth].start), int(loc[nth].end)
    else:
        raise NotImplementedError(f"Unknown location type: {type(loc)}")
    return start, end


def topology(rec: SeqRecord) -> str:
    return rec.annotations.get("topology", "linear")
