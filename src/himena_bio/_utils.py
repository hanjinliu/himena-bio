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


def slice_seq_record(rec: SeqRecord, index: slice):
    from Bio.SeqFeature import SimpleLocation

    parent_length = len(rec)
    answer = rec._from_validated(
        rec.seq[index],
        id=rec.id,
        name=rec.name,
        description=rec.description,
    )
    if "molecule_type" in rec.annotations:
        # This will still apply, and we need it for GenBank/EMBL etc output
        answer.annotations["molecule_type"] = rec.annotations["molecule_type"]

    start, stop, step = index.indices(parent_length)
    if step == 1:
        # Select relevant features, add them with shifted locations
        # assert str(self.seq)[index] == str(self.seq)[start:stop]
        for feat in rec.features:
            try:
                if feat.location.start < start < feat.location.end:
                    feat_clipped = type(feat)(
                        location=SimpleLocation(0, feat.location.end - start),
                        type=feat.type,
                        id=feat.id,
                        qualifiers=feat.qualifiers,
                    )
                    answer.features.append(feat_clipped)
                elif feat.location.start < stop < feat.location.end:
                    feat_clipped = type(feat)(
                        location=SimpleLocation(
                            max(feat.location.start - start, 0), stop - start
                        ),
                        type=feat.type,
                        id=feat.id,
                        qualifiers=feat.qualifiers,
                    )
                    answer.features.append(feat_clipped)
                elif start <= feat.location.start and feat.location.end <= stop:
                    answer.features.append(feat._shift(-start))
            except TypeError:
                # Will fail on UnknownPosition
                pass

    # Slice all the values to match the sliced sequence
    # (this should also work with strides, even negative strides):
    for key, value in rec.letter_annotations.items():
        answer.letter_annotations[key] = value[index]

    return answer
