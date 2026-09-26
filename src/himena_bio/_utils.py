from __future__ import annotations

from typing import TYPE_CHECKING
import re
from qtpy import QtGui
from cmap import Color

if TYPE_CHECKING:
    from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
    from Bio.SeqRecord import SeqRecord

_GRAY_PATTERN = re.compile(r"gray(\d+)")


def parse_ape_color(color: str) -> QtGui.QColor:
    if match := _GRAY_PATTERN.match(color):
        val = round(255 * int(match.group(1)) / 100)
        return QtGui.QColor(val, val, val)
    return QtGui.QColor(Color(color).hex)


def feature_color(feature: SeqFeature) -> QtGui.QColor | None:
    """Get the ApE color of the feature (reverse color for the minus strand)."""
    from himena_bio.consts import ApeAnnotation

    fw = feature.qualifiers.get(ApeAnnotation.FWCOLOR)
    rv = feature.qualifiers.get(ApeAnnotation.RVCOLOR)
    colors = (rv or fw) if feature.location.strand == -1 else (fw or rv)
    if not colors:
        return None
    try:
        return parse_ape_color(colors[0])
    except Exception:
        return None


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
        part = loc.parts[nth]
        start, end = int(part.start), int(part.end)
    else:
        raise NotImplementedError(f"Unknown location type: {type(loc)}")
    return start, end


def topology(rec: SeqRecord) -> str:
    return rec.annotations.get("topology", "linear")


def slice_seq_record(rec: SeqRecord, index: slice) -> SeqRecord:
    """Slice a SeqRecord, keeping the features partially covered by the slice.

    Unlike `SeqRecord.__getitem__`, features that overlap with the slice boundaries
    are clipped instead of being dropped.
    """
    start, stop, step = index.indices(len(rec))
    if step != 1:
        raise ValueError("Only slices with step 1 are supported.")
    return _slice_segments(rec, [(start, max(start, stop))])


def slice_circular(rec: SeqRecord, start: int, stop: int) -> SeqRecord:
    """Slice a circular SeqRecord from `start` to `stop`, going across the origin.

    If `start >= stop`, the returned record is `rec[start:] + rec[:stop]`. Features
    split by the origin are joined together.
    """
    if start < stop:
        return _slice_segments(rec, [(start, stop)])
    return _slice_segments(rec, [(start, len(rec)), (0, stop)])


def rotate_seq_record(rec: SeqRecord, origin: int) -> SeqRecord:
    """Rotate a circular SeqRecord so that the position `origin` becomes 0."""
    origin %= max(len(rec), 1)
    if origin == 0:
        return _slice_segments(rec, [(0, len(rec))])
    return slice_circular(rec, origin, origin)


def _slice_segments(rec: SeqRecord, segments: list[tuple[int, int]]) -> SeqRecord:
    """Concatenate the given segments of a record, mapping features properly."""
    from Bio.SeqRecord import SeqRecord

    seq = rec.seq[0:0]
    for start, stop in segments:
        seq += rec.seq[start:stop]
    out = SeqRecord(
        seq,
        id=rec.id,
        name=rec.name,
        description=rec.description,
        dbxrefs=list(rec.dbxrefs),
    )
    if "molecule_type" in rec.annotations:
        # we need this for GenBank/EMBL etc output
        out.annotations["molecule_type"] = rec.annotations["molecule_type"]
    for feat in rec.features:
        try:
            loc = map_location(feat.location, segments)
        except TypeError:
            # Will fail on UnknownPosition
            continue
        if loc is not None:
            out.features.append(copy_feature(feat, loc))
    for key, value in rec.letter_annotations.items():
        new_value = value[0:0]
        for start, stop in segments:
            new_value += value[start:stop]
        out.letter_annotations[key] = new_value
    return out


def map_location(
    loc: SimpleLocation | CompoundLocation,
    segments: list[tuple[int, int]],
) -> SimpleLocation | CompoundLocation | None:
    """Map a location on the source sequence to the concatenated segments.

    Each part of the location is clipped to each segment. Parts that became
    contiguous in the new coordinates (such as features split by the origin of a
    circular sequence) are merged. Returns None if nothing is left.
    """
    from Bio.SeqFeature import SimpleLocation, CompoundLocation

    offsets = []
    _offset = 0
    for start, stop in segments:
        offsets.append(_offset)
        _offset += stop - start

    new_parts: list[SimpleLocation] = []
    for part in loc.parts:
        p_start, p_end = int(part.start), int(part.end)
        pieces: list[tuple[int, int, int]] = []  # (source start, new start, new end)
        for (s_start, s_stop), offset in zip(segments, offsets):
            c_start, c_end = max(p_start, s_start), min(p_end, s_stop)
            if c_start < c_end:
                pieces.append(
                    (c_start, c_start - s_start + offset, c_end - s_start + offset)
                )
        # keep the biological order of the pieces
        pieces.sort(key=lambda x: x[0], reverse=part.strand == -1)
        for _, new_start, new_end in pieces:
            new_part = SimpleLocation(new_start, new_end, strand=part.strand)
            if new_parts and _is_contiguous(new_parts[-1], new_part):
                prev = new_parts.pop()
                new_part = SimpleLocation(
                    min(prev.start, new_part.start),
                    max(prev.end, new_part.end),
                    strand=part.strand,
                )
            new_parts.append(new_part)
    if len(new_parts) == 0:
        return None
    elif len(new_parts) == 1:
        return new_parts[0]
    return CompoundLocation(new_parts, operator=getattr(loc, "operator", "join"))


def _is_contiguous(prev: SimpleLocation, next: SimpleLocation) -> bool:
    if prev.strand != next.strand:
        return False
    if next.strand == -1:
        return next.end == prev.start
    return prev.end == next.start


def copy_feature(
    feature: SeqFeature, loc: SimpleLocation | CompoundLocation | None = None
) -> SeqFeature:
    """Shallow copy of a SeqFeature, optionally with a new location."""
    from Bio.SeqFeature import SeqFeature

    if loc is None:
        loc = feature.location
    return SeqFeature(
        location=loc,
        type=feature.type,
        id=feature.id,
        qualifiers={k: v for k, v in feature.qualifiers.items()},
    )
