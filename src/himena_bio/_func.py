from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from Bio.Seq import Seq
from himena_bio._utils import (
    topology,
    slice_seq_record,
    slice_circular,
    rotate_seq_record,
    copy_feature,
)


def _find_all(ref: str, pattern: str) -> list[int]:
    """Find all the (possibly overlapping) occurrences of the pattern."""
    out = []
    pos = ref.find(pattern)
    while pos >= 0:
        out.append(pos)
        pos = ref.find(pattern, pos + 1)
    return out


def find_match(
    vec: Seq,
    seq: Seq,
    min_match: int = 15,
    circular: bool = False,
) -> list[SimpleLocation]:
    r"""Find all the primer binding sites.

    _____ full match  ... OK
    ____/ flanking region contained ... OK
    __/\_ mismatch ... NG

    Parameters
    ----------
    vec : Seq
        The template sequence.
    seq : str or Seq
        Sequence of primer.
    min_match : int, optional
        The minimun length of match, by default 15.
    circular : bool, default False
        If true, the template is considered circular and primers binding across the
        origin are also found. In this case, the `end` of the returned locations may
        exceed the length of the template.

    Returns
    -------
    list of SimpleLocation
        Binding sites on the template. `strand` is 1 for the forward binding and -1
        for the reverse binding.
    """
    if min_match <= 0:
        raise ValueError("`min_match` must be positive value")

    vec_str = str(vec).upper()
    seq_str = str(seq).upper()
    seq_rc = str(Seq(seq_str).reverse_complement())
    min_match = min(min_match, len(seq_str))
    nvec = len(vec_str)
    if nvec == 0 or min_match == 0:
        return []
    if circular:
        # the middle copy is the canonical one, others are used to extend matches
        ref = vec_str * 3
        lo, hi = nvec, 2 * nvec
    else:
        ref = vec_str
        lo, hi = 0, nvec

    matches: list[SimpleLocation] = []

    # 1: forward check (the 3' end of the primer is the right end)
    for pos in _find_all(ref, seq_str[-min_match:]):
        if not lo <= pos < hi:
            continue
        start = pos
        prpos = len(seq_str) - min_match  # position on seq
        while start > 0 and prpos > 0 and ref[start - 1] == seq_str[prpos - 1]:
            start -= 1
            prpos -= 1
        matches.append(_make_location(start - lo, pos + min_match - lo, 1, nvec))

    # 2: reverse check (the 3' end of the primer is the left end)
    for pos in _find_all(ref, seq_rc[:min_match]):
        if not lo <= pos < hi:
            continue
        end = pos + min_match
        prpos = min_match  # position on seq_rc
        while end < len(ref) and prpos < len(seq_rc) and ref[end] == seq_rc[prpos]:
            end += 1
            prpos += 1
        matches.append(_make_location(pos - lo, end - lo, -1, nvec))

    return matches


def _make_location(start: int, end: int, strand: int, size: int) -> SimpleLocation:
    if start < 0:
        start, end = start + size, end + size
    return SimpleLocation(start, end, strand=strand)


def _do_pcr(
    f_match: tuple[Seq, SimpleLocation],
    r_match: tuple[Seq, SimpleLocation],
    rec: SeqRecord,
) -> SeqRecord:
    f_seq, f_loc = f_match
    r_seq, r_loc = r_match
    size = len(rec)
    start, end = int(f_loc.start), int(r_loc.end)
    if start <= int(r_loc.start) and end <= size:
        product_seq = slice_seq_record(rec, slice(start, end))
    else:
        if topology(rec) == "linear":
            raise ValueError(
                "No PCR product obtained. Maybe the template should be circular?"
            )
        product_seq = slice_circular(rec, start, end % size)

    # deal with flanking regions
    out = len(f_seq) - len(f_loc)
    if out > 0:
        product_seq = f_seq[:out] + product_seq
    out = len(r_seq) - len(r_loc)
    if out > 0:
        product_seq = product_seq + r_seq.reverse_complement()[-out:]

    return product_seq


def pcr(rec: SeqRecord, forward: str | Seq, reverse: str | Seq, min_match: int = 15):
    """Conduct PCR using `rec` as the template DNA.

    If the template is circular and the primers face outward (inverse PCR), the
    product goes across the origin.

    Parameters
    ----------
    rec : SeqRecord
        The template DNA.
    forward : str or Seq
        Sequence of forward primer
    reverse : str or Seq
        Sequence of reverse primer
    min_match : int, optional
        The minimum length of base match, by default 15
    """
    forward = Seq(forward).upper()
    reverse = Seq(reverse).upper()
    circular = topology(rec) == "circular"
    match_f = find_match(rec.seq, forward, min_match, circular=circular)
    match_r = find_match(rec.seq, reverse, min_match, circular=circular)

    if len(match_f) + len(match_r) == 0:
        raise ValueError("No PCR product obtained. No match found.")
    elif len(match_f) == 0:
        raise ValueError("No PCR product obtained. Only reverse primer matched.")
    elif len(match_r) == 0:
        raise ValueError("No PCR product obtained. Only forward primer matched.")
    elif len(match_f) > 1 or len(match_r) > 1:
        raise ValueError(
            f"Too many matches: {len(match_f)} matches found for the forward primer, "
            f"and {len(match_r)} matches found for the reverse primer."
        )
    elif match_f[0].strand == match_r[0].strand:
        raise ValueError("Each primer binds to the template in the same direction.")
    elif match_f[0].strand == 1 and match_r[0].strand == -1:
        out = _do_pcr((forward, match_f[0]), (reverse, match_r[0]), rec)
    else:
        out = _do_pcr((reverse, match_r[0]), (forward, match_f[0]), rec)

    return _as_product(out, "linear")


def gibson_assembly_single(
    seq: SeqRecord,
    overlap_range: tuple[int, int] = (15, 25),
) -> SeqRecord:
    """Simulated self-Gibson Assembly.

    Parameters
    ----------
    seq : SeqRecord
        The sequence to be assembled.

    Returns
    -------
    SeqRecord
        The product of self-Gibson Assembly.
    """
    if topology(seq) == "circular":
        raise ValueError("The input sequence must be linear DNA.")
    if len(seq) < overlap_range[0] * 2:
        raise ValueError(f"`{seq.name}` is too short.")

    overlap = _find_gibson_overlap(seq, seq, overlap_range)
    out = slice_seq_record(seq, slice(overlap, None))
    return _as_product(out, "circular")


def gibson_assembly(vec: SeqRecord, insert: SeqRecord):
    """Simulated Gibson Assembly.

    Parameters
    ----------
    vec : SeqRecord
        The (linearized) vector to assemble into.
    insert : SeqRecord
        The insert to be assembled into the vector. The 5' end of the insert must
        overlap with the 3' end of the vector, and vice versa.

    Returns
    -------
    SeqRecord
        The circular product of Gibson Assembly.
    """
    if topology(vec) == "circular" or topology(insert) == "circular":
        raise ValueError("Both vector and insert must be linear DNA.")
    if len(vec) < 30:
        raise ValueError(f"`{vec.name}` is too short.")
    if len(insert) < 30:
        raise ValueError(f"`{insert.name}` is too short.")

    ov_vec_start = _find_gibson_overlap(insert, vec)
    ov_vec_end = _find_gibson_overlap(vec, insert)
    vec_trimmed = slice_seq_record(vec, slice(ov_vec_start, len(vec) - ov_vec_end))
    return _as_product(vec_trimmed + insert, "circular")


def _find_gibson_overlap(
    left: SeqRecord,
    right: SeqRecord,
    overlap_range: tuple[int, int] = (15, 25),
) -> int:
    left_seq = left.seq.upper()
    right_seq = right.seq.upper()
    for overlap in range(overlap_range[0], overlap_range[1] + 1):
        if right_seq[:overlap] == left_seq[-overlap:]:
            break
    else:
        raise ValueError("No overlap found to perform Gibson Assembly.")
    return overlap


def self_ligation(rec: SeqRecord) -> SeqRecord:
    """Simulate self-ligation (circularization) of a linear DNA.

    The two ends are assumed to be compatible (blunt ends or phosphorylated PCR
    products such as those from inverse PCR). Features split at the two ends are
    joined into one feature across the origin.
    """
    if topology(rec) == "circular":
        raise ValueError("The input sequence is already circular.")
    if len(rec) == 0:
        raise ValueError("The input sequence is empty.")
    out = slice_seq_record(rec, slice(None))
    out.features = _join_features_at_ends(out.features, len(out))
    return _as_product(out, "circular")


def _join_features_at_ends(features: list[SeqFeature], size: int) -> list[SeqFeature]:
    """Join features that are split at the ends of a sequence being circularized."""

    # Location parts are in the biological order. On the plus strand, a feature
    # across the origin is [..., size) -> [0, ...); on the minus strand, it is
    # [0, ...) -> [..., size).
    def _touches_end(f: SeqFeature) -> bool:
        parts = f.location.parts
        part = parts[0] if parts[0].strand == -1 else parts[-1]
        return int(part.end) == size

    def _touches_start(f: SeqFeature) -> bool:
        parts = f.location.parts
        part = parts[-1] if parts[0].strand == -1 else parts[0]
        return int(part.start) == 0

    features = list(features)
    for feat_end in [f for f in features if _touches_end(f)]:
        for feat_start in features:
            if (
                feat_start is feat_end
                or not _touches_start(feat_start)
                or feat_start.type != feat_end.type
                or feat_start.qualifiers != feat_end.qualifiers
                or feat_start.location.strand != feat_end.location.strand
            ):
                continue
            if feat_end.location.parts[0].strand == -1:
                parts = feat_start.location.parts + feat_end.location.parts
            else:
                parts = feat_end.location.parts + feat_start.location.parts
            joined = copy_feature(feat_end, CompoundLocation(parts))
            features[features.index(feat_end)] = joined
            features.remove(feat_start)
            break
    return features


def _as_product(rec: SeqRecord, topo: str) -> SeqRecord:
    rec.annotations["topology"] = topo
    rec.annotations.setdefault("molecule_type", "DNA")
    return rec


def is_circular_equal(seq1: Seq, seq2: Seq) -> bool:
    """Check if two circular DNA sequences are equal."""
    if len(seq1) != len(seq2):
        return False
    return str(seq1 * 2).find(str(seq2)) >= 0


def sequencing(vec: SeqRecord, primer: str | Seq, length: int = 1000) -> SeqRecord:
    """Simulate Sanger sequencing by primer.

    Parameters
    ----------
    vec : SeqRecord
        The vector to be sequenced.
    primer : Seq
        The primer to be used for sequencing.

    Returns
    -------
    SeqRecord
        The sequenced product.
    """
    circular = topology(vec) == "circular"
    matches = find_match(vec.seq, Seq(primer).upper(), circular=circular)
    if not matches:
        raise ValueError("No match found for the primer.")
    if len(matches) > 1:
        raise ValueError("Multiple matches found for the primer.")
    m0 = matches[0]
    if m0.strand == 1:
        if circular:
            all_seq = rotate_seq_record(vec, int(m0.start))
        else:
            all_seq = slice_seq_record(vec, slice(int(m0.start), None))
    else:
        if circular:
            all_seq = rotate_seq_record(vec, int(m0.end))
        else:
            all_seq = slice_seq_record(vec, slice(None, int(m0.end)))
        all_seq = all_seq.reverse_complement(
            id=True, name=True, description=True, annotations=True, dbxrefs=True
        )
    return _as_product(slice_seq_record(all_seq, slice(0, length)), "linear")
