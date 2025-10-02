from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SimpleLocation
from Bio.Seq import Seq
from himena_bio._utils import topology, slice_seq_record


def _find_all_match(vec: Seq, primer: Seq):
    return (
        list(pos for pos, _ in Seq(vec).search([primer])),
        list(pos for pos, _ in Seq(vec).search([primer.reverse_complement()])),
    )


def find_match(vec: Seq, seq: Seq, min_match: int = 15) -> list["SimpleLocation"]:
    r"""Find all the primer binding sites.

    _____ full match  ... OK
    ____/ flanking region contained ... OK
    __/\_ mismatch ... NG

    Parameters
    ----------
    seq : str or DNA
        Sequence of primer.
    min_match : int, optional
        The minimun length of match, by default 15.

    Returns
    -------
    list of DNAmatch objects
    """
    if min_match <= 0:
        raise ValueError("`min_match` must be positive value")

    matches: list["SimpleLocation"] = []
    fw_pos, rv_pos = _find_all_match(vec, seq[-min_match:])

    if min_match > len(seq):
        min_match = len(seq)

    # 1: forward check
    for pos in fw_pos:
        prpos = len(seq) - min_match  # position on seq
        while pos > 0 and prpos > 0 and vec[pos - 1] == seq[prpos - 1]:
            pos -= 1
            prpos -= 1

        matches.append(SimpleLocation(pos, len(seq) - prpos + pos, strand=1))

    # 2: reverse check
    seq_rc = seq.reverse_complement()
    for pos in rv_pos:
        prpos = min_match  # position on seq
        pos3 = pos + min_match
        while (
            pos3 < len(vec)
            and prpos < len(seq)
            and pos3 > min_match - 1
            and prpos > 0
            and vec[pos3] == seq_rc[prpos]
        ):
            pos3 += 1
            prpos += 1

        matches.append(SimpleLocation(pos, pos3, strand=-1))

    return matches


def _do_pcr(
    f_match: tuple["Seq", "SimpleLocation"],
    r_match: tuple["Seq", "SimpleLocation"],
    rec: SeqRecord,
) -> SeqRecord:
    f_seq, f_loc = f_match
    r_seq, r_loc = r_match
    if f_loc.start < r_loc.start:
        product_seq = slice_seq_record(rec, slice(f_loc.start, r_loc.end))
    else:
        if topology(rec) == "linear":
            raise ValueError("No PCR product obtained.")
        product_seq = slice_seq_record(
            rec, slice(f_loc.start, None)
        ) + slice_seq_record(rec, slice(None, r_loc.end))

    # deal with flanking regions
    out = len(f_seq) - len(f_loc)
    if out > 0:
        product_seq = f_seq[:out] + product_seq
    out = len(r_seq) - len(r_loc)
    if out > 0:
        product_seq = product_seq + r_seq.reverse_complement()[-out:]

    return product_seq


def pcr(self: SeqRecord, forward: str | Seq, reverse: str | Seq, min_match: int = 15):
    """Conduct PCR using 'self' as the template DNA.

    Parameters
    ----------
    forward : str or DNA
        Sequence of forward primer
    reverse : str or DNA
        Sequence of reverse primer
    min_match : int, optional
        The minimum length of base match, by default 15
    """
    forward = Seq(forward).upper()
    reverse = Seq(reverse).upper()
    seq = self.seq.upper()
    match_f = find_match(seq, forward, min_match)
    match_r = find_match(seq, reverse, min_match)

    # print result
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
        out = _do_pcr((forward, match_f[0]), (reverse, match_r[0]), self)
    else:
        out = _do_pcr((reverse, match_r[0]), (forward, match_f[0]), self)

    out.annotations["topology"] = "linear"
    return out


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
    num = len(seq) // 2
    return slice_seq_record(seq, slice(num, None)) + slice_seq_record(
        seq, slice(overlap, num)
    )


def gibson_assembly(vec: SeqRecord, insert: SeqRecord | None = None):
    """Simulated Gibson Assembly.

    Parameters
    ----------
    vec : SeqRecord
        The vector to assemble into.
    insert : SeqRecord, optional
        The insert to be assembled into the vector.

    Returns
    -------
    DNA object
        The product of In-Fusion
    """
    if topology(vec) == "circular" or topology(insert) == "circular":
        raise ValueError("Both vector and insert must be linear DNA.")
    if len(vec) < 30:
        raise ValueError(f"`{vec.name}` is too short.")
    if len(insert) < 30:
        raise ValueError("insert is too short.")

    pos0 = len(vec) // 2
    frag_l = slice_seq_record(vec, slice(None, pos0))
    frag_r = slice_seq_record(vec, slice(pos0, None))

    ov0 = _find_gibson_overlap(insert, frag_l)
    ov1 = _find_gibson_overlap(frag_r, insert)
    return (
        slice_seq_record(frag_r, slice(None, len(frag_r) - ov1))
        + insert
        + slice_seq_record(frag_l, slice(ov0, None))
    )


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
    matches = find_match(vec.seq, Seq(primer))
    if not matches:
        raise ValueError("No match found for the primer.")
    if len(matches) > 1:
        raise ValueError("Multiple matches found for the primer.")
    m0 = matches[0]
    if m0.strand == 1:
        if topology(vec) == "circular":
            all_seq = slice_seq_record(vec, slice(m0.start, None)) + slice_seq_record(
                vec, slice(None, m0.start)
            )
        else:
            all_seq = slice_seq_record(vec, slice(m0.start, None))
    else:
        if topology(vec) == "circular":
            all_seq = slice_seq_record(vec, slice(m0.end, None)) + slice_seq_record(
                vec, slice(None, m0.end)
            )
        else:
            all_seq = slice_seq_record(vec, slice(None, m0.end))
        all_seq = all_seq.reverse_complement()
    return slice_seq_record(all_seq, slice(0, length))
