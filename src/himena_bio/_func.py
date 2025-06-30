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
    forward = Seq(forward)
    reverse = Seq(reverse)
    match_f = find_match(self.seq, forward, min_match)
    match_r = find_match(self.seq, reverse, min_match)

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


def in_fusion(vec: SeqRecord, insert: SeqRecord):
    """Simulated In-Fusion.

    Parameters
    ----------
    insert : str or DNA
        The DNA fragment to insert.

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
    frag_l, frag_r = vec[:pos0], vec[pos0:]

    if frag_l[:15].seq != insert[-15:].seq:
        raise ValueError(
            "Mismatch! Check carefully:\n"
            f"--{insert[-20:].seq}\n"
            f"       {frag_l[:20].seq}--"
        )
    if frag_r[-15:].seq != insert[:15].seq:
        raise ValueError(
            "Mismatch! Check carefully:\n"
            f"       {insert[:20].seq}--\n"
            f"--{frag_r[-20:].seq}"
        )

    frag_r = frag_r[: len(frag_r) - 15]
    frag_l = frag_l[15:]
    out = frag_r + insert + frag_l
    return out


def is_circular_equal(seq1: Seq, seq2: Seq) -> bool:
    """Check if two circular DNA sequences are equal."""
    if len(seq1) != len(seq2):
        return False
    return str(seq1 * 2).find(str(seq2)) >= 0
