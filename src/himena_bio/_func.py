from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SimpleLocation
from Bio.Seq import Seq
from himena_bio._utils import topology


def _find_all_match(vec: Seq, primer: Seq):
    return (
        list(pos for pos, _ in Seq(vec).search([primer])),
        list(pos for pos, _ in Seq(vec).search([primer.reverse_complement()])),
    )


def find_match(vec: Seq, seq: Seq, min_match: int = 15) -> list["SimpleLocation"]:
    r"""
    Find all the primer binding sites.
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
):
    # def _do_pcr(match: tuple["SimpleLocation", "SimpleLocation"], rec: SeqRecord):
    f_seq, f_loc = f_match
    r_seq, r_loc = r_match
    if f_loc.start > r_loc.start:
        f_loc, r_loc = r_loc, f_loc
        f_seq, r_seq = r_seq, f_seq

    if f_loc.strand == 1 and r_loc.strand == -1:
        mf, mr = f_loc, r_loc
        product_seq = rec[mf.start : mr.end]

    elif topology(rec) == "linear":
        raise ValueError("No PCR product obtained.")

    else:
        mr, mf = r_loc, f_loc
        product_seq = (rec + rec)[mf.start : mr.end + len(rec)]

    # deal with flanking regions
    out = len(f_seq) - len(mf)
    if out > 0:
        product_seq = f_seq[:out] + product_seq
    out = len(r_seq) - len(mr)
    if out > 0:
        product_seq = product_seq + r_seq.reverse_complement()[-out:]
        # product_seq = product_seq + mr.ref.r[-out:]

    return product_seq


def pcr(self: SeqRecord, forward: str | Seq, reverse: str | Seq, min_match: int = 15):
    """
    Conduct PCR using 'self' as the template DNA.

    Parameters
    ----------
    forward : str or DNA
        Sequence of forward primer
    reverse : str or DNA
        Sequence of reverse primer
    min_match : int, optional
        The minimum length of base match, by default 15

    Returns
    -------
    DNA
        The PCR product.
    """
    forward = Seq(forward)
    reverse = Seq(reverse)
    matches = find_match(self.seq, forward, min_match) + find_match(
        self.seq, reverse, min_match
    )

    # print result
    if len(matches) == 0:
        raise ValueError("No PCR product obtained. No match found.")
    elif len(matches) == 1:
        raise ValueError("No PCR product obtained. Only one match found.")
    elif len(matches) > 2:
        raise ValueError(f"Too many matches ({len(matches)} matches found).")
    elif matches[0].strand == matches[1].strand:
        raise ValueError("Each primer binds to the template in the same direction.")
    else:
        ans = _do_pcr((forward, matches[0]), (reverse, matches[1]), self)

    return ans


def in_fusion(vec: SeqRecord, insert: SeqRecord):
    """
    Simulated In-Fusion.

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

    if frag_l[:15] != insert[-15:]:
        raise ValueError(
            "Mismatch! Check carefully:\n"
            f"--{insert[-20:]}\n"
            f"       {frag_l[:20]}--"
        )
    if frag_r[-15:] != insert[:15]:
        raise ValueError(
            "Mismatch! Check carefully:\n"
            f"       {insert[:20]}--\n"
            f"--{frag_r[-20:]}"
        )

    frag_r = frag_r[: len(frag_r) - 15]
    frag_l = frag_l[15:]
    out = frag_r + insert + frag_l
    return out
