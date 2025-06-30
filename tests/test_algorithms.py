from himena_bio._func import gibson_assembly, pcr, is_circular_equal, gibson_assembly_single
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

SEQ_EGFP = "GTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGG"

@pytest.mark.parametrize(
    "template, forward, reverse, expected",
    [
        (SEQ_EGFP, "GTGAGCAAG", "CCAGGATGGGCACC", SEQ_EGFP),
        (SEQ_EGFP, "GTGAGCAAG", "CCAGGATGGGCAC", SEQ_EGFP),
        (SEQ_EGFP, "GAGCAAGGG", "CCAGGATGGGC", SEQ_EGFP[2:]),
        (SEQ_EGFP, "GAGCAAGGG", "GGATGGGCACC", SEQ_EGFP[2:-3]),
        (SEQ_EGFP, "ATATAATGTGAGCAAG", "ATTATTCCAGGATGGGCAC", f"ATATAAT{SEQ_EGFP}AATAAT"),
    ]
)
def test_pcr_linear(template: str, forward: str, reverse: str, expected: str):
    rec = SeqRecord(id="test", seq=Seq(template))
    rec.annotations["topology"] = "linear"
    out = pcr(rec, forward, reverse, min_match=8)
    assert str(out.seq) == expected

@pytest.mark.parametrize(
    "template, forward, reverse, expected",
    [
        (SEQ_EGFP, "GTGAGCAAG", "CCAGGATGGGCACC", SEQ_EGFP),
        (SEQ_EGFP, "GTGAGCAAG", "CCAGGATGGGCAC", SEQ_EGFP),
        (SEQ_EGFP, "GGTGGTGCCCA", "TCCTCGCCCTTG", "GGTGGTGCCCATCCTGGGTGAGCAAGGGCGAGGA"),
        (SEQ_EGFP, "ATATAATGGTGGTGCCCA", "ATTATTTCCTCGCCCTTG", "ATATAATGGTGGTGCCCATCCTGGGTGAGCAAGGGCGAGGAAATAAT"),
    ],
)
def test_pcr_circular(template: str, forward: str, reverse: str, expected: str):
    rec = SeqRecord(id="test", seq=Seq(template))
    rec.annotations["topology"] = "circular"
    out = pcr(rec, forward, reverse, min_match=8)
    assert str(out.seq) == expected

@pytest.mark.parametrize(
    "seq1, seq2, expected",
    [
        ("ATATGC", "ATATGC", True),
        ("ATAT", "ATATGC", False),
        ("ATATGC", "ATGCAT", True),
        ("GGCTAATTGACTCT", "ATTGACTCTGGCTA", True),
        ("GGCTAATTGACTCT", "CTGGCTAATTGACT", True),
        ("GGCTAATTGACTCT", "CTGGCTAATTGAGT", False),
    ],
)
def test_circular_equal(seq1, seq2, expected):
    assert is_circular_equal(Seq(seq1), Seq(seq2)) == expected

@pytest.mark.parametrize("ov0, ov1", [(15, 15), (16, 18), (19, 17)])
def test_gibson(ov0: int, ov1: int):
    vec = SeqRecord(seq=Seq(SEQ_EGFP))
    insert = SeqRecord(seq=Seq(SEQ_EGFP[-ov0:] + "ATATATATAT" + SEQ_EGFP[:ov1]))
    out = gibson_assembly(vec, insert)
    assert is_circular_equal(out.seq, Seq(SEQ_EGFP + "ATATATATAT"))

@pytest.mark.parametrize(
    "ov_seq",
    [
        "GTGAAGTTCCTCAGT",
        "CGTGAAGTTCCTCAGTC",
    ]
)
def test_gibson_single(ov_seq: str):
    vec = SeqRecord(seq=Seq(ov_seq + SEQ_EGFP + ov_seq))
    out = gibson_assembly_single(vec)
    assert is_circular_equal(out.seq, Seq(SEQ_EGFP + ov_seq))

def test_sanger_sequencing():
    from himena_bio._func import sequencing

    rec = SeqRecord(seq=Seq(SEQ_EGFP))
    out = sequencing(rec, SEQ_EGFP[15:25])
    assert str(out.seq) == SEQ_EGFP[15:]

    assert len(SEQ_EGFP) > 37
    out = sequencing(rec, Seq(SEQ_EGFP[27:37]).reverse_complement())
    assert str(out.seq) == str(Seq(SEQ_EGFP[:37]).reverse_complement())
