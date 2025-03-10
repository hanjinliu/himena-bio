from himena_bio._func import in_fusion, pcr, is_circular_equal
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

def test_in_fusion():
    vec = SeqRecord(seq=Seq(SEQ_EGFP))
    insert = SeqRecord(seq=Seq(SEQ_EGFP[-15:] + "ATATATATAT" + SEQ_EGFP[:15]))
    out = in_fusion(vec, insert)
    assert is_circular_equal(out.seq, Seq(SEQ_EGFP + "ATATATATAT"))
