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
        (SEQ_EGFP, "GAGcAAGgG", "CcAGGaTGGGC", SEQ_EGFP[2:]),
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

def _circular(seq: str) -> SeqRecord:
    rec = SeqRecord(id="test", seq=Seq(seq))
    rec.annotations["topology"] = "circular"
    return rec

def test_inverse_pcr_features():
    from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation

    rec = _circular(SEQ_EGFP)
    # feature across the origin
    rec.features.append(
        SeqFeature(
            CompoundLocation([SimpleLocation(40, 46, 1), SimpleLocation(0, 5, 1)]),
            type="wrap",
        )
    )
    # feature in the deleted region
    rec.features.append(SeqFeature(SimpleLocation(20, 25, -1), type="deleted"))
    # feature at the edge of the product
    rec.features.append(SeqFeature(SimpleLocation(10, 20, -1), type="edge"))
    out = pcr(rec, "GGTGGTGCCCA", "TCCTCGCCCTTG", min_match=8)
    assert str(out.seq) == "GGTGGTGCCCATCCTGGGTGAGCAAGGGCGAGGA"
    feats = {f.type: f for f in out.features}
    assert set(feats) == {"wrap", "edge"}
    assert (feats["wrap"].location.start, feats["wrap"].location.end) == (11, 22)
    assert feats["wrap"].location.strand == 1
    assert (feats["edge"].location.start, feats["edge"].location.end) == (27, 34)
    assert feats["edge"].location.strand == -1

@pytest.mark.parametrize("shift", [0, 5, 20, 40, 45])
def test_pcr_primer_across_origin(shift: int):
    template = SEQ_EGFP[shift:] + SEQ_EGFP[:shift]
    out = pcr(_circular(template), "GTGAGCAAG", "CCAGGATGGGCAC", min_match=8)
    assert str(out.seq) == SEQ_EGFP

def test_self_ligation():
    from himena_bio._func import self_ligation
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    rec = SeqRecord(seq=Seq(SEQ_EGFP))
    rec.annotations["topology"] = "linear"
    n = len(SEQ_EGFP)
    rec.features = [
        SeqFeature(SimpleLocation(n - 5, n, 1), type="f", qualifiers={"label": ["a"]}),
        SeqFeature(SimpleLocation(0, 3, 1), type="f", qualifiers={"label": ["a"]}),
        SeqFeature(SimpleLocation(0, 3, 1), type="f", qualifiers={"label": ["b"]}),
    ]
    out = self_ligation(rec)
    assert out.annotations["topology"] == "circular"
    assert str(out.seq) == SEQ_EGFP
    assert len(out.features) == 2
    joined = out.features[0].location
    assert [(p.start, p.end) for p in joined.parts] == [(n - 5, n), (0, 3)]
    with pytest.raises(ValueError):
        self_ligation(out)

def test_inverse_pcr_then_self_ligation():
    from himena_bio._func import self_ligation

    rec = _circular(SEQ_EGFP)
    product = pcr(rec, "GGTGGTGCCCA", "TCCTCGCCCTTG", min_match=8)
    out = self_ligation(product)
    # region between the primers are deleted
    assert is_circular_equal(out.seq, Seq(SEQ_EGFP[:17] + SEQ_EGFP[29:]))

def test_sequencing_circular_lowercase():
    from himena_bio._func import sequencing

    rec = _circular(SEQ_EGFP.lower())
    out = sequencing(rec, SEQ_EGFP[30:40], length=20)
    assert str(out.seq).upper() == (SEQ_EGFP[30:] + SEQ_EGFP)[:20]
    out = sequencing(rec, Seq(SEQ_EGFP[5:15]).reverse_complement(), length=20)
    expected = Seq(SEQ_EGFP + SEQ_EGFP[:15]).reverse_complement()[:20]
    assert str(out.seq).upper() == str(expected)
