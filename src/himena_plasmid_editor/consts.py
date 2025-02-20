from types import SimpleNamespace

class Type(SimpleNamespace):
    SEQS = "bio-seqs"

class Keys(SimpleNamespace):
    DNA = frozenset(["A", "T", "C", "G", "N"])
    RNA = frozenset(["A", "U", "C", "G", "N"])

class ApeAnnotation(SimpleNamespace):
    FWCOLOR = "ApEinfo_fwdcolor"
    RVCOLOR = "ApEinfo_revcolor"
    COMMENT = "comment"
