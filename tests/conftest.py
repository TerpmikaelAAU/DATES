import random

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

TELOMERE = "TTAGGG" * 20


def _cds(start, end, strand, note=None):
    q = {"locus_tag": [f"G{start}"]}
    if note:
        q["note"] = [note]
    return SeqFeature(FeatureLocation(start, end, strand), type="CDS", qualifiers=q)


@pytest.fixture(scope="session")
def toy_gbk(tmp_path_factory):
    """A made-up two-chromosome genome plus a small contig, in the shape a
    funannotate2/FIGS GenBank file has (notes joined with '; ')."""
    rng = random.Random(20261021)

    def seq(n, gc=0.5):
        return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT")
                       for _ in range(n))

    # chrA: telomeres at both ends, an AT-rich block in the middle
    chr_a = ("CCCTAA" * 20 + seq(40_000) + seq(20_000, gc=0.2) + seq(39_760) + TELOMERE)
    a = SeqRecord(Seq(chr_a), id="chrA", annotations={"molecule_type": "DNA"})
    a.features = [
        _cds(1_000, 2_000, 1, "antiSMASH: cluster_1 [T1PKS]; InterPro:IPR000001 Kringle"),
        _cds(3_000, 4_500, -1, "antiSMASH: cluster_1 [T1PKS]"),
        _cds(10_000, 11_000, 1, "CAZy:GH18; InterPro:IPR001223 Glycoside hydrolase"),
        SeqFeature(CompoundLocation([FeatureLocation(20_000, 20_500, -1),
                                     FeatureLocation(20_700, 21_000, -1)]),
                   type="CDS", qualifiers={"note": ["CAZy:AA7; CAZy:AA7"]}),
        _cds(50_000, 51_000, 1, "antiSMASH: cluster_2 [NRPS-like]"),
        SeqFeature(FeatureLocation(60_000, 60_072, 1), type="tRNA",
                   qualifiers={"product": ["tRNA-Met"]}),
        SeqFeature(FeatureLocation(70_000, 70_120, -1), type="rRNA",
                   qualifiers={"product": ["5S ribosomal RNA"]}),
    ]
    # chrB: a telomere at the start only
    b = SeqRecord(Seq("CCCTAA" * 20 + seq(49_880)), id="chrB",
                  annotations={"molecule_type": "DNA"})
    b.features = [_cds(5_000, 6_000, -1), _cds(7_000, 9_000, 1, "antiSMASH: cluster_3 [terpene]")]
    # a small contig with nothing on it
    c = SeqRecord(Seq(seq(5_000, gc=0.2)), id="small", annotations={"molecule_type": "DNA"})

    path = tmp_path_factory.mktemp("toy") / "toy.gbk"
    SeqIO.write([a, b, c], path, "genbank")
    return path
