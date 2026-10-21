import pytest
from Bio.Seq import reverse_complement

from dates.genome import bgc_class, find_telomere, gc_windows, has_telomere, read_genome


@pytest.fixture(scope="module")
def genome(toy_gbk):
    return {c.id: c for c in read_genome(toy_gbk)}


def test_records_in_file_order(toy_gbk):
    assert [c.id for c in read_genome(toy_gbk)] == ["chrA", "chrB", "small"]


def test_cds_strands(genome):
    a = genome["chrA"]
    assert len(a.cds) == 5
    assert sorted(s.strand for s in a.cds) == [-1, -1, 1, 1, 1]


def test_rna_genes(genome):
    a = genome["chrA"]
    assert [(s.start, s.end, s.label) for s in a.trna] == [(60_000, 60_072, "tRNA-Met")]
    assert [(s.start, s.end, s.label) for s in a.rrna] == [(70_000, 70_120, "5S ribosomal RNA")]


def test_bgc_spans_first_to_last_gene(genome):
    a = genome["chrA"]
    assert [(b.start, b.end, b.label) for b in a.bgcs] == [
        (1_000, 4_500, "T1PKS"), (50_000, 51_000, "NRPS-like")]
    assert [b.label for b in genome["chrB"].bgcs] == ["terpene"]


def test_cazymes_one_per_cds(genome):
    a = genome["chrA"]
    assert [(c.start, c.label) for c in a.cazymes] == [(10_000, "GH18"), (20_000, "AA7")]


def test_telomeres(genome):
    assert (genome["chrA"].telomere_start, genome["chrA"].telomere_end) == (0, 0)
    assert (genome["chrB"].telomere_start, genome["chrB"].telomere_end) == (0, None)
    assert (genome["small"].telomere_start, genome["small"].telomere_end) == (None, None)


def test_has_telomere_either_strand():
    assert has_telomere("TTAGGG" * 10, "TTAGGG", 10)
    assert has_telomere("CCCTAA" * 10, "TTAGGG", 10)
    assert not has_telomere("TTAGGG" * 9, "TTAGGG", 10)


@pytest.mark.parametrize("bgc_type, expected", [
    ("NRPS", "NRPS"), ("NRPS-like", "NRPS"), ("NRP-metallophore", "NRPS"),
    ("T1PKS", "PKS"), ("T3PKS", "PKS"), ("terpene", "terpene"),
    ("terpene-precursor", "terpene"), ("indole", "other"), ("", "other"),
    ("NRPS,T1PKS", "hybrid"),
])
def test_bgc_class(bgc_type, expected):
    assert bgc_class(bgc_type) == expected


def test_gc_windows():
    windows = gc_windows("GGCC" * 25 + "AATT" * 25 + "GA", 100)
    assert windows == [(50.0, 1.0), (150.0, 0.0), (201.0, 0.5)]


def test_empty_file(tmp_path):
    empty = tmp_path / "empty.gbk"
    empty.write_text("")
    with pytest.raises(ValueError):
        read_genome(empty)


def test_find_telomere_away_from_the_end():
    # an inverted copy of the subtelomere folded back past the telomere
    sub = "ACGTTGCA" * 1000  # 8 kb
    end = sub + "TTAGGG" * 20 + "CCCTAA" * 20 + sub[::-1]
    assert find_telomere(end, "TTAGGG", 1000, 10, from_end=True) == 8000
    assert find_telomere(reverse_complement(end), "TTAGGG", 1000, 10) == 8000
    assert find_telomere(sub * 3, "TTAGGG", 1000, 10) is None
    # shorter than the window: the whole sequence is searched
    assert find_telomere("TTAGGG" * 10, "TTAGGG", 1000, 10) == 0
