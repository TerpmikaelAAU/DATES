import pytest

from dates.circos import CircosOptions, plot_circos
from dates.cli import main
from dates.genome import read_genome


def test_plot_writes_file_and_counts(toy_gbk, tmp_path):
    out = tmp_path / "plot.png"
    counts = plot_circos(read_genome(toy_gbk), str(out), CircosOptions(title="Genus species"),
                         dpi=50)
    assert out.stat().st_size > 0
    assert counts == {"cds": 7, "trna": 1, "rrna": 1, "bgc": 3, "cazy": 2, "telomeres": 3}


@pytest.mark.parametrize("suffix", [".svg", ".pdf"])
def test_vector_formats(toy_gbk, tmp_path, suffix):
    out = tmp_path / f"plot{suffix}"
    plot_circos(read_genome(toy_gbk), str(out), CircosOptions(tracks=("cds",)), dpi=50)
    assert out.stat().st_size > 0


def test_unknown_track(toy_gbk, tmp_path):
    with pytest.raises(ValueError, match="unknown track"):
        plot_circos(read_genome(toy_gbk), str(tmp_path / "x.png"),
                    CircosOptions(tracks=("cds", "methylation")))


def test_cli_min_length_and_labels(toy_gbk, tmp_path, capsys):
    out = tmp_path / "plot.png"
    main(["circos", str(toy_gbk), "-o", str(out), "--min-length", "10000", "--dpi", "50"])
    text = capsys.readouterr().out
    assert "2 sequences" in text
    assert "chr_1      chrA" in text and "chr_2      chrB" in text
    assert "left out: small" in text
    assert "telomere: start + end" in text and "telomere: start\n" in text


def test_cli_exclude_everything(toy_gbk, tmp_path):
    with pytest.raises(SystemExit, match="no chromosomes left"):
        main(["circos", str(toy_gbk), "-o", str(tmp_path / "x.png"),
              "--exclude", "chrA", "chrB", "small"])
