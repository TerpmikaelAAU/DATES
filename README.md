# DATES

**D**iagrams of **A**nnotated **T**2T **E**ukaryotic **S**equences.

Publication plots of annotated fungal genomes, made straight from the
GenBank file an annotation pipeline writes. Its companion pipelines are
[T2T_fungal_genome_pipeline](https://github.com/TerpmikaelAAU/T2T_fungal_genome_pipeline)
(assembly) and [FIGS](https://github.com/TerpmikaelAAU/FIGS) (annotation).

## Circos plot

```bash
dates circos genome.gbk -o genome_circos.png --title "Genus species STRAIN" --min-length 100000
```

From the outside in:

| ring | what | read from |
|---|---|---|
| chromosomes | length, with telomeres in orange | the sequence: copies of `TTAGGG` (either strand) near each end |
| `cds` | forward (red) and reverse (blue) strand CDSs | `CDS` features |
| `rna` | tRNA (pink) and rRNA (orange) genes | `tRNA` / `rRNA` features |
| `bgc` | antiSMASH gene clusters, coloured NRPS / PKS / hybrid / terpene / other | `antiSMASH: cluster_N [type]` CDS notes, drawn first to last gene |
| `cazy` | CAZymes | `CAZy:` CDS notes |
| `gc` | GC content per 50 kb, above (green) or below (red) the genome mean | the sequence |

It writes PNG, SVG or PDF, depending on the extension of `-o`, and prints a
short report: which chromosome ends have a telomere, what was left out, and
how many features were drawn.

Telomeres are looked for up to 20 kb in from each end, not only at the very
end, since an assembly can carry sequence past a telomere: for example, an
inverted copy of the subtelomere folded back at the end of a contig. The
report gives the distance for such a telomere (`end (2.4 kb in)`). It is
still drawn at the chromosome end, but the sequence past it may be worth
trimming from the assembly.

```
genome_circos.png: 7 sequences, 36.0 Mb
  chr_1      contig_1                7.20 Mb  telomere: start + end
  ...
  chr_6      contig_6                3.90 Mb  telomere: start + end (2.4 kb in)
  chr_7      contig_7                3.10 Mb  telomere: start
  left out: contig_8
  drawn: 12000 cds, 150 trna, 40 rrna, 50 bgc, 400 cazy, 13 telomeres
```

(made-up numbers)

Useful options (`dates circos --help` lists them all):

| option | default | |
|---|---|---|
| `--tracks` | `cds,rna,bgc,cazy,gc` | which rings to draw, outside in |
| `--min-length` | 0 | leave out shorter sequences, e.g. the mitochondrion |
| `--exclude ID ...` | | leave out these records |
| `--labels` | `number` | `chr_1, chr_2, ...` in file order, or `id` for the record IDs |
| `--telomere-repeat` | `TTAGGG` | the telomeric repeat unit |
| `--telomere-min-repeats` | 10 | copies within a `--telomere-window` (1000 bp) that make a telomere |
| `--telomere-search` | 20000 | how far in from each end (bp) to look for a telomere |
| `--gc-window` | 50000 | GC window (bp) |
| `--bgc-single-colour` | | one colour for all BGCs |

Any whole-genome GenBank file works. Functional rings stay empty for a file
without antiSMASH or CAZy notes, e.g. one straight from NCBI.

## Installation

```bash
pip install git+https://github.com/TerpmikaelAAU/DATES.git
```

or, from a clone, `uv run dates ...` (dependencies: Biopython, pyCirclize,
matplotlib). Tests: `uv run --extra test pytest`.

## Data

This repository is public and holds no genome data: the tests build a
made-up genome, and `.gitignore` keeps sequence, annotation and image files
out of commits. Keep plots of unpublished genomes outside the clone.
