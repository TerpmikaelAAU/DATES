"""Read an annotated genome GenBank file into what the plots draw.

The GenBank file is the one an annotation pipeline writes for the whole
genome (funannotate, funannotate2, FIGS): one record per chromosome or
contig, with gene/mRNA/CDS, tRNA and rRNA features. Functional annotation
is read from the CDS notes, where funannotate writes it:

    antiSMASH: cluster_12 [NRPS]      biosynthetic gene cluster membership
    CAZy:GH18                         carbohydrate-active enzyme family

A cluster is drawn from the first to the last of its genes. Telomeres are
found in the sequence itself, as tandem copies of the telomeric repeat near a
chromosome end, so no telomere table is needed. They are looked for some way
in from the end (20 kb by default), because an assembly can carry sequence
past the telomere, e.g. an inverted copy of the subtelomere folded back at a
contig end; the report says how far in such a telomere is.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import reverse_complement

ANTISMASH_RE = re.compile(r"antiSMASH:\s*([Cc]luster_?\w+)\s*(?:\[([^\]]*)\])?")
CAZY_RE = re.compile(r"CAZy:\s*([A-Za-z]+\d*)")


@dataclass
class Span:
    start: int  # 0-based, end exclusive
    end: int
    strand: int = 0
    label: str = ""


@dataclass
class Chromosome:
    id: str
    sequence: str
    cds: list[Span] = field(default_factory=list)
    trna: list[Span] = field(default_factory=list)
    rrna: list[Span] = field(default_factory=list)
    bgcs: list[Span] = field(default_factory=list)  # label = BGC type, e.g. "NRPS"
    cazymes: list[Span] = field(default_factory=list)  # label = CAZy family
    # bp from that end to the telomere's outer edge; None: no telomere found
    telomere_start: int | None = None
    telomere_end: int | None = None

    @property
    def length(self) -> int:
        return len(self.sequence)


def bgc_class(bgc_type: str) -> str:
    """Group antiSMASH region types for colouring: NRPS, PKS, terpene, other.

    Hybrids (several types, e.g. "NRPS,T1PKS") are "hybrid"."""
    types = [t.strip() for t in re.split(r"[,;/]", bgc_type) if t.strip()]
    classes = set()
    for t in types:
        if t.startswith("NRP"):  # NRPS, NRPS-like, NRP-metallophore
            classes.add("NRPS")
        elif "PKS" in t:  # T1PKS, T3PKS, PKS-like, ...
            classes.add("PKS")
        elif t.startswith("terpene"):  # terpene, terpene-precursor
            classes.add("terpene")
        else:
            classes.add("other")
    if len(classes) > 1:
        return "hybrid"
    return classes.pop() if classes else "other"


def has_telomere(end_sequence: str, repeat: str, min_repeats: int) -> bool:
    """True if the repeat (either strand) occurs min_repeats times in it."""
    s = end_sequence.upper()
    rc = reverse_complement(repeat.upper())
    return max(s.count(repeat.upper()), s.count(rc)) >= min_repeats


def find_telomere(end_sequence: str, repeat: str, window: int, min_repeats: int,
                  from_end: bool = False) -> int | None:
    """Distance (bp) from the outer edge of `end_sequence` to the telomere
    nearest it, or None. The outer edge is the start of the sequence, or its
    end with from_end=True. The telomere is the outermost `window` holding
    `min_repeats` copies of the repeat; its edge is the outermost copy there."""
    s = end_sequence.upper()
    units = (repeat.upper(), reverse_complement(repeat.upper()))
    n = len(s)
    window = min(window, n)
    for offset in range(0, n - window + 1, max(window // 10, 1)):
        start = n - offset - window if from_end else offset
        w = s[start:start + window]
        if not has_telomere(w, repeat, min_repeats):
            continue
        if from_end:
            return n - start - max(w.rfind(u) + len(u) for u in units if u in w)
        return start + min(w.find(u) for u in units if u in w)
    return None


def gc_windows(sequence: str, window: int) -> list[tuple[float, float]]:
    """(window midpoint, GC fraction of A/C/G/T) for each window."""
    out = []
    seq = sequence.upper()
    for start in range(0, len(seq), window):
        w = seq[start:start + window]
        gc = w.count("G") + w.count("C")
        acgt = gc + w.count("A") + w.count("T")
        if acgt:
            out.append((start + len(w) / 2, gc / acgt))
    return out


def _span(feature, label: str = "") -> Span:
    loc = feature.location
    return Span(int(loc.start), int(loc.end), loc.strand or 0, label)


def read_genome(path: str | Path, telomere_repeat: str = "TTAGGG",
                telomere_window: int = 1000, telomere_min_repeats: int = 10,
                telomere_search: int = 20_000) -> list[Chromosome]:
    """Every record of a GenBank file, in file order.

    A chromosome end has a telomere if a `telomere_window` within
    `telomere_search` bp of it holds `telomere_min_repeats` copies of the
    repeat."""
    chromosomes = []
    for record in SeqIO.parse(str(path), "genbank"):
        seq = str(record.seq)
        chrom = Chromosome(record.id, seq)
        clusters: dict[str, list] = {}
        for f in record.features:
            if f.type == "CDS":
                chrom.cds.append(_span(f))
                note = " ".join(f.qualifiers.get("note", []))
                for m in ANTISMASH_RE.finditer(note):
                    cid = m.group(1)
                    entry = clusters.setdefault(cid, [int(f.location.start), int(f.location.end), ""])
                    entry[0] = min(entry[0], int(f.location.start))
                    entry[1] = max(entry[1], int(f.location.end))
                    entry[2] = entry[2] or (m.group(2) or "")
                families = sorted(set(CAZY_RE.findall(note)))
                if families:
                    chrom.cazymes.append(_span(f, ",".join(families)))
            elif f.type == "tRNA":
                chrom.trna.append(_span(f, f.qualifiers.get("product", [""])[0]))
            elif f.type == "rRNA":
                chrom.rrna.append(_span(f, f.qualifiers.get("product", [""])[0]))
        chrom.bgcs = [Span(s, e, 0, t) for s, e, t in
                      sorted(clusters.values(), key=lambda c: c[0])]
        if seq:
            search = max(telomere_search, telomere_window)
            chrom.telomere_start = find_telomere(seq[:search], telomere_repeat,
                                                 telomere_window, telomere_min_repeats)
            chrom.telomere_end = find_telomere(seq[-search:], telomere_repeat,
                                               telomere_window, telomere_min_repeats,
                                               from_end=True)
        chromosomes.append(chrom)
    if not chromosomes:
        raise ValueError(f"{path}: no GenBank records")
    return chromosomes
