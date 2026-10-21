"""Circos plot of an annotated genome (pyCirclize).

From the outside in: the chromosomes, with their length and any telomeres
in orange, then the tracks asked for, in the order given:

    cds    forward and reverse strand CDSs (two rings)
    rna    tRNA (outer half) and rRNA (inner half) genes
    bgc    antiSMASH biosynthetic gene clusters, coloured by type
    cazy   CAZymes
    gc     GC content per window, as the difference from the genome mean
"""
from __future__ import annotations

from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
from matplotlib.patches import Patch  # noqa: E402
from pycirclize import Circos  # noqa: E402

from .genome import Chromosome, bgc_class, gc_windows  # noqa: E402

TRACKS = ("cds", "rna", "bgc", "cazy", "gc")
HEIGHTS = {"cds": 12, "rna": 5, "bgc": 5, "cazy": 5, "gc": 10}  # radius units of 100
GAP = 1.5

COLOURS = {
    "fwd": "#d73027",
    "rev": "#4575b4",
    "trna": "#e7298a",
    "rrna": "#ff7f00",
    "cazy": "#762a83",
    "telomere": "orange",
    "gc_high": "#1a9850",
    "gc_low": "#a50026",
    "ideogram": "lightgrey",
    "bgc": "#1b7837",  # all BGCs, with bgc_by_type off
}
# Okabe-Ito colours: told apart with colour blindness too, and none of them
# the colour of another ring.
BGC_COLOURS = {
    "NRPS": "#0072b2",  # blue
    "PKS": "#f0e442",  # yellow
    "hybrid": "#cc79a7",  # reddish purple
    "terpene": "#009e73",  # bluish green
    "other": "#999999",  # grey
}
# A sector this small (fraction of the genome) gets no name or length labels,
# which would only overlap their neighbours.
MIN_LABELLED_FRACTION = 0.015


@dataclass
class CircosOptions:
    title: str = ""
    tracks: tuple[str, ...] = TRACKS
    labels: list[str] | None = None  # one per chromosome; default: the record IDs
    gc_window: int = 50_000
    bgc_by_type: bool = True
    space: float = 2  # degrees between chromosomes
    figsize: float = 10  # inches
    legend: bool = True


def _ideogram(sector, chrom: Chromosome, label: str, labelled: bool) -> None:
    track = sector.add_track((97, 100))
    track.axis(fc=COLOURS["ideogram"], ec="black", lw=1)
    if labelled:
        track.xticks([0, sector.size], ["0 Mb", f"{sector.size / 1e6:.1f} Mb"],
                     label_orientation="vertical", tick_length=1, label_size=7)
        sector.text(label, r=104, size=10)
    tel = max(sector.size * 0.004, 1)
    kw = dict(fc=COLOURS["telomere"], ec=COLOURS["telomere"], lw=3)
    if chrom.telomere_start is not None:
        track.rect(0, min(tel, sector.size), **kw)
    if chrom.telomere_end is not None:
        track.rect(max(sector.size - tel, 0), sector.size, **kw)


def _rects(track, spans, colour, r_lim=None, ec=None, lw=0.3) -> None:
    for s in spans:
        kw = dict(fc=colour, ec=ec or colour, lw=lw)
        if r_lim:
            kw["r_lim"] = r_lim
        track.rect(s.start, max(s.end, s.start + 1), **kw)


def plot_circos(chromosomes: list[Chromosome], out: str, options: CircosOptions | None = None,
                dpi: int = 300) -> dict[str, int]:
    """Draw the plot to `out` (format from the extension) and return the
    number of features drawn per kind."""
    options = options or CircosOptions()
    unknown = set(options.tracks) - set(TRACKS)
    if unknown:
        raise ValueError(f"unknown track(s) {sorted(unknown)}; choose from {list(TRACKS)}")
    labels = options.labels or [c.id for c in chromosomes]
    total = sum(c.length for c in chromosomes)

    circos = Circos({c.id: c.length for c in chromosomes},
                    space=0 if len(chromosomes) == 1 else options.space)
    if options.title:
        circos.text(options.title, size=14, r=8)

    # Track radii, outside in.
    radii = {}
    top = 97 - GAP
    for name in options.tracks:
        radii[name] = (top - HEIGHTS[name], top)
        top -= HEIGHTS[name] + GAP
    if top < 10:
        raise ValueError("too many tracks to fit")

    gc = {}
    if "gc" in options.tracks:
        gc = {c.id: gc_windows(c.sequence, options.gc_window) for c in chromosomes}
        values = [v for windows in gc.values() for _, v in windows]
        gc_mean = sum(values) / len(values) if values else 0
        gc_lim = max((abs(v - gc_mean) for v in values), default=0) or 1

    counts = dict.fromkeys(("cds", "trna", "rrna", "bgc", "cazy", "telomeres"), 0)
    bgc_classes_seen = set()
    for chrom, label in zip(chromosomes, labels):
        sector = circos.get_sector(chrom.id)
        _ideogram(sector, chrom, label, chrom.length / total >= MIN_LABELLED_FRACTION)
        counts["telomeres"] += (chrom.telomere_start is not None) + (chrom.telomere_end is not None)
        for name in options.tracks:
            lo, hi = radii[name]
            if name == "cds":
                mid = (lo + hi) / 2
                track = sector.add_track((lo, hi))
                fwd = [s for s in chrom.cds if s.strand >= 0]
                rev = [s for s in chrom.cds if s.strand < 0]
                _rects(track, fwd, COLOURS["fwd"], r_lim=(mid, hi), lw=0.1)
                _rects(track, rev, COLOURS["rev"], r_lim=(lo, mid), lw=0.1)
                counts["cds"] += len(chrom.cds)
            elif name == "rna":
                mid = (lo + hi) / 2
                track = sector.add_track((lo, hi))
                _rects(track, chrom.trna, COLOURS["trna"], r_lim=(mid, hi), lw=0.4)
                _rects(track, chrom.rrna, COLOURS["rrna"], r_lim=(lo, mid), lw=0.4)
                counts["trna"] += len(chrom.trna)
                counts["rrna"] += len(chrom.rrna)
            elif name == "bgc":
                track = sector.add_track((lo, hi))
                for b in chrom.bgcs:
                    if options.bgc_by_type:
                        cls = bgc_class(b.label)
                        bgc_classes_seen.add(cls)
                        colour = BGC_COLOURS[cls]
                    else:
                        colour = COLOURS["bgc"]
                    _rects(track, [b], colour, ec="black")
                counts["bgc"] += len(chrom.bgcs)
            elif name == "cazy":
                track = sector.add_track((lo, hi))
                _rects(track, chrom.cazymes, COLOURS["cazy"])
                counts["cazy"] += len(chrom.cazymes)
            elif name == "gc":
                track = sector.add_track((lo, hi))
                windows = gc[chrom.id]
                if windows:
                    xs = [x for x, _ in windows]
                    dev = [v - gc_mean for _, v in windows]
                    track.fill_between(xs, [max(d, 0) for d in dev], vmin=-gc_lim, vmax=gc_lim,
                                       color=COLOURS["gc_high"], lw=0)
                    track.fill_between(xs, [min(d, 0) for d in dev], vmin=-gc_lim, vmax=gc_lim,
                                       color=COLOURS["gc_low"], lw=0)

    fig = circos.plotfig(dpi=dpi, figsize=(options.figsize, options.figsize))
    if options.legend:
        handles = []
        if "cds" in options.tracks:
            handles += [Patch(color=COLOURS["fwd"], label="Forward CDS"),
                        Patch(color=COLOURS["rev"], label="Reverse CDS")]
        if "rna" in options.tracks:
            handles += [Patch(color=COLOURS["trna"], label="tRNA"),
                        Patch(color=COLOURS["rrna"], label="rRNA")]
        if "bgc" in options.tracks:
            if options.bgc_by_type:
                handles += [Patch(color=c, label=f"BGC: {k}") for k, c in BGC_COLOURS.items()
                            if k in bgc_classes_seen]
            else:
                handles.append(Patch(color=COLOURS["bgc"], label="BGC"))
        if "cazy" in options.tracks:
            handles.append(Patch(color=COLOURS["cazy"], label="CAZyme"))
        if "gc" in options.tracks:
            handles += [Patch(color=COLOURS["gc_high"], label="GC above mean"),
                        Patch(color=COLOURS["gc_low"], label="GC below mean")]
        handles.append(Patch(color=COLOURS["telomere"], label="Telomere"))
        circos.ax.legend(handles=handles, bbox_to_anchor=(0.5, 0.42), loc="center",
                         fontsize=8, ncol=2, frameon=False)
    fig.savefig(out, dpi=dpi)
    return counts
