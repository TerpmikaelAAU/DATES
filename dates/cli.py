"""Command line interface.

    dates circos  -- circos plot of an annotated genome GenBank file

The output format follows the extension of -o (.png, .svg, .pdf).
"""
from __future__ import annotations

import argparse
import sys

from . import __version__
from .circos import TRACKS, CircosOptions, plot_circos
from .genome import read_genome


def _select(chromosomes, min_length: int, exclude: list[str]):
    kept = [c for c in chromosomes if c.length >= min_length and c.id not in exclude]
    dropped = [c.id for c in chromosomes if c not in kept]
    if not kept:
        sys.exit("ERROR: no chromosomes left to plot after --min-length/--exclude")
    return kept, dropped


def _telomere(end: str, distance: int, window: int) -> str:
    """'end', or 'end (2.4 kb in)' for a telomere away from the very end."""
    return end if distance < window else f"{end} ({distance / 1e3:.1f} kb in)"


def cmd_circos(args) -> None:
    chromosomes = read_genome(args.genbank, telomere_repeat=args.telomere_repeat,
                              telomere_window=args.telomere_window,
                              telomere_min_repeats=args.telomere_min_repeats,
                              telomere_search=args.telomere_search)
    chromosomes, dropped = _select(chromosomes, args.min_length, args.exclude)
    if args.labels == "number":
        labels = [f"{args.label_prefix}{i}" for i in range(1, len(chromosomes) + 1)]
    else:
        labels = [c.id for c in chromosomes]
    tracks = tuple(t.strip() for t in args.tracks.split(",") if t.strip())
    options = CircosOptions(title=args.title, tracks=tracks, labels=labels,
                            gc_window=args.gc_window, bgc_by_type=not args.bgc_single_colour,
                            legend=not args.no_legend, figsize=args.size)
    try:
        counts = plot_circos(chromosomes, args.out, options, dpi=args.dpi)
    except ValueError as err:
        sys.exit(f"ERROR: {err}")

    print(f"{args.out}: {len(chromosomes)} sequences, "
          f"{sum(c.length for c in chromosomes) / 1e6:.1f} Mb")
    for chrom, label in zip(chromosomes, labels):
        ends = [_telomere(end, d, args.telomere_window)
                for end, d in (("start", chrom.telomere_start), ("end", chrom.telomere_end))
                if d is not None]
        print(f"  {label:10} {chrom.id:20} {chrom.length / 1e6:6.2f} Mb  "
              f"telomere: {' + '.join(ends) or 'none'}")
    if dropped:
        print(f"  left out: {', '.join(dropped)}")
    print("  drawn: " + ", ".join(f"{n} {k}" for k, n in counts.items()))


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="dates", description=__doc__.split("\n\n")[0])
    ap.add_argument("--version", action="version", version=f"dates {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    c = sub.add_parser("circos", help="circos plot of an annotated genome GenBank file",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    c.add_argument("genbank", help="whole-genome GenBank file (funannotate, funannotate2, FIGS)")
    c.add_argument("-o", "--out", required=True, help="output plot (.png, .svg or .pdf)")
    c.add_argument("--title", default="", help="text in the middle, e.g. the species and strain")
    c.add_argument("--tracks", default=",".join(TRACKS),
                   help=f"tracks to draw, outside in, from: {', '.join(TRACKS)}")
    c.add_argument("--labels", choices=("number", "id"), default="number",
                   help="name chromosomes by number (chr_1, chr_2, ...) in file order, "
                        "or by their GenBank record ID")
    c.add_argument("--label-prefix", default="chr_", help="prefix for --labels number")
    c.add_argument("--min-length", type=int, default=0,
                   help="leave out sequences shorter than this (bp), e.g. the mitochondrion")
    c.add_argument("--exclude", nargs="+", default=[], metavar="ID",
                   help="record IDs to leave out")
    c.add_argument("--gc-window", type=int, default=50_000, help="GC content window (bp)")
    c.add_argument("--telomere-repeat", default="TTAGGG",
                   help="telomeric repeat unit (either strand is counted)")
    c.add_argument("--telomere-window", type=int, default=1000,
                   help="window (bp) that must hold --telomere-min-repeats copies")
    c.add_argument("--telomere-search", type=int, default=20_000,
                   help="how far in from each chromosome end to look for that window (bp)")
    c.add_argument("--telomere-min-repeats", type=int, default=10,
                   help="copies of the repeat in that window that make a telomere")
    c.add_argument("--bgc-single-colour", action="store_true",
                   help="draw all BGCs in one colour instead of by type")
    c.add_argument("--no-legend", action="store_true", help="leave out the legend")
    c.add_argument("--size", type=float, default=10, help="figure width and height (inches)")
    c.add_argument("--dpi", type=int, default=300)
    c.set_defaults(func=cmd_circos)
    return ap


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
