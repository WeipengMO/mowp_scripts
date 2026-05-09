#!/usr/bin/env python3
import argparse
import math
import subprocess
from pathlib import Path

import pandas as pd


def mapped_reads(bam: str) -> int:
    cmd = ["samtools", "view", "-c", "-F", "260", bam]
    return int(subprocess.check_output(cmd, text=True).strip())


def coverage_counts(peaks: str, bam: str) -> pd.Series:
    cmd = ["bedtools", "coverage", "-a", peaks, "-b", bam, "-counts"]
    out = subprocess.check_output(cmd, text=True)
    counts = []
    for line in out.splitlines():
        if not line.strip():
            continue
        counts.append(int(float(line.rstrip().split("\t")[-1])))
    return pd.Series(counts, dtype="float64")


def main():
    parser = argparse.ArgumentParser(
        description="Filter S9.6/HBD R-loop CUT&Tag candidate peaks by background and RNase H depletion."
    )
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--signal-bam", required=True)
    parser.add_argument("--background-bam", required=True)
    parser.add_argument("--rnaseh-bam", required=True)
    parser.add_argument("--out-bed", required=True)
    parser.add_argument("--out-metrics", required=True)
    parser.add_argument("--min-log2fc-background", type=float, default=1.0)
    parser.add_argument("--min-log2fc-rnaseh", type=float, default=1.0)
    parser.add_argument("--min-signal-cpm", type=float, default=0.0)
    parser.add_argument("--pseudocount", type=float, default=1.0)
    args = parser.parse_args()

    peaks = pd.read_csv(args.peaks, sep="\t", header=None, comment="#")
    if peaks.shape[0] == 0:
        Path(args.out_bed).write_text("")
        Path(args.out_metrics).write_text("")
        return

    peaks = peaks.iloc[:, :min(peaks.shape[1], 10)].copy()
    base_cols = [
        "chrom", "start", "end", "name", "score", "strand",
        "fold_enrichment", "minus_log10_p", "minus_log10_q", "summit"
    ]
    peaks.columns = base_cols[:peaks.shape[1]]

    sig_total = mapped_reads(args.signal_bam)
    bg_total = mapped_reads(args.background_bam)
    rh_total = mapped_reads(args.rnaseh_bam)

    peaks["signal_count"] = coverage_counts(args.peaks, args.signal_bam).values
    peaks["background_count"] = coverage_counts(args.peaks, args.background_bam).values
    peaks["rnaseh_count"] = coverage_counts(args.peaks, args.rnaseh_bam).values

    peaks["signal_cpm"] = peaks["signal_count"] / max(sig_total, 1) * 1e6
    peaks["background_cpm"] = peaks["background_count"] / max(bg_total, 1) * 1e6
    peaks["rnaseh_cpm"] = peaks["rnaseh_count"] / max(rh_total, 1) * 1e6

    pc = args.pseudocount
    peaks["log2fc_vs_background"] = (
        (peaks["signal_cpm"] + pc) / (peaks["background_cpm"] + pc)
    ).map(math.log2)
    peaks["log2fc_vs_rnaseh"] = (
        (peaks["signal_cpm"] + pc) / (peaks["rnaseh_cpm"] + pc)
    ).map(math.log2)

    keep = (
        (peaks["signal_cpm"] >= args.min_signal_cpm)
        & (peaks["log2fc_vs_background"] >= args.min_log2fc_background)
        & (peaks["log2fc_vs_rnaseh"] >= args.min_log2fc_rnaseh)
    )

    peaks.to_csv(args.out_metrics, sep="\t", index=False)

    out = peaks.loc[keep].copy()
    bed_cols = ["chrom", "start", "end"]
    if "name" in out.columns:
        bed_cols.append("name")
    bed_cols += [
        "signal_cpm", "background_cpm", "rnaseh_cpm",
        "log2fc_vs_background", "log2fc_vs_rnaseh"
    ]
    out[bed_cols].to_csv(args.out_bed, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
