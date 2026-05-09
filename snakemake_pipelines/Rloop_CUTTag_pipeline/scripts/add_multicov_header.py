#!/usr/bin/env python3
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--samples", required=True, help="Comma-separated sample names matching bedtools multicov BAM order.")
parser.add_argument("--out", required=True)
args = parser.parse_args()

samples = [s for s in args.samples.split(",") if s]
header = ["chr", "start", "end"] + samples

with open(args.out, "w") as out:
    out.write("\t".join(header) + "\n")
    for line in sys.stdin:
        if not line.strip():
            continue
        out.write(line)
