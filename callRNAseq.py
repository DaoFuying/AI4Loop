#!/usr/bin/env python3
"""Backward-compatible wrapper for RNA-seq feature extraction.

Preferred usage:
    python extract_rnaseq_features.py --pairs PAIRS.csv --rnaseq_bed RNA.tsv --sample K562 --output FEATURES.csv

Legacy usage is still supported:
    python callRNAseq.py k562_ctcf out_dir
"""

import sys
from pathlib import Path

from extract_rnaseq_features import main as extract_main


def legacy_main():
    if len(sys.argv) == 3 and not sys.argv[1].startswith("--"):
        eachcell = sys.argv[1]
        out_dir = sys.argv[2]
        sample = eachcell.split("_")[0].upper()
        if sample == "HELAS3":
            sample = "HeLaS3"
        pair_file = Path(out_dir) / f"{eachcell}_distance_matched.csv"
        out_file = Path(out_dir) / f"{eachcell}_distance_matched.csv_winGEfea.csv"
        sys.argv = [
            sys.argv[0],
            "--pairs", str(pair_file),
            "--rnaseq_bed", "data/allRNAseq.tsv",
            "--sample", sample,
            "--output", str(out_file),
        ]
    extract_main()


if __name__ == "__main__":
    legacy_main()
