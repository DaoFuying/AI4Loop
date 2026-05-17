#!/usr/bin/env python3
"""Backward-compatible wrapper for AI4Loop prediction.

Preferred usage:
    python predict.py --model MODEL.h5 --gene_pairs GENE_PAIRS.bedpe --rnaseq_bed RNA.bed --column_names RNA.columns.txt --output predictions.csv

Legacy usage is still supported by running this file without arguments.
"""

import sys

from predict import main as predict_main


def legacy_main():
    if len(sys.argv) == 1:
        sys.argv = [
            sys.argv[0],
            "--model", "models/K562_RandomSplit.model.h5",
            "--gene_pairs", "prediction/Genepairs_1000.bedpe",
            "--rnaseq_bed", "prediction/Samples_RNASeqdata.bed",
            "--column_names", "prediction/Samples_RNASeqdata.bed.columns.txt",
            "--output", "prediction/Samples_RNASeqdata.bed.Pre.csv",
            "--bedpe_format", "ai4loop",
        ]
    predict_main()


if __name__ == "__main__":
    legacy_main()
