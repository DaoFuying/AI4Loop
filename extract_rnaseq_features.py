#!/usr/bin/env python3
"""Extract multi-scale RNA-seq features for AI4Loop candidate gene pairs.

This script replaces the hard-coded callRNAseq.py workflow with a command-line interface.
It supports coordinate files with columns: label, chr1, x1, x2, chr2, y1, y2.
"""

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


def normalize_chr(value) -> str:
    value = str(value).replace(".0", "")
    if value.startswith("chr"):
        return value
    if value in {"23", "X"}:
        return "chrX"
    if value in {"24", "Y"}:
        return "chrY"
    return f"chr{value}"


def read_rnaseq(path: str, sample_columns: List[str] = None) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    df.columns = df.columns.astype(str).str.replace("\ufeff", "", regex=False).str.strip()

    # Harmonize common first-column names.
    rename = {}
    for c in df.columns:
        if c in {"#chr", "chrom", "Chromosome"}:
            rename[c] = "chr"
        elif c.lower() == "start":
            rename[c] = "start"
        elif c.lower() == "end":
            rename[c] = "end"
    df = df.rename(columns=rename)

    required = {"chr", "start", "end"}
    if not required.issubset(set(df.columns)):
        # If no header was detected, assume first columns are chr/start/end/gene.
        if df.shape[1] >= 5:
            cols = list(df.columns)
            df = pd.read_csv(path, sep=None, engine="python", header=None)
            df.columns = ["chr", "start", "end", "gene"] + [f"sample_{i}" for i in range(1, df.shape[1] - 3)]
        else:
            raise ValueError(f"RNA-seq file must contain chr/start/end columns. Found: {list(df.columns)}")

    df["chr"] = df["chr"].map(normalize_chr)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    if sample_columns is None:
        sample_columns = [c for c in df.columns if c not in {"chr", "start", "end", "gene", "gene_id", "gene_name"}]
    missing = [c for c in sample_columns if c not in df.columns]
    if missing:
        raise ValueError(f"Sample columns not found in RNA-seq file: {missing}")

    for c in sample_columns:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    return df[["chr", "start", "end"] + sample_columns]


def read_pairs(path: str) -> pd.DataFrame:
    pairs = pd.read_csv(path, sep=None, engine="python")
    pairs.columns = pairs.columns.astype(str).str.replace("\ufeff", "", regex=False).str.strip()
    required = ["chr1", "x1", "x2", "chr2", "y1", "y2"]
    missing = [c for c in required if c not in pairs.columns]
    if missing:
        raise ValueError(f"Pair file must contain {required}. Missing: {missing}; found: {list(pairs.columns)}")
    pairs["chr1"] = pairs["chr1"].map(normalize_chr)
    pairs["chr2"] = pairs["chr2"].map(normalize_chr)
    for c in ["x1", "x2", "y1", "y2"]:
        pairs[c] = pairs[c].astype(int)
    return pairs


def group_genes_by_chrom(rnaseq: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    return {chrom: sub.sort_values("start").reset_index(drop=True) for chrom, sub in rnaseq.groupby("chr", sort=False)}


def mean_overlap_expression(genes: pd.DataFrame, chrom: str, start: int, end: int, sample_cols: List[str]) -> np.ndarray:
    if chrom not in genes_by_chr_global:
        return np.zeros(len(sample_cols), dtype=float)
    sub = genes_by_chr_global[chrom]
    # Overlap: gene.start < bin.end and gene.end > bin.start
    hit = sub[(sub["start"] < end) & (sub["end"] > start)]
    if hit.empty:
        return np.zeros(len(sample_cols), dtype=float)
    return hit[sample_cols].mean(axis=0).values.astype(float)


def anchor_windows(chrom: str, start: int, end: int, flank: int) -> Tuple[str, int, int]:
    center = int(start + (end - start) / 2)
    return chrom, max(0, center - flank), center + flank


def extract_features(pairs: pd.DataFrame, rnaseq: pd.DataFrame, sample: str, bins: Iterable[int], flank: int) -> pd.DataFrame:
    global genes_by_chr_global
    genes_by_chr_global = group_genes_by_chrom(rnaseq)
    sample_cols = [sample]

    features = {}
    for win in bins:
        region_len = 2 * flank
        n_bins = int(region_len / win)
        for side, chrom_col, start_col, end_col, prefix in [
            ("left", "chr1", "x1", "x2", "L"),
            ("right", "chr2", "y1", "y2", "R"),
        ]:
            anchors = [anchor_windows(row[chrom_col], row[start_col], row[end_col], flank) for _, row in pairs.iterrows()]
            for i in range(n_bins):
                values = []
                for chrom, a_start, _ in anchors:
                    b_start = a_start + i * win
                    b_end = a_start + (i + 1) * win
                    values.append(mean_overlap_expression(genes_by_chr_global, chrom, b_start, b_end, sample_cols)[0])
                features[f"{prefix}{i}_{win}"] = values

    feature_df = pd.DataFrame(features)
    # Preserve metadata for evaluation and reproducibility.
    for col in ["label", "chr1", "x1", "x2", "chr2", "y1", "y2"]:
        if col in pairs.columns:
            feature_df.insert(0 if col == "label" else len(feature_df.columns), col, pairs[col].values)
    return feature_df


def main():
    parser = argparse.ArgumentParser(description="Extract AI4Loop multi-scale RNA-seq features from coordinate pairs.")
    parser.add_argument("--pairs", required=True, help="Candidate pair CSV with chr1/x1/x2/chr2/y1/y2 and optional label")
    parser.add_argument("--rnaseq_bed", required=True, help="RNA-seq BED/TSV file with chr/start/end and expression columns")
    parser.add_argument("--sample", required=True, help="RNA-seq sample column to use, e.g. K562")
    parser.add_argument("--output", required=True, help="Output feature CSV")
    parser.add_argument("--bins", nargs="+", type=int, default=[1000, 2000, 3000, 4000, 5000, 6000, 7000])
    parser.add_argument("--flank", type=int, default=30000)
    args = parser.parse_args()

    pairs = read_pairs(args.pairs)
    rnaseq = read_rnaseq(args.rnaseq_bed, [args.sample])
    out = extract_features(pairs, rnaseq, args.sample, args.bins, args.flank)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"Saved {out.shape[0]} rows and {out.shape[1]} columns to {args.output}")


if __name__ == "__main__":
    genes_by_chr_global = {}
    main()
