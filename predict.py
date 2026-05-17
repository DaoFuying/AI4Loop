#!/usr/bin/env python3
"""Predict gene-centered chromatin interactions from RNA-seq using a pretrained AI4Loop model."""

import argparse
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

from extract_rnaseq_features import extract_features, normalize_chr, read_rnaseq


def read_column_names(path: str) -> List[str]:
    names = [line.strip() for line in open(path) if line.strip()]
    return names


def read_rnaseq_for_prediction(path: str, column_names_file: str = None):
    if column_names_file:
        names = read_column_names(column_names_file)
        df = pd.read_csv(path, sep="\t", header=None)
        if len(names) != df.shape[1]:
            raise ValueError(f"Column-name file contains {len(names)} names, but RNA-seq file has {df.shape[1]} columns.")
        df.columns = names
        tmp = Path(path).with_suffix(Path(path).suffix + ".with_header.tmp.tsv")
        df.to_csv(tmp, sep="\t", index=False)
        try:
            rnaseq = read_rnaseq(str(tmp))
        finally:
            tmp.unlink(missing_ok=True)
        sample_cols = [c for c in names if c not in {"0", "1", "2", "chr", "start", "end", "gene", "gene_id", "gene_id2"}]
        # Harmonize common first-column names.
        if {"0", "1", "2"}.issubset(set(names)):
            rnaseq = df.rename(columns={"0": "chr", "1": "start", "2": "end", "gene_id2": "gene"})
            rnaseq["chr"] = rnaseq["chr"].map(normalize_chr)
            rnaseq["start"] = rnaseq["start"].astype(int)
            rnaseq["end"] = rnaseq["end"].astype(int)
            for s in sample_cols:
                rnaseq[s] = pd.to_numeric(rnaseq[s], errors="coerce").fillna(0.0)
            rnaseq = rnaseq[["chr", "start", "end"] + sample_cols]
        return rnaseq, sample_cols

    rnaseq = read_rnaseq(path)
    sample_cols = [c for c in rnaseq.columns if c not in {"chr", "start", "end"}]
    return rnaseq, sample_cols


def read_gene_pairs(path: str, bedpe_format: str = "auto"):
    raw = pd.read_csv(path, sep="\t", header=None, comment="#")
    if raw.shape[1] < 6:
        raise ValueError("Gene-pair file must contain at least six BEDPE columns.")

    if bedpe_format == "auto":
        # AI4Loop example: chr1 start1 end1 geneid1 genename1 chr2 start2 end2 geneid2 genename2 score
        if raw.shape[1] >= 8 and str(raw.iloc[0, 5]).startswith("chr"):
            bedpe_format = "ai4loop"
        else:
            bedpe_format = "standard"

    if bedpe_format == "ai4loop":
        pairs = pd.DataFrame({
            "chr1": raw.iloc[:, 0].map(normalize_chr),
            "x1": raw.iloc[:, 1].astype(int),
            "x2": raw.iloc[:, 2].astype(int),
            "chr2": raw.iloc[:, 5].map(normalize_chr),
            "y1": raw.iloc[:, 6].astype(int),
            "y2": raw.iloc[:, 7].astype(int),
        })
        metadata = raw.copy()
        metadata.columns = [f"input_col_{i}" for i in range(raw.shape[1])]
    else:
        pairs = pd.DataFrame({
            "chr1": raw.iloc[:, 0].map(normalize_chr),
            "x1": raw.iloc[:, 1].astype(int),
            "x2": raw.iloc[:, 2].astype(int),
            "chr2": raw.iloc[:, 3].map(normalize_chr),
            "y1": raw.iloc[:, 4].astype(int),
            "y2": raw.iloc[:, 5].astype(int),
        })
        metadata = raw.copy()
        metadata.columns = [f"input_col_{i}" for i in range(raw.shape[1])]

    pairs.insert(0, "label", 0)
    return pairs, metadata


def predict_feature_df(feature_df: pd.DataFrame, model) -> np.ndarray:
    left_cols = [c for c in feature_df.columns if c.startswith("L")]
    right_cols = [c for c in feature_df.columns if c.startswith("R")]
    left = feature_df[left_cols].replace([np.inf, -np.inf], np.nan).fillna(0).astype("float32").values.reshape(-1, len(left_cols), 1)
    right = feature_df[right_cols].replace([np.inf, -np.inf], np.nan).fillna(0).astype("float32").values.reshape(-1, len(right_cols), 1)
    return model.predict([left, right], verbose=0).ravel()


def main():
    parser = argparse.ArgumentParser(description="Predict GCIs from RNA-seq using a pretrained AI4Loop model.")
    parser.add_argument("--model", required=True, help="Pretrained AI4Loop .h5 model")
    parser.add_argument("--gene_pairs", required=True, help="BEDPE-like candidate gene-pair file")
    parser.add_argument("--rnaseq_bed", required=True, help="RNA-seq BED-like file")
    parser.add_argument("--column_names", default=None, help="Optional one-column-per-line name file for RNA-seq BED columns")
    parser.add_argument("--output", required=True, help="Output CSV containing prediction scores")
    parser.add_argument("--bedpe_format", choices=["auto", "standard", "ai4loop"], default="auto")
    parser.add_argument("--samples", nargs="+", default=None, help="Optional subset of RNA-seq sample columns")
    parser.add_argument("--bins", nargs="+", type=int, default=[1000, 2000, 3000, 4000, 5000, 6000, 7000])
    parser.add_argument("--flank", type=int, default=30000)
    parser.add_argument("--threads", type=int, default=1, help="Reserved for compatibility; feature extraction is currently single-process.")
    args = parser.parse_args()

    from tensorflow.keras.models import load_model

    pairs, metadata = read_gene_pairs(args.gene_pairs, args.bedpe_format)
    rnaseq, sample_cols = read_rnaseq_for_prediction(args.rnaseq_bed, args.column_names)
    if args.samples:
        missing = [s for s in args.samples if s not in sample_cols]
        if missing:
            raise ValueError(f"Requested samples not found in RNA-seq file: {missing}")
        sample_cols = args.samples

    model = load_model(args.model)
    output = metadata.copy()

    for sample in sample_cols:
        feature_df = extract_features(pairs, rnaseq[["chr", "start", "end", sample]], sample, args.bins, args.flank)
        output[f"AI4Loop_score_{sample}"] = predict_feature_df(feature_df, model)
        print(f"Predicted {len(feature_df)} pairs for sample {sample}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    print(f"Saved predictions to {args.output}")


if __name__ == "__main__":
    main()
