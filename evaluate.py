#!/usr/bin/env python3
"""Evaluate a trained AI4Loop model on an extracted feature matrix."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)


def load_feature_matrix(path: str, label_col: str = "label"):
    data = pd.read_csv(path)
    data.columns = data.columns.astype(str).str.strip()
    if label_col not in data.columns:
        raise ValueError(f"Label column '{label_col}' was not found in {path}.")
    left_cols = [c for c in data.columns if c.startswith("L")]
    right_cols = [c for c in data.columns if c.startswith("R")]
    if not left_cols or not right_cols:
        raise ValueError("Feature matrix must contain columns beginning with L and R.")
    left = data[left_cols].replace([np.inf, -np.inf], np.nan).fillna(0).astype("float32").values.reshape(-1, len(left_cols), 1)
    right = data[right_cols].replace([np.inf, -np.inf], np.nan).fillna(0).astype("float32").values.reshape(-1, len(right_cols), 1)
    y = data[label_col].astype(int).values
    return data, left, right, y, left_cols, right_cols


def calculate_metrics(y_true, prob, threshold: float):
    pred = (prob >= threshold).astype(int)
    return {
        "n_samples": int(len(y_true)),
        "positive_fraction": float(np.mean(y_true)),
        "threshold": float(threshold),
        "auroc": float(roc_auc_score(y_true, prob)),
        "auprc": float(average_precision_score(y_true, prob)),
        "accuracy": float(accuracy_score(y_true, pred)),
        "f1": float(f1_score(y_true, pred)),
        "mcc": float(matthews_corrcoef(y_true, pred)),
    }


def save_curves(y_true, prob, prefix: Path):
    fpr, tpr, _ = roc_curve(y_true, prob)
    precision, recall, _ = precision_recall_curve(y_true, prob)
    pd.DataFrame({"FPR": fpr, "TPR": tpr}).to_csv(f"{prefix}.roc.csv", index=False)
    pd.DataFrame({"Precision": precision, "Recall": recall}).to_csv(f"{prefix}.pr.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description="Evaluate an AI4Loop .h5 model.")
    parser.add_argument("--model", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--metrics_out", required=True)
    parser.add_argument("--predictions_out", default=None)
    parser.add_argument("--label_col", default="label")
    parser.add_argument("--threshold", type=float, default=0.5)
    args = parser.parse_args()

    from tensorflow.keras.models import load_model

    data, left, right, y, left_cols, right_cols = load_feature_matrix(args.features, args.label_col)
    model = load_model(args.model)
    prob = model.predict([left, right], verbose=0).ravel()

    metrics = calculate_metrics(y, prob, args.threshold)
    metrics.update({
        "model_file": args.model,
        "feature_file": args.features,
        "left_feature_count": len(left_cols),
        "right_feature_count": len(right_cols),
    })

    out = Path(args.metrics_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as handle:
        json.dump(metrics, handle, indent=2)
    save_curves(y, prob, out.with_suffix(""))

    if args.predictions_out:
        pred_df = data.copy()
        pred_df["AI4Loop_score"] = prob
        pred_df.to_csv(args.predictions_out, index=False)

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
