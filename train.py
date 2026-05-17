#!/usr/bin/env python3
"""Train AI4Loop models with reproducible random or chromosome splits."""

import argparse
import json
import os
import random
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
from sklearn.model_selection import train_test_split
from sklearn.utils.class_weight import compute_class_weight


def set_seed(seed: int) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    try:
        import tensorflow as tf
        tf.random.set_seed(seed)
    except Exception:
        pass


def normalize_chr(value) -> str:
    value = str(value)
    if value.startswith("chr"):
        return value
    if value in {"23", "X"}:
        return "chrX"
    if value in {"24", "Y"}:
        return "chrY"
    return f"chr{value}"


def build_model(rna_feature_length: int, lstm_units=None, dense_units=None, dropout: float = 0.3):
    from tensorflow.keras.layers import Bidirectional, Concatenate, Dense, Dropout, Input, LSTM
    from tensorflow.keras.models import Model

    if lstm_units is None:
        lstm_units = [16, 16]
    if dense_units is None:
        dense_units = [64, 32]

    input_left = Input(shape=(rna_feature_length, 1), dtype="float32", name="left_anchor")
    input_right = Input(shape=(rna_feature_length, 1), dtype="float32", name="right_anchor")

    x_left = input_left
    x_right = input_right

    for units in lstm_units[:-1]:
        layer = Bidirectional(LSTM(units, return_sequences=True))
        x_left = layer(x_left)
        x_right = layer(x_right)

    layer = Bidirectional(LSTM(lstm_units[-1], return_sequences=False))
    x_left = layer(x_left)
    x_right = layer(x_right)

    x = Concatenate()([x_left, x_right])
    for units in dense_units:
        x = Dense(units, activation="relu")(x)
        if dropout and dropout > 0:
            x = Dropout(dropout)(x)

    output = Dense(1, activation="sigmoid", name="interaction_probability")(x)
    return Model(inputs=[input_left, input_right], outputs=output)


def load_feature_matrix(path: str, label_col: str = "label"):
    data = pd.read_csv(path)
    data.columns = data.columns.astype(str).str.strip()
    if label_col not in data.columns:
        raise ValueError(f"Label column '{label_col}' was not found in {path}.")

    left_cols = [c for c in data.columns if c.startswith("L")]
    right_cols = [c for c in data.columns if c.startswith("R")]
    if not left_cols or not right_cols:
        raise ValueError("Feature matrix must contain left-anchor columns beginning with 'L' and right-anchor columns beginning with 'R'.")

    left = data[left_cols].replace([np.inf, -np.inf], np.nan).fillna(0).astype("float32").values.reshape(-1, len(left_cols), 1)
    right = data[right_cols].replace([np.inf, -np.inf], np.nan).fillna(0).astype("float32").values.reshape(-1, len(right_cols), 1)
    y = data[label_col].astype(int).values
    return data, left, right, y, left_cols, right_cols


def make_split(data, y, split: str, test_size: float, val_size: float, seed: int, test_chromosomes):
    indices = np.arange(len(y))
    if split == "random":
        train_val_idx, test_idx = train_test_split(indices, test_size=test_size, random_state=seed, stratify=y)
        y_train_val = y[train_val_idx]
        val_fraction = val_size / (1.0 - test_size)
        train_idx, val_idx = train_test_split(train_val_idx, test_size=val_fraction, random_state=seed, stratify=y_train_val)
        return train_idx, val_idx, test_idx

    required = {"chr1", "chr2"}
    if not required.issubset(set(data.columns)):
        raise ValueError("Chromosome split requires columns 'chr1' and 'chr2' in the feature matrix.")
    test_set = {normalize_chr(c) for c in test_chromosomes}
    chrom1 = data["chr1"].map(normalize_chr)
    chrom2 = data["chr2"].map(normalize_chr)
    test_mask = chrom1.isin(test_set) | chrom2.isin(test_set)
    test_idx = indices[test_mask.values]
    train_val_idx = indices[~test_mask.values]
    if len(test_idx) == 0:
        raise ValueError("Chromosome split produced an empty test set. Check chromosome names and --test_chromosomes.")
    train_idx, val_idx = train_test_split(train_val_idx, test_size=val_size, random_state=seed, stratify=y[train_val_idx])
    return train_idx, val_idx, test_idx


def calculate_metrics(y_true, prob, threshold: float = 0.5):
    pred = (prob >= threshold).astype(int)
    metrics = {
        "n_samples": int(len(y_true)),
        "positive_fraction": float(np.mean(y_true)),
        "threshold": float(threshold),
        "auroc": float(roc_auc_score(y_true, prob)),
        "auprc": float(average_precision_score(y_true, prob)),
        "accuracy": float(accuracy_score(y_true, pred)),
        "f1": float(f1_score(y_true, pred)),
        "mcc": float(matthews_corrcoef(y_true, pred)),
    }
    return metrics


def save_curves(y_true, prob, prefix: Path):
    fpr, tpr, _ = roc_curve(y_true, prob)
    precision, recall, _ = precision_recall_curve(y_true, prob)
    pd.DataFrame({"FPR": fpr, "TPR": tpr}).to_csv(f"{prefix}.roc.csv", index=False)
    pd.DataFrame({"Precision": precision, "Recall": recall}).to_csv(f"{prefix}.pr.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description="Train an AI4Loop Bi-LSTM model from extracted RNA-seq features.")
    parser.add_argument("--features", required=True, help="Feature CSV generated by extract_rnaseq_features.py or callRNAseq.py")
    parser.add_argument("--output_model", required=True, help="Output .h5 model path")
    parser.add_argument("--metrics_out", required=True, help="Output JSON metrics path")
    parser.add_argument("--label_col", default="label")
    parser.add_argument("--split", choices=["random", "chromosome"], default="random")
    parser.add_argument("--test_chromosomes", nargs="+", default=["chr4", "chr7", "chr8", "chr11"])
    parser.add_argument("--test_size", type=float, default=0.2)
    parser.add_argument("--val_size", type=float, default=0.1)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch_size", type=int, default=30)
    parser.add_argument("--learning_rate", type=float, default=0.001)
    parser.add_argument("--dropout", type=float, default=0.3)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--use_class_weight", action="store_true", help="Use class weights during model fitting")
    args = parser.parse_args()

    set_seed(args.seed)
    from tensorflow.keras.optimizers import Adam

    data, left, right, y, left_cols, right_cols = load_feature_matrix(args.features, args.label_col)
    train_idx, val_idx, test_idx = make_split(data, y, args.split, args.test_size, args.val_size, args.seed, args.test_chromosomes)

    model = build_model(left.shape[1], dropout=args.dropout)
    model.compile(optimizer=Adam(learning_rate=args.learning_rate), loss="binary_crossentropy", metrics=["accuracy"])

    class_weight = None
    if args.use_class_weight:
        classes = np.unique(y[train_idx])
        weights = compute_class_weight(class_weight="balanced", classes=classes, y=y[train_idx])
        class_weight = {int(c): float(w) for c, w in zip(classes, weights)}

    history = model.fit(
        [left[train_idx], right[train_idx]],
        y[train_idx],
        validation_data=([left[val_idx], right[val_idx]], y[val_idx]),
        epochs=args.epochs,
        batch_size=args.batch_size,
        class_weight=class_weight,
        verbose=2,
    )

    Path(args.output_model).parent.mkdir(parents=True, exist_ok=True)
    model.save(args.output_model)

    test_prob = model.predict([left[test_idx], right[test_idx]], verbose=0).ravel()
    val_prob = model.predict([left[val_idx], right[val_idx]], verbose=0).ravel()

    metrics = {
        "feature_file": args.features,
        "model_file": args.output_model,
        "split": args.split,
        "seed": args.seed,
        "left_feature_count": len(left_cols),
        "right_feature_count": len(right_cols),
        "n_train": int(len(train_idx)),
        "n_validation": int(len(val_idx)),
        "n_test": int(len(test_idx)),
        "validation": calculate_metrics(y[val_idx], val_prob, args.threshold),
        "test": calculate_metrics(y[test_idx], test_prob, args.threshold),
        "class_weight": class_weight,
    }

    metrics_path = Path(args.metrics_out)
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    with open(metrics_path, "w") as handle:
        json.dump(metrics, handle, indent=2)

    prefix = metrics_path.with_suffix("")
    save_curves(y[test_idx], test_prob, prefix)
    pd.DataFrame(history.history).to_csv(f"{prefix}.history.csv", index=False)
    print(json.dumps(metrics["test"], indent=2))


if __name__ == "__main__":
    main()
