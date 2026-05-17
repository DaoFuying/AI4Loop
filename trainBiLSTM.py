#!/usr/bin/env python3
"""Backward-compatible entry point for AI4Loop training.

Preferred usage:
    python train.py --features FEATURES.csv --output_model MODEL.h5 --metrics_out METRICS.json

Legacy usage is still supported:
    python trainBiLSTM.py k562_ctcf out_dir
"""

import os
import sys
from pathlib import Path

from train import main as train_main


def legacy_main():
    if len(sys.argv) == 3 and not sys.argv[1].startswith("--"):
        eachcell = sys.argv[1]
        out_dir = sys.argv[2]
        feature_file = Path(out_dir) / f"{eachcell}_distance_matched.csv_winGEfea.csv"
        model_file = Path(out_dir) / f"{eachcell}_RandomSplit.model.h5"
        metrics_file = Path(out_dir) / f"{eachcell}_RandomSplit.metrics.json"
        sys.argv = [
            sys.argv[0],
            "--features", str(feature_file),
            "--output_model", str(model_file),
            "--metrics_out", str(metrics_file),
            "--split", "random",
        ]
    train_main()


if __name__ == "__main__":
    legacy_main()
