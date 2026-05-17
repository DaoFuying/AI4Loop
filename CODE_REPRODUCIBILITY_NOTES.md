# Code reproducibility notes

This repository revision was reorganized to improve transparency and reproducibility for peer review.

## Main changes

1. Replaced hard-coded prediction and training paths with command-line arguments.
2. Added `train.py` with explicit random/chromosome split options, validation set, random seed, metrics output, ROC/PR curve export and optional class weighting.
3. Added `evaluate.py` for independent evaluation of pretrained models.
4. Added `predict.py` for command-line prediction from RNA-seq and gene-pair files.
5. Added `extract_rnaseq_features.py` for multi-scale RNA-seq feature extraction without editing source code.
6. Kept `AI4LoopPrediction.py`, `trainBiLSTM.py` and `callRNAseq.py` as backward-compatible wrappers.
7. Added `requirements.txt`, `environment.yml`, `.gitignore` and a new README.

## Split definitions

- Random split: stratified 80/20 train/test split with an internal validation split.
- Chromosome split: chromosomes 4, 7, 8 and 11 are held out for testing by default.

## Metrics

AUROC and AUPRC are calculated from continuous prediction probabilities. Accuracy, F1 and MCC use a default threshold of 0.5 unless specified otherwise.
