#!/bin/bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CHROMS_DIR="${REPO_ROOT}/data/processed/evn/GRCh38/chroms"
FOLDS_DIR="${REPO_ROOT}/data/processed/evn/GRCh38/folds"
TRAINING_DIR="${REPO_ROOT}/data/processed/evn/GRCh38/training"
MASTER_SCRIPT="${REPO_ROOT}/src/make_GRCh38_master_training.py"
MASTER_TARGET=202493
MASTER_SEED=42

echo "Sampling one chromosome-proportional master dataset and deriving all folds/outers"
python "${MASTER_SCRIPT}" \
    --chroms-dir "${CHROMS_DIR}" \
    --folds-dir "${FOLDS_DIR}" \
    --output-dir "${TRAINING_DIR}" \
    --target-count "${MASTER_TARGET}" \
    --seed "${MASTER_SEED}" \
    --overwrite

echo "Finished building the sampled master dataset, folds, and AlphaGenome outer sets"
