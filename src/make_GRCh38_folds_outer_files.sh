#!/bin/bash

set -euo pipefail

if [[ "$#" -ne 1 ]]; then
    echo "Usage: $0 <evn|svn>" >&2
    exit 1
fi

dataset="$1"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MASTER_SEED=42

case "${dataset}" in
    evn)
        MASTER_SCRIPT="${REPO_ROOT}/src/make_GRCh38_master_training_evn.py"
        CHROMS_DIR="${REPO_ROOT}/data/processed/evn/GRCh38/chroms"
        FOLDS_DIR="${REPO_ROOT}/data/processed/evn/GRCh38/folds"
        TRAINING_DIR="${REPO_ROOT}/data/processed/evn/GRCh38/training"
        echo "Sampling the chromosome-proportional EVN master dataset and deriving folds/outers"
        python "${MASTER_SCRIPT}" \
            --chroms-dir "${CHROMS_DIR}" \
            --folds-dir "${FOLDS_DIR}" \
            --output-dir "${TRAINING_DIR}" \
            --target-count 202493 \
            --seed "${MASTER_SEED}" \
            --overwrite
        ;;
    svn)
        MASTER_SCRIPT="${REPO_ROOT}/src/make_GRCh38_master_training_svn.py"
        CHROMS_DIR="${REPO_ROOT}/data/processed/svn/GRCh38/chroms"
        FOLDS_DIR="${REPO_ROOT}/data/processed/svn/GRCh38/folds"
        TRAINING_DIR="${REPO_ROOT}/data/processed/svn/GRCh38/training"
        echo "Building the native-distribution SVN master dataset and deriving folds/outers"
        python "${MASTER_SCRIPT}" \
            --chroms-dir "${CHROMS_DIR}" \
            --folds-dir "${FOLDS_DIR}" \
            --output-dir "${TRAINING_DIR}" \
            --seed "${MASTER_SEED}" \
            --overwrite
        ;;
    *)
        echo "ERROR: dataset must be 'evn' or 'svn'; received '${dataset}'" >&2
        exit 1
        ;;
esac

printf 'Finished building %s GRCh38 master, folds, and outer sets\n' "${dataset}"
