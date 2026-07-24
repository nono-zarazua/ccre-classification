#!/bin/bash

set -euo pipefail

TEST_FOLD="${1:?Usage: make_outer_files.sh <fold_number>}"

DATA_DIR="data/processed/folds"
TRAIN_DIR="data/processed/ls-gkm/training/outer${TEST_FOLD}"

mkdir -p "${TRAIN_DIR}"

TRAIN_POS="${TRAIN_DIR}/train_pos.fa"
TRAIN_NEG="${TRAIN_DIR}/train_neg.fa"

: > "${TRAIN_POS}"
: > "${TRAIN_NEG}"

case "${TEST_FOLD}" in
    1) HELDOUT_REGEX='chr1|chr20|chr14|chr22' ;;
    2) HELDOUT_REGEX='chr2|chr12|chr9|chr13' ;;
    3) HELDOUT_REGEX='chr17|chr11|chr7|chr15|chrX' ;;
    4) HELDOUT_REGEX='chr6|chr19|chr5|chr4|chr21' ;;
    5) HELDOUT_REGEX='chr3|chr16|chr10|chr8|chr18' ;;
    *) echo "ERROR: fold must be 1–5" >&2; exit 1 ;;
esac

for FOLD in 1 2 3 4 5; do
    if [[ "${FOLD}" -ne "${TEST_FOLD}" ]]; then
        cat "${DATA_DIR}/fold${FOLD}/fold${FOLD}.fa" >> "${TRAIN_POS}"
        cat "${DATA_DIR}/fold${FOLD}/neg1x_fold${FOLD}.fa" >> "${TRAIN_NEG}"
    fi
done

if grep '^>' "${TRAIN_POS}" | grep -Eq "^>(${HELDOUT_REGEX})_"; then
    echo "ERROR: held-out chromosomes found in ${TRAIN_POS}" >&2
    exit 1
fi

echo "Check passed: held-out chromosomes absent"
echo "Training positive sequences: $(grep -c '^>' "${TRAIN_POS}")"
echo "Training negative sequences: $(grep -c '^>' "${TRAIN_NEG}")"
