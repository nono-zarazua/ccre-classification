#!/bin/bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CHROMS_DIR="${REPO_ROOT}/data/processed/evn/GRCh38/chroms"
FOLDS_DIR="${REPO_ROOT}/data/processed/evn/GRCh38/folds"
TRAINING_DIR="${REPO_ROOT}/data/processed/evn/GRCh38/training"

# Fold assignments from data/processed/evn/GRCh38/folds/README.
fold_chroms_1=(1 7 16 X)
fold_chroms_2=(2 12 17 19)
fold_chroms_3=(3 8 9 18 22)
fold_chroms_4=(6 4 11 13 21)
fold_chroms_5=(5 10 15 14 20)

# AlphaGenome outer assignments from
# data/processed/evn/GRCh38/training/README.
# Outer N holds fold N out for evaluation. The listed chromosome is also
# removed from the remaining folds and used for early-stopping validation.
outer_folds_1=(2 3 4)
outer_folds_2=(1 3 4)
outer_folds_3=(1 2 4)
outer_folds_4=(1 2 3)

declare -A validation_chrom=(
    [1]=6
    [2]=6
    [3]=6
    [4]=3
)

source_file() {
    local chromosome="$1"
    local class_name="$2"
    local extension="$3"
    local chrom_dir="${CHROMS_DIR}/chrom${chromosome}"

    if [[ "${class_name}" == "pos" ]]; then
        printf '%s\n' "${chrom_dir}/chrom${chromosome}_GRCh38_ELS.${extension}"
    else
        printf '%s\n' "${chrom_dir}/neg1x_chrom${chromosome}_GRCh38_ELS.${extension}"
    fi
}

# Check every chromosome artifact before creating any merged output.
all_chromosomes=({1..22} X Y)
for chromosome in "${all_chromosomes[@]}"; do
    for class_name in pos neg; do
        for extension in bed fa; do
            input=$(source_file "${chromosome}" "${class_name}" "${extension}")
            if [[ ! -s "${input}" ]]; then
                echo "ERROR: required input is missing or empty: ${input}" >&2
                exit 1
            fi
        done
    done
done

merge_chromosomes() {
    local output="$1"
    local class_name="$2"
    local extension="$3"
    shift 3
    local chromosomes=("$@")
    local temporary="${output}.tmp.$$"

    mkdir -p "$(dirname "${output}")"
    : > "${temporary}"

    for chromosome in "${chromosomes[@]}"; do
        input=$(source_file "${chromosome}" "${class_name}" "${extension}")
        cat "${input}" >> "${temporary}"
    done

    mv "${temporary}" "${output}"
}

echo "Building folds 1-5"
for fold in {1..5}; do
    declare -n fold_chromosomes="fold_chroms_${fold}"
    fold_dir="${FOLDS_DIR}/fold${fold}"

    merge_chromosomes "${fold_dir}/fold${fold}.bed" pos bed "${fold_chromosomes[@]}"
    merge_chromosomes "${fold_dir}/fold${fold}.fa" pos fa "${fold_chromosomes[@]}"
    merge_chromosomes "${fold_dir}/neg1x_fold${fold}.bed" neg bed "${fold_chromosomes[@]}"
    merge_chromosomes "${fold_dir}/neg1x_fold${fold}.fa" neg fa "${fold_chromosomes[@]}"

    echo "Fold ${fold}: ${fold_chromosomes[*]}"
    unset -n fold_chromosomes
done

echo "Building AlphaGenome outer training and validation sets"
for outer in {1..4}; do
    declare -n included_folds="outer_folds_${outer}"
    validation="${validation_chrom[${outer}]}"
    training_chromosomes=()

    for fold in "${included_folds[@]}"; do
        declare -n fold_chromosomes="fold_chroms_${fold}"
        for chromosome in "${fold_chromosomes[@]}"; do
            if [[ "${chromosome}" != "${validation}" ]]; then
                training_chromosomes+=("${chromosome}")
            fi
        done
        unset -n fold_chromosomes
    done

    outer_dir="${TRAINING_DIR}/outer${outer}"

    merge_chromosomes "${outer_dir}/train_pos.bed" pos bed "${training_chromosomes[@]}"
    merge_chromosomes "${outer_dir}/train_pos.fa" pos fa "${training_chromosomes[@]}"
    merge_chromosomes "${outer_dir}/train_neg.bed" neg bed "${training_chromosomes[@]}"
    merge_chromosomes "${outer_dir}/train_neg.fa" neg fa "${training_chromosomes[@]}"

    merge_chromosomes "${outer_dir}/validation_pos.bed" pos bed "${validation}"
    merge_chromosomes "${outer_dir}/validation_pos.fa" pos fa "${validation}"
    merge_chromosomes "${outer_dir}/validation_neg.bed" neg bed "${validation}"
    merge_chromosomes "${outer_dir}/validation_neg.fa" neg fa "${validation}"

    echo "Outer ${outer}: train chromosomes ${training_chromosomes[*]}; validation chr${validation}; held-out fold ${outer}"
    unset -n included_folds
done

echo "Finished building chromosome folds and AlphaGenome outer sets"
