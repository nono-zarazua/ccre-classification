#!/bin/bash

MAIN="/work/zarazuanav/workspace/repos/ccre-classification"
PROCESSED="$MAIN/data/processed"
BED="$PROCESSED/hepg2_ELS.trimmed.bed"
FOLDS="$PROCESSED/ls-gkm/folds"


echo "Building fold 1 with chroms: 1, 20, 14, 22"
awk '$1=="chr1" || $1=="chr20" || $1=="chr14" || $1=="chr22"' \
  $BED > "$FOLDS/fold1/fold1.bed"

echo "Building fold 2 with chroms: 2, 12, 9, 13"
awk '$1=="chr2" || $1=="chr12" || $1=="chr9" || $1=="chr13"' \
  $BED > "$FOLDS/fold2/fold2.bed"

echo "Building fold 1 with chroms: 17, 11, 7, 15, X"
awk '$1=="chr17" || $1=="chr11" || $1=="chr7" || $1=="chr15" || $1=="chrX"' \
  $BED > "$FOLDS/fold3/fold3.bed"

echo "Building fold 1 with chroms: 6, 19, 5, 4, 21"
awk '$1=="chr6" || $1=="chr19" || $1=="chr5" || $1=="chr4" || $1=="chr21"' \
  $BED > "$FOLDS/fold4/fold4.bed"

echo "Building fold 1 with chroms: 3, 16, 10, 8, 18"
awk '$1=="chr3" || $1=="chr16" || $1=="chr10" || $1=="chr8" || $1=="chr18"' \
  $BED > "$FOLDS/fold5/fold5.bed"
