# Supervisor handoff: full cleaned-EVN LS-GKM benchmark

This portable package trains one LS-GKM model on the largest balanced cleaned-EVN dataset supported by the available real-genomic backgrounds. It is separate from the matched-size EVN/SVN chromosome-learning benchmark.

## Data design

EVN was deduplicated, all loci overlapping SVN were removed, and chrY was excluded. Positive and real-negative membership was balanced per chromosome, then exact negative coordinates were deduplicated globally across splits.

| Split | Chromosomes | Records per class |
|---|---|---:|
| Training | chr3–22 and chrX | 1,105,087 |
| Validation | chr1 | 123,034 |
| Test | chr2 | 128,663 |
| **Total** | chr1–22/X | **1,356,784** |

All splits are balanced 1:1. Validation and test are not used to fit the model.

## Contents

- `positive/` and `negative/`: BED/FASTA split files.
- `scripts/train_full_lsgkm.sh`: training plus training/validation/test prediction.
- `ls-gkm.yml`: Conda environment.
- `manifest.tsv`, `summary.json`, `PROVENANCE.md`, and `checksums.sha256`: provenance and integrity metadata.

LS-GKM itself is not included.

## Verify, install, and submit

```bash
sha256sum -c full_evn_benchmark.tar.zst.sha256
tar --zstd -xf full_evn_benchmark.tar.zst
cd full-benchmark
sha256sum -c checksums.sha256
conda env create -f ls-gkm.yml
conda activate ls-gkm
git clone https://github.com/Dongwon-Lee/lsgkm.git
cd lsgkm && git checkout d7ef90f7cb326b190619c51a13843da48c99e239
cd src && make && make install && cd ../..
test -x lsgkm/bin/gkmtrain
test -x lsgkm/bin/gkmpredict
sbatch scripts/train_full_lsgkm.sh "$PWD/lsgkm/bin"
```

Review the Slurm header first. The supplied request is 16 CPUs, 500 GB RAM, and 14 days; `gkmtrain` uses a 400,000 MB cache. Outputs are written under `outputs/lsgkm/`: `full_evn.model.txt` and `{100_train,validation,test}_{pos,neg}_scores.txt`. Retain these and the Slurm logs, including job ID, elapsed time, and maximum resident memory.
