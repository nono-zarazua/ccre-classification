#!/usr/bin/env bash
## Run one GkmExplain mode over a FASTA using independent sequence shards.
set -euo pipefail

if (( $# != 6 )); then
  echo "Usage: $0 FASTA MODEL OUTPUT MODE N_SHARDS GKMEXPLAIN_BIN" >&2
  exit 2
fi

input_fasta="$1"
model="$2"
output="$3"
mode="$4"
n_shards="$5"
gkmexplain_bin="$6"

[[ -s "${input_fasta}" ]] || { echo "ERROR: missing/empty ${input_fasta}" >&2; exit 1; }
[[ -s "${model}" ]] || { echo "ERROR: missing/empty ${model}" >&2; exit 1; }
[[ -x "${gkmexplain_bin}" ]] || { echo "ERROR: not executable: ${gkmexplain_bin}" >&2; exit 1; }
[[ "${mode}" =~ ^[0-5]$ ]] || { echo "ERROR: MODE must be 0-5" >&2; exit 1; }
[[ "${n_shards}" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: N_SHARDS must be positive" >&2; exit 1; }
[[ ! -s "${output}" ]] || { echo "SKIP: ${output}"; exit 0; }

work_dir="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/gkmexplain.${SLURM_JOB_ID:-local}.XXXXXX")"
trap 'rm -rf -- "${work_dir}"' EXIT

norc_lines="$(grep -c '^norc ' "${model}" || true)"
if (( norc_lines > 0 )); then
  (( norc_lines == 1 )) || { echo "ERROR: multiple norc headers in ${model}" >&2; exit 1; }
  [[ "$(grep '^norc ' "${model}")" == "norc 0" ]] || {
    echo "ERROR: only a non-operative 'norc 0' header can be converted" >&2
    exit 1
  }
  sed '/^norc 0$/d' "${model}" > "${work_dir}/compatible.model.txt"
else
  cp "${model}" "${work_dir}/compatible.model.txt"
fi

awk -v out="${work_dir}" -v n="${n_shards}" '
  /^>/ { shard=(records++ % n); file=sprintf("%s/input.%03d.fa", out, shard) }
  { print >> file }
' "${input_fasta}"

pids=()
for ((shard=0; shard<n_shards; shard++)); do
  shard_fasta="$(printf '%s/input.%03d.fa' "${work_dir}" "${shard}")"
  shard_output="$(printf '%s/output.%03d.txt' "${work_dir}" "${shard}")"
  [[ -s "${shard_fasta}" ]] || continue
  "${gkmexplain_bin}" -v 0 -m "${mode}" \
    "${shard_fasta}" "${work_dir}/compatible.model.txt" "${shard_output}" &
  pids+=("$!")
done

failed=0
for pid in "${pids[@]}"; do
  wait "${pid}" || failed=1
done
(( failed == 0 )) || { echo "ERROR: one or more GkmExplain shards failed" >&2; exit 1; }

temporary_output="${work_dir}/aggregated.txt"
: > "${temporary_output}"
for ((shard=0; shard<n_shards; shard++)); do
  shard_output="$(printf '%s/output.%03d.txt' "${work_dir}" "${shard}")"
  [[ ! -e "${shard_output}" ]] || cat "${shard_output}" >> "${temporary_output}"
done

expected="$(grep -c '^>' "${input_fasta}")"
observed="$(wc -l < "${temporary_output}")"
(( expected > 0 && observed == expected )) || {
  echo "ERROR: expected ${expected} rows, observed ${observed}" >&2
  exit 1
}
mkdir -p "$(dirname "${output}")"
mv "${temporary_output}" "${output}"
echo "PASS: mode=${mode}; sequences=${observed}; output=${output}"
