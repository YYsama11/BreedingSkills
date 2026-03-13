#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <bfile_prefix> <sig_snp_file> <group_dir> <outdir>" >&2
  exit 1
fi

bfile_prefix="$1"
sig_snp_file="$2"
group_dir="$3"
outdir="$4"
mkdir -p "${outdir}"
plink_bin="${PLINK_BIN:-$(command -v plink 2>/dev/null || true)}"
if [[ -z "${plink_bin}" ]]; then
  echo "PLINK_BIN is not set and plink was not found in PATH." >&2
  exit 1
fi

for keep_file in "${group_dir}"/Group*.keep.txt; do
  if [[ ! -f "${keep_file}" ]]; then
    continue
  fi
  group_name="$(basename "${keep_file}" .keep.txt)"
  "${plink_bin}" \
    --bfile "${bfile_prefix}" \
    --extract "${sig_snp_file}" \
    --keep "${keep_file}" \
    --freq \
    --allow-no-sex \
    --threads "${GWAS_THREADS:-2}" \
    --out "${outdir}/${group_name}"
done
