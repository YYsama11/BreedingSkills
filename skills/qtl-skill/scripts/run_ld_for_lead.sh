#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <lead_snp> <bfile_prefix> <out_prefix>" >&2
  exit 1
fi

lead_snp="$1"
bfile_prefix="$2"
out_prefix="$3"
plink_bin="${PLINK_BIN:-$(command -v plink 2>/dev/null || true)}"
if [[ -z "${plink_bin}" ]]; then
  echo "PLINK_BIN is not set and plink was not found in PATH." >&2
  exit 1
fi

if [[ -f "${out_prefix}.ld" ]]; then
  exit 0
fi

"${plink_bin}" \
  --bfile "${bfile_prefix}" \
  --ld-snp "${lead_snp}" \
  --r2 \
  --ld-window-kb 500 \
  --ld-window 99999 \
  --ld-window-r2 0.2 \
  --allow-no-sex \
  --threads "${GWAS_THREADS:-2}" \
  --out "${out_prefix}"
