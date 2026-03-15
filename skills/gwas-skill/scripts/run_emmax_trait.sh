#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <trait_id> <pheno_file> <tped_prefix> <kinship_file> <cov_file>" >&2
  exit 1
fi

trait_id="$1"
pheno_file="$2"
tped_prefix="$3"
kinship_file="$4"
cov_file="$5"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace_dir="${GWAS_WORKSPACE_DIR:-$(cd "${script_dir}/.." && pwd)}"
outdir="${GWAS_EMMAX_RESULTS_DIR:-${workspace_dir}/analysis/gwas/emmax/results}"
mkdir -p "${outdir}"
emmax_bin="${EMMAX_BIN:-$(command -v emmax-intel64 2>/dev/null || true)}"
keep_tool_logs="$(printf '%s' "${KEEP_TOOL_LOGS:-false}" | tr '[:upper:]' '[:lower:]')"
if [[ -z "${emmax_bin}" ]]; then
  echo "EMMAX_BIN is not set and emmax-intel64 was not found in PATH." >&2
  exit 1
fi

if [[ -f "${outdir}/${trait_id}.ps.gz" && -f "${outdir}/${trait_id}.reml" ]]; then
  exit 0
fi

if [[ -f "${outdir}/${trait_id}.ps" && -f "${outdir}/${trait_id}.reml" ]]; then
  exit 0
fi

"${emmax_bin}" \
  -t "${tped_prefix}" \
  -p "${pheno_file}" \
  -k "${kinship_file}" \
  -c "${cov_file}" \
  -d 12 \
  -o "${outdir}/${trait_id}"

if [[ "${keep_tool_logs}" != "true" ]]; then
  rm -f "${outdir}/${trait_id}.log"
fi
