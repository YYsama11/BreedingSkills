#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <log-prefix> <command...>" >&2
  exit 1
fi

log_prefix="$1"
shift

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace_dir="${GWAS_WORKSPACE_DIR:-$(cd "${script_dir}/.." && pwd)}"
log_dir="${workspace_dir}/analysis/logs/runtime"
mkdir -p "${log_dir}"

timestamp="$(date '+%Y%m%d_%H%M%S')"
top_log="${log_dir}/${timestamp}_${log_prefix}.top.log"
cmd_log="${log_dir}/${timestamp}_${log_prefix}.cmd.log"
manifest="${log_dir}/command_manifest.tsv"

top -bn1 | head -n 8 | tee "${top_log}"

idle="$(awk -F',' '/%Cpu/{for(i=1;i<=NF;i++){if($i ~ / id/){gsub(/[^0-9.]/,\"\",$i); print $i; exit}}}' "${top_log}")"
if [[ -z "${idle}" ]]; then
  idle="0"
fi

if awk "BEGIN{exit !(${idle} < 20)}"; then
  threads=1
elif awk "BEGIN{exit !(${idle} < 40)}"; then
  threads=2
elif awk "BEGIN{exit !(${idle} < 60)}"; then
  threads=4
elif awk "BEGIN{exit !(${idle} < 75)}"; then
  threads=6
else
  threads=8
fi

if [[ "${threads}" -gt 30 ]]; then
  threads=30
fi

export OMP_NUM_THREADS="${threads}"
export OPENBLAS_NUM_THREADS="${threads}"
export MKL_NUM_THREADS="${threads}"
export NUMEXPR_MAX_THREADS="${threads}"
export VECLIB_MAXIMUM_THREADS="${threads}"
export GWAS_THREADS="${threads}"

if [[ ! -f "${manifest}" ]]; then
  printf "timestamp\tlog_prefix\tcpu_idle\trecommended_threads\tcommand\n" > "${manifest}"
fi

printf "%s\t%s\t%s\t%s\t%s\n" \
  "${timestamp}" "${log_prefix}" "${idle}" "${threads}" "$*" >> "${manifest}"

{
  echo "[run_with_top] timestamp=${timestamp}"
  echo "[run_with_top] cpu_idle=${idle}"
  echo "[run_with_top] recommended_threads=${threads}"
  echo "[run_with_top] command=$*"
  "$@"
} 2>&1 | tee "${cmd_log}"
