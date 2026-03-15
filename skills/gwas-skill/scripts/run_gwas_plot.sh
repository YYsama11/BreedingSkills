#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

trait="${GWAS_PLOT_TRAIT:-GWAS trait}"
gwas_table="${GWAS_PLOT_TABLE:-}"
chrom_sizes="${GWAS_PLOT_CHROM_SIZES:-}"
out_prefix="${GWAS_PLOT_OUT_PREFIX:-gwas_plot}"

cmd=(Rscript "${script_dir}/plot_manhattan_qq.R" --trait "${trait}" --gwas-table "${gwas_table}" --out-prefix "${out_prefix}")
if [[ -n "${chrom_sizes}" ]]; then
  cmd+=(--chrom-sizes "${chrom_sizes}")
fi
"${cmd[@]}"
