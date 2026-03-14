#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

trait_label="${QTL_PLOT_TRAIT_LABEL:-QTL summary}"
chrom_sizes="${QTL_PLOT_CHROM_SIZES:-}"
global_manhattan="${QTL_PLOT_GLOBAL_MANHATTAN:-}"
qtl_regions="${QTL_PLOT_QTL_REGIONS:-}"
local_manhattan="${QTL_PLOT_LOCAL_MANHATTAN:-}"
gene_models="${QTL_PLOT_GENE_MODELS:-}"
highlight_genes="${QTL_PLOT_HIGHLIGHT_GENES:-}"
out_path="${QTL_PLOT_OUT:-qtl_summary_plot.pdf}"
y_min="${QTL_PLOT_Y_MIN:-3}"
max_loci="${QTL_PLOT_MAX_LOCI:-0}"

cmd=(Rscript "${script_dir}/plot_ld_qtl_summary.R" --trait-label "${trait_label}" --global-manhattan "${global_manhattan}" --qtl-regions "${qtl_regions}" --local-manhattan "${local_manhattan}" --out "${out_path}" --y-min "${y_min}" --max-loci "${max_loci}")

if [[ -n "${chrom_sizes}" ]]; then
  cmd+=(--chrom-sizes "${chrom_sizes}")
fi
if [[ -n "${gene_models}" ]]; then
  cmd+=(--gene-models "${gene_models}")
fi
if [[ -n "${highlight_genes}" ]]; then
  cmd+=(--highlight-genes "${highlight_genes}")
fi

"${cmd[@]}"
