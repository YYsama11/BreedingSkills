#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace_dir="${GWAS_WORKSPACE_DIR:-$(cd "${script_dir}/.." && pwd)}"
wrapper="${script_dir}/run_with_top.sh"

data_dir="${GWAS_DATA_DIR:-${workspace_dir}/data}"
gwas_dir="${GWAS_GWAS_DIR:-${workspace_dir}/analysis/gwas/emmax}"
qtl_dir="${GWAS_QTL_DIR:-${workspace_dir}/analysis/qtl}"
annot_dir="${GWAS_ANNOT_DIR:-${workspace_dir}/analysis/annotation}"
geno_prefix="${GWAS_BFILE_PREFIX:-${data_dir}/genotype_panel}"
master_sig="${GWAS_MASTER_SIG:-${gwas_dir}/master_significant_snps.tsv.gz}"
trait_summaries="${GWAS_TRAIT_SUMMARIES:-${gwas_dir}/trait_summaries.tsv}"
annotation_gff="${GWAS_GFF_FILE:-}"
gene_annotation_file="${GWAS_GENE_ANNOTATION_FILE:-}"

mkdir -p "${qtl_dir}"

"${wrapper}" qtl_leads \
  python3 "${script_dir}/prepare_qtl_leads.py" \
  --master-significant "${master_sig}" \
  --outdir "${qtl_dir}" \
  --clump-kb "${GWAS_QTL_CLUMP_KB:-500}"

if [[ -s "${qtl_dir}/global_lead_snps.txt" ]]; then
  "${wrapper}" qtl_ld_batch \
    python3 "${script_dir}/run_ld_batch.py" \
    --lead-file "${qtl_dir}/global_lead_snps.txt" \
    --bfile-prefix "${geno_prefix}" \
    --workspace "${workspace_dir}" \
    --parallel "${GWAS_QTL_LD_PARALLEL:-2}"
fi

integrate_cmd=(python3 "${script_dir}/integrate_qtl_candidates.py" --lead-snps "${qtl_dir}/global_lead_snps.tsv" --ld-dir "${qtl_dir}/ld" --master-significant "${master_sig}" --outdir "${qtl_dir}")

if [[ -f "${trait_summaries}" ]]; then
  integrate_cmd+=(--trait-summaries "${trait_summaries}")
fi

if [[ -n "${annotation_gff}" && -f "${annotation_gff}" ]]; then
  "${wrapper}" prepare_annotation \
    python3 "${script_dir}/prepare_gene_annotation.py" \
    --gff "${annotation_gff}" \
    --outdir "${annot_dir}"
  integrate_cmd+=(--annotation-dir "${annot_dir}")
fi

if [[ -n "${gene_annotation_file}" && -f "${gene_annotation_file}" ]]; then
  integrate_cmd+=(--gene-annotation-file "${gene_annotation_file}")
fi

"${wrapper}" qtl_integrate "${integrate_cmd[@]}"
