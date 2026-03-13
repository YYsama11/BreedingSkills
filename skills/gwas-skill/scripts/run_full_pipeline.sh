#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace_dir="${GWAS_WORKSPACE_DIR:-$(cd "${script_dir}/.." && pwd)}"
data_dir="${GWAS_DATA_DIR:-${workspace_dir}/data}"
wrapper="${script_dir}/run_with_top.sh"

trait_matrix="${GWAS_TRAIT_MATRIX:-${data_dir}/trait_matrix.tsv}"
trait_meta="${GWAS_TRAIT_META:-}"
covariates="${GWAS_COVARIATES:-${data_dir}/covariates.tsv}"
tped_prefix="${GWAS_TPED_PREFIX:-${data_dir}/genotype_panel}"
kinship="${GWAS_KINSHIP:-${data_dir}/kinship_matrix.tsv}"
reference_fai="${GWAS_FASTA_FAI_FILE:-}"

prepare_cmd=(python3 "${script_dir}/prepare_emmax_inputs.py" --trait-matrix "${trait_matrix}" --covariates "${covariates}" --outdir "${workspace_dir}/analysis/emmax_inputs")
if [[ -n "${trait_meta}" && -f "${trait_meta}" ]]; then
  prepare_cmd+=(--trait-meta "${trait_meta}")
fi

"${wrapper}" prepare_emmax_inputs "${prepare_cmd[@]}"

"${wrapper}" emmax_batch \
  python3 "${script_dir}/run_emmax_batch.py" \
  --manifest "${workspace_dir}/analysis/emmax_inputs/trait_manifest.tsv" \
  --tped-prefix "${tped_prefix}" \
  --kinship "${kinship}" \
  --covariates "${workspace_dir}/analysis/emmax_inputs/covariates_emmax.tsv" \
  --parallel "${GWAS_EMMAX_PARALLEL:-10}" \
  --workspace "${workspace_dir}"

ref_cmd=(python3 "${script_dir}/prepare_snp_reference.py" --tped "${tped_prefix}.tped" --outdir "${workspace_dir}/analysis/gwas/reference")
if [[ -n "${reference_fai}" && -f "${reference_fai}" ]]; then
  ref_cmd+=(--fai "${reference_fai}")
fi
"${wrapper}" prepare_snp_reference "${ref_cmd[@]}"

"${wrapper}" summarize_emmax \
  python3 "${script_dir}/summarize_emmax_results.py" \
  --manifest "${workspace_dir}/analysis/emmax_inputs/trait_manifest.tsv" \
  --results-dir "${workspace_dir}/analysis/gwas/emmax/results" \
  --reference-dir "${workspace_dir}/analysis/gwas/reference" \
  --trait-meta "${trait_meta:-${workspace_dir}/analysis/emmax_inputs/trait_manifest.tsv}" \
  --outdir "${workspace_dir}/analysis/gwas/emmax" \
  --workers "${GWAS_THREADS:-2}"
