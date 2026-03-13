#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace_dir="${GWAS_WORKSPACE_DIR:-$(cd "${script_dir}/.." && pwd)}"
data_dir="${GWAS_DATA_DIR:-${workspace_dir}/data}"
wrapper="${script_dir}/run_with_top.sh"

id_file="${GWAS_ID_FILE:-${data_dir}/sample_ids.txt}"
lipid_file="${GWAS_LIPID_FILE:-${data_dir}/phenotype_matrix.tsv}"
gff_file="${GWAS_GFF_FILE:-${data_dir}/annotation.gff3}"

"${wrapper}" prepare_phenotypes \
  python3 "${script_dir}/prepare_phenotypes.py" \
  --id-file "${id_file}" \
  --lipid-file "${lipid_file}" \
  --outdir "${workspace_dir}/analysis/phenotypes"

"${wrapper}" prepare_annotation \
  python3 "${script_dir}/prepare_gene_annotation.py" \
  --gff "${gff_file}" \
  --outdir "${workspace_dir}/analysis/annotation"

"${wrapper}" prepare_genotypes \
  "${script_dir}/prepare_genotypes.sh"

"${wrapper}" prepare_emmax_inputs \
  python3 "${script_dir}/prepare_emmax_inputs.py" \
  --trait-matrix "${workspace_dir}/analysis/phenotypes/all_traits_int.tsv" \
  --trait-meta "${workspace_dir}/analysis/phenotypes/trait_metadata.tsv" \
  --covariates "${workspace_dir}/analysis/genotype/covariates_subset.tsv" \
  --outdir "${workspace_dir}/analysis/emmax_inputs"

"${wrapper}" emmax_batch \
  python3 "${script_dir}/run_emmax_batch.py" \
  --manifest "${workspace_dir}/analysis/emmax_inputs/trait_manifest.tsv" \
  --tped-prefix "${workspace_dir}/analysis/genotype/target_subset_tped" \
  --kinship "${workspace_dir}/analysis/genotype/kinship_subset.tsv" \
  --covariates "${workspace_dir}/analysis/emmax_inputs/covariates_emmax.tsv" \
  --parallel "${GWAS_EMMAX_PARALLEL:-10}" \
  --workspace "${workspace_dir}"

"${script_dir}/run_post_emmax_pipeline.sh"
