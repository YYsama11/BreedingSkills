#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace_dir="${GWAS_WORKSPACE_DIR:-$(cd "${script_dir}/.." && pwd)}"
wrapper="${script_dir}/run_with_top.sh"

results_dir="${workspace_dir}/analysis/gwas/emmax/results"
gwas_dir="${workspace_dir}/analysis/gwas/emmax"
ref_dir="${workspace_dir}/analysis/gwas/reference"
pheno_dir="${workspace_dir}/analysis/phenotypes"
geno_dir="${workspace_dir}/analysis/genotype"
annot_dir="${workspace_dir}/analysis/annotation"
network_dir="${workspace_dir}/analysis/network_population"
qtl_dir="${workspace_dir}/analysis/qtl"
pop_dir="${workspace_dir}/analysis/population_genetics"

mkdir -p "${gwas_dir}" "${qtl_dir}" "${network_dir}" "${pop_dir}"

"${wrapper}" prepare_snp_reference \
  python3 "${script_dir}/prepare_snp_reference.py" \
  --bim "${geno_dir}/target_subset.bim" \
  --fai "${GWAS_FASTA_FAI_FILE:-${workspace_dir}/data/reference.fa.fai}" \
  --outdir "${ref_dir}"

"${wrapper}" summarize_emmax \
  python3 "${script_dir}/summarize_emmax_results.py" \
  --manifest "${workspace_dir}/analysis/emmax_inputs/trait_manifest.tsv" \
  --results-dir "${results_dir}" \
  --reference-dir "${ref_dir}" \
  --trait-meta "${pheno_dir}/trait_metadata.tsv" \
  --outdir "${gwas_dir}" \
  --workers "${GWAS_THREADS:-2}"

"${wrapper}" phenotype_network \
  python3 "${script_dir}/phenotype_network_population_analysis.py" \
  --raw-phenotypes "${pheno_dir}/all_traits_raw.tsv" \
  --int-phenotypes "${pheno_dir}/all_traits_int.tsv" \
  --trait-meta "${pheno_dir}/trait_metadata.tsv" \
  --trait-summaries "${gwas_dir}/trait_summaries.tsv" \
  --pca "${geno_dir}/target_subset_pca.eigenvec" \
  --covariates "${geno_dir}/covariates_subset.tsv" \
  --outdir "${network_dir}"

"${wrapper}" qtl_leads \
  python3 "${script_dir}/prepare_qtl_leads.py" \
  --master-significant "${gwas_dir}/master_significant_snps.tsv.gz" \
  --outdir "${qtl_dir}" \
  --clump-kb 500

if [[ -s "${qtl_dir}/global_lead_snps.txt" ]]; then
  "${wrapper}" qtl_ld_batch \
    python3 "${script_dir}/run_ld_batch.py" \
    --lead-file "${qtl_dir}/global_lead_snps.txt" \
    --bfile-prefix "${geno_dir}/target_subset" \
    --workspace "${workspace_dir}" \
    --parallel 2
fi

"${wrapper}" qtl_integrate \
  python3 "${script_dir}/integrate_qtl_candidates.py" \
  --lead-snps "${qtl_dir}/global_lead_snps.tsv" \
  --ld-dir "${qtl_dir}/ld" \
  --master-significant "${gwas_dir}/master_significant_snps.tsv.gz" \
  --trait-summaries "${gwas_dir}/trait_summaries.tsv" \
  --annotation-dir "${annot_dir}" \
  --module-membership "${network_dir}/lipid_module_membership.tsv" \
  --outdir "${qtl_dir}"

export GWAS_WORKSPACE_DIR="${workspace_dir}"
python3 - <<'PY'
import os
import pandas as pd
from pathlib import Path
workspace = Path(os.environ["GWAS_WORKSPACE_DIR"])
sig = pd.read_csv(workspace / "analysis/gwas/emmax/master_significant_snps.tsv.gz", sep="\t")
sig["snp_id"].drop_duplicates().to_csv(workspace / "analysis/population_genetics/significant_snp_ids.txt", index=False, header=False)
PY

for keep_file in "${network_dir}"/Group*.keep.txt; do
  if [[ ! -f "${keep_file}" ]]; then
    continue
  fi
  group_name="$(basename "${keep_file}" .keep.txt)"
  plink_bin="${PLINK_BIN:-$(command -v plink 2>/dev/null || true)}"
  if [[ -z "${plink_bin}" ]]; then
    echo "PLINK_BIN is not set and plink was not found in PATH." >&2
    exit 1
  fi
  "${wrapper}" "freq_${group_name}" \
    "${plink_bin}" \
    --bfile "${geno_dir}/target_subset" \
    --extract "${pop_dir}/significant_snp_ids.txt" \
    --keep "${keep_file}" \
    --freq \
    --allow-no-sex \
    --threads "${GWAS_THREADS:-2}" \
    --out "${pop_dir}/${group_name}"
done

"${wrapper}" population_freq \
  python3 "${script_dir}/summarize_population_frequencies.py" \
  --master-significant "${gwas_dir}/master_significant_snps.tsv.gz" \
  --overall-frq "${geno_dir}/target_subset_freq.frq" \
  --group-frq-dir "${pop_dir}" \
  --qtl-regions "${qtl_dir}/qtl_regions.tsv" \
  --shared-snps "${gwas_dir}/shared_snp_stats.tsv" \
  --trait-summaries "${gwas_dir}/trait_summaries.tsv" \
  --outdir "${pop_dir}"

"${wrapper}" cleanup_ps \
  bash -lc 'find "'"${results_dir}"'" -name "*.ps" -type f -print0 | xargs -0 -r -n 1 -P 2 gzip -f'
