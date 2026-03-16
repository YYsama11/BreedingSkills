#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: bash scripts/run_qtl_candidate_gene.sh --config path/to/qtl_candidate_gene.env
EOF
}

config=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      config="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

[[ -n "$config" ]] || { usage; exit 1; }
[[ -f "$config" ]] || { echo "Config not found: $config" >&2; exit 1; }

set -a
source "$config"
set +a

required_vars=(GWAS_TABLE GENE_ANNOTATION OUTDIR P_THRESHOLD WINDOW_BP)
for var_name in "${required_vars[@]}"; do
  [[ -n "${!var_name:-}" ]] || { echo "Missing variable: $var_name" >&2; exit 1; }
done

mkdir -p "$OUTDIR/qtl_regions" "$OUTDIR/candidates" "$OUTDIR/logs"

awk -F'\t' -v p_threshold="$P_THRESHOLD" '
NR==1 {
  if ($1 != "trait" || $2 != "chrom" || $3 != "pos" || $4 != "snp_id" || $5 != "pvalue") {
    print "gwas table header must start with: trait chrom pos snp_id pvalue" > "/dev/stderr"
    exit 2
  }
  print $0
  next
}
$5 + 0 <= p_threshold {
  print $0
}
' "$GWAS_TABLE" > "$OUTDIR/significant_hits.tsv"

cat > "$OUTDIR/workflow_plan.txt" <<EOF
1. Sort ${OUTDIR}/significant_hits.tsv by trait, chrom, and pvalue.
2. Collapse nearby hits within ${WINDOW_BP} bp to provisional lead SNPs.
3. If ${LD_TABLE:-<not_provided>} exists, replace fixed windows with LD-supported boundaries.
4. Intersect QTL intervals with ${GENE_ANNOTATION} to create candidate genes.
5. Export QTL summary and candidate gene tables under ${OUTDIR}/qtl_regions and ${OUTDIR}/candidates.
EOF

printf 'Prepared %s and %s\n' "$OUTDIR/significant_hits.tsv" "$OUTDIR/workflow_plan.txt"
