#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare EMMAX phenotype and covariate inputs.")
    parser.add_argument("--trait-matrix", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--trait-meta", required=False)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    pheno_dir = outdir / "phenotypes"
    pheno_dir.mkdir(parents=True, exist_ok=True)

    matrix = pd.read_csv(args.trait_matrix, sep="\t")
    cov = pd.read_csv(args.covariates, sep="\t")

    required_cols = {"FID", "IID"}
    if not required_cols.issubset(matrix.columns):
        raise ValueError("Trait matrix must contain FID and IID columns.")
    if not required_cols.issubset(cov.columns):
        raise ValueError("Covariate file must contain FID and IID columns.")

    if args.trait_meta and Path(args.trait_meta).exists():
        trait_meta = pd.read_csv(args.trait_meta, sep="\t")
    else:
        trait_ids = [col for col in matrix.columns if col not in {"FID", "IID"}]
        trait_meta = pd.DataFrame(
            {
                "trait_id": trait_ids,
                "trait_level": ["trait"] * len(trait_ids),
                "display_name": trait_ids,
            }
        )

    sample_order = matrix[["FID", "IID"]].copy()
    cov.to_csv(outdir / "covariates_emmax.tsv", sep="\t", index=False, header=False)
    sample_order.to_csv(outdir / "sample_order.tsv", sep="\t", index=False)

    trait_ids = [col for col in matrix.columns if col not in {"FID", "IID"}]
    manifest = []
    for trait_id in trait_ids:
        pheno = sample_order.copy()
        pheno["trait"] = matrix[trait_id]
        pheno_path = pheno_dir / f"{trait_id}.tsv"
        pheno.to_csv(pheno_path, sep="\t", index=False, header=False)
        meta_match = trait_meta.loc[trait_meta["trait_id"] == trait_id]
        meta_row = meta_match.iloc[0].to_dict() if not meta_match.empty else {"trait_level": "trait", "display_name": trait_id}
        manifest.append(
            {
                "trait_id": trait_id,
                "trait_level": meta_row.get("trait_level", "trait"),
                "display_name": meta_row.get("display_name", trait_id),
                "phenotype_file": str(pheno_path),
            }
        )

    pd.DataFrame(manifest).to_csv(outdir / "trait_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
