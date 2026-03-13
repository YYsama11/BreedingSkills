#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare EMMAX phenotype and covariate inputs.")
    parser.add_argument("--trait-matrix", required=True)
    parser.add_argument("--trait-meta", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    pheno_dir = outdir / "phenotypes"
    pheno_dir.mkdir(parents=True, exist_ok=True)

    matrix = pd.read_csv(args.trait_matrix, sep="\t")
    trait_meta = pd.read_csv(args.trait_meta, sep="\t")
    cov = pd.read_csv(args.covariates, sep="\t")

    sample_order = matrix[["FID", "IID"]].copy()
    cov_3col = cov.copy()
    cov_3col.to_csv(outdir / "covariates_emmax.tsv", sep="\t", index=False, header=False)
    sample_order.to_csv(outdir / "sample_order.tsv", sep="\t", index=False)

    trait_ids = [col for col in matrix.columns if col not in {"FID", "IID"}]
    manifest = []
    for trait_id in trait_ids:
        pheno = sample_order.copy()
        pheno["trait"] = matrix[trait_id]
        pheno_path = pheno_dir / f"{trait_id}.tsv"
        pheno.to_csv(pheno_path, sep="\t", index=False, header=False)
        meta_row = trait_meta.loc[trait_meta["trait_id"] == trait_id].iloc[0].to_dict()
        manifest.append(
            {
                "trait_id": trait_id,
                "trait_level": meta_row.get("trait_level", ""),
                "display_name": meta_row.get("display_name", ""),
                "phenotype_file": str(pheno_path),
            }
        )

    pd.DataFrame(manifest).to_csv(outdir / "trait_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
