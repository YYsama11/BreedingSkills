#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare trait-wise lead SNPs from significant GWAS hits.")
    parser.add_argument("--master-significant", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--clump-kb", type=int, default=500)
    return parser.parse_args()


def select_leads(df: pd.DataFrame, clump_bp: int) -> pd.DataFrame:
    selected = []
    for chrom, sub in df.sort_values("p").groupby("chrom"):
        chosen_positions = []
        for _, row in sub.iterrows():
            pos = int(row["pos"])
            if any(abs(pos - existing) <= clump_bp for existing in chosen_positions):
                continue
            selected.append(row)
            chosen_positions.append(pos)
    if not selected:
        return pd.DataFrame(columns=df.columns)
    return pd.DataFrame(selected)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sig = pd.read_csv(args.master_significant, sep="\t")
    if sig.empty:
        empty = pd.DataFrame(columns=sig.columns)
        empty.to_csv(outdir / "global_lead_snps.tsv", sep="\t", index=False)
        empty.to_csv(outdir / "unique_significant_snps.tsv", sep="\t", index=False)
        return

    unique = (
        sig.groupby(["snp_id", "chrom", "pos"], as_index=False)
        .agg(min_p=("p", "min"), trait_count=("trait_id", "nunique"), trait_ids=("trait_id", lambda values: ",".join(sorted(set(values)))))
        .sort_values(["min_p", "chrom", "pos"])
    )
    unique = unique.sort_values("min_p")
    unique.to_csv(outdir / "unique_significant_snps.tsv", sep="\t", index=False)

    clump_bp = args.clump_kb * 1000
    global_leads = []
    for chrom, sub in unique.groupby("chrom"):
        chosen_positions = []
        for _, row in sub.iterrows():
            pos = int(row["pos"])
            if any(abs(pos - existing) <= clump_bp for existing in chosen_positions):
                continue
            chosen_positions.append(pos)
            global_leads.append(row)
    global_lead_df = pd.DataFrame(global_leads).sort_values(["chrom", "pos"])
    global_lead_df.to_csv(outdir / "global_lead_snps.tsv", sep="\t", index=False)
    global_lead_df["snp_id"].to_csv(outdir / "global_lead_snps.txt", sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
