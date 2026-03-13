#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize significant SNP allele frequencies by subgroup.")
    parser.add_argument("--master-significant", required=True)
    parser.add_argument("--overall-frq", required=True)
    parser.add_argument("--group-frq-dir", required=True)
    parser.add_argument("--qtl-regions", required=True)
    parser.add_argument("--shared-snps", required=True)
    parser.add_argument("--trait-summaries", required=False)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def read_overall_freq(freq_file: Path, snps: set[str]) -> pd.DataFrame:
    chunks = []
    for chunk in pd.read_csv(freq_file, sep=r"\s+", chunksize=500000):
        sub = chunk[chunk["SNP"].isin(snps)].copy()
        if not sub.empty:
            chunks.append(sub)
    if not chunks:
        return pd.DataFrame(columns=["SNP", "MAF", "NCHROBS"])
    out = pd.concat(chunks, ignore_index=True)
    return out.loc[:, ["SNP", "MAF", "NCHROBS"]]


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    fig_dir = outdir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    sig = pd.read_csv(args.master_significant, sep="\t", usecols=["snp_id", "chrom", "pos", "trait_id", "trait_level", "p"])
    if args.trait_summaries and Path(args.trait_summaries).exists():
        summaries = pd.read_csv(args.trait_summaries, sep="\t")
        trait_level_map = summaries.set_index("trait_id")["trait_level"].to_dict()
        sig["trait_level"] = sig["trait_id"].map(trait_level_map).fillna(sig.get("trait_level")).fillna("unknown")
    shared = pd.read_csv(args.shared_snps, sep="\t")
    qtl = pd.read_csv(args.qtl_regions, sep="\t")
    sig_unique = (
        sig.groupby(["snp_id", "chrom", "pos"], as_index=False)
        .agg(
            min_p=("p", "min"),
            trait_count=("trait_id", "nunique"),
            trait_levels=("trait_level", lambda values: ",".join(sorted(set(map(str, values))))),
        )
    )
    unique_snps = set(sig_unique["snp_id"].unique())

    overall = read_overall_freq(Path(args.overall_frq), unique_snps).rename(columns={"SNP": "snp_id", "MAF": "overall_maf"})
    merged = sig_unique.merge(overall, on="snp_id", how="left")

    group_tables = []
    for frq_file in sorted(Path(args.group_frq_dir).glob("Group*.frq")):
        group_name = frq_file.stem
        frq = pd.read_csv(frq_file, sep=r"\s+")
        frq = frq.loc[:, ["SNP", "MAF"]].rename(columns={"SNP": "snp_id", "MAF": f"{group_name}_maf"})
        group_tables.append(frq)
    for table in group_tables:
        merged = merged.merge(table, on="snp_id", how="left")

    group_cols = [col for col in merged.columns if col.endswith("_maf") and col != "overall_maf"]
    if group_cols:
        merged["max_group_maf_diff"] = merged[group_cols].max(axis=1) - merged[group_cols].min(axis=1)
    else:
        merged["max_group_maf_diff"] = np.nan
    merged.to_csv(outdir / "significant_snp_allele_frequencies.tsv", sep="\t", index=False)

    if not merged.empty and group_cols:
        threshold = merged["max_group_maf_diff"].quantile(0.95)
        diff_proxy = merged[merged["max_group_maf_diff"] >= threshold].copy()
    else:
        threshold = np.nan
        diff_proxy = merged.iloc[0:0].copy()
    diff_proxy.to_csv(outdir / "selection_proxy_snps.tsv", sep="\t", index=False)

    if not diff_proxy.empty and not qtl.empty:
        overlap_rows = []
        for _, snp in diff_proxy.iterrows():
            overlap = qtl[(qtl["chrom"] == snp["chrom"]) & (qtl["qtl_start"] <= snp["pos"]) & (qtl["qtl_end"] >= snp["pos"])]
            if overlap.empty:
                continue
            cols = ["qtl_id", "chrom", "qtl_start", "qtl_end"]
            if "trait_id" in overlap.columns:
                cols.append("trait_id")
            if "trait_level" in overlap.columns:
                cols.append("trait_level")
            if "trait_ids" in overlap.columns:
                cols.append("trait_ids")
            if "trait_count" in overlap.columns:
                cols.append("trait_count")
            tmp = overlap.loc[:, cols].copy()
            tmp["snp_id"] = snp["snp_id"]
            tmp["pos"] = snp["pos"]
            tmp["max_group_maf_diff"] = snp["max_group_maf_diff"]
            overlap_rows.append(tmp)
        overlap_df = pd.concat(overlap_rows, ignore_index=True) if overlap_rows else pd.DataFrame(columns=["qtl_id", "chrom", "qtl_start", "qtl_end", "snp_id", "pos", "max_group_maf_diff"])
    else:
        overlap_df = pd.DataFrame(columns=["qtl_id", "trait_id", "trait_level", "chrom", "qtl_start", "qtl_end", "snp_id", "pos", "max_group_maf_diff"])
    overlap_df.to_csv(outdir / "selection_proxy_qtl_overlap.tsv", sep="\t", index=False)

    if group_cols and not shared.empty:
        top_shared = shared.sort_values(["trait_count", "min_p"], ascending=[False, True]).head(50)["snp_id"].tolist()
        heatmap_df = merged[merged["snp_id"].isin(top_shared)].drop_duplicates("snp_id").set_index("snp_id")[group_cols]
        if not heatmap_df.empty:
            fig, ax = plt.subplots(figsize=(10, max(4, len(heatmap_df) * 0.2)))
            sns.heatmap(heatmap_df, cmap="RdYlGn", ax=ax)
            ax.set_title("Top pleiotropic SNP subgroup MAF")
            fig.tight_layout()
            fig.savefig(fig_dir / "top_pleiotropic_snp_group_maf_heatmap.png")
            plt.close(fig)

    check = {
        "unique_significant_snp_count": int(sig["snp_id"].nunique()),
        "selection_proxy_threshold": None if pd.isna(threshold) else float(threshold),
        "selection_proxy_overlap_count": int(overlap_df.shape[0]),
    }
    (outdir / "population_frequency_check.json").write_text(json.dumps(check, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
