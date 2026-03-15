#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build plotting tables for GWAS-QTL summary figures.")
    parser.add_argument("--qtl-regions", required=True)
    parser.add_argument("--gwas-significant", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--gwas-all", required=False)
    parser.add_argument("--fai", required=False)
    parser.add_argument("--ld-dir", required=False)
    parser.add_argument("--candidate-genes", required=False)
    parser.add_argument("--representative-genes", required=False)
    return parser.parse_args()


def read_table(path: str | None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    file_path = Path(path)
    if not file_path.exists() or file_path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(file_path, sep="\t")


def load_chrom_sizes(args: argparse.Namespace, qtl: pd.DataFrame, global_df: pd.DataFrame) -> pd.DataFrame:
    if args.fai and Path(args.fai).exists():
        fai = pd.read_csv(args.fai, sep="\t", header=None, names=["chrom", "length_bp", "offset", "linebases", "linewidth"])
        return fai.loc[:, ["chrom", "length_bp"]]
    candidates = []
    if not global_df.empty:
        candidates.append(global_df.groupby("chrom", as_index=False)["bp"].max().rename(columns={"bp": "length_bp"}))
    if not qtl.empty:
        candidates.append(qtl.groupby("chrom", as_index=False)["panel_end"].max().rename(columns={"panel_end": "length_bp"}))
    if not candidates:
        raise ValueError("Unable to infer chromosome sizes without FAI or coordinate tables.")
    merged = pd.concat(candidates, ignore_index=True).groupby("chrom", as_index=False)["length_bp"].max()
    return merged


def build_local(qtl: pd.DataFrame, sig: pd.DataFrame, ld_dir: str | None) -> pd.DataFrame:
    rows = []
    for _, locus in qtl.iterrows():
        sub = sig[(sig["chrom"].astype(str) == str(locus["chrom"])) & (sig["pos"] >= locus["panel_start"]) & (sig["pos"] <= locus["panel_end"])].copy()
        if sub.empty:
            continue
        sub["locus_id"] = locus["locus_id"]
        sub["bp"] = sub["pos"]
        sub["snp_id"] = sub["snp_id"]
        sub["mlog10_p"] = -np.log10(sub["p"].clip(lower=1e-300))
        sub["r2"] = pd.NA
        sub["is_lead"] = sub["snp_id"] == locus["lead_snp_id"]
        if ld_dir:
            ld_file = Path(ld_dir) / f"{str(locus['lead_snp_id']).replace(':', '_')}.ld"
            if ld_file.exists() and ld_file.stat().st_size > 0:
                ld = pd.read_csv(ld_file, sep=r"\s+")
                if {"SNP_B", "R2"}.issubset(ld.columns):
                    r2_map = ld.set_index("SNP_B")["R2"].to_dict()
                    sub["r2"] = sub["snp_id"].map(r2_map)
                    sub.loc[sub["snp_id"] == locus["lead_snp_id"], "r2"] = 1.0
        rows.append(sub.loc[:, ["locus_id", "snp_id", "chrom", "bp", "mlog10_p", "r2", "is_lead"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["locus_id", "snp_id", "chrom", "bp", "mlog10_p", "r2", "is_lead"])


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    qtl = read_table(args.qtl_regions)
    sig = read_table(args.gwas_significant)
    gwas_all = read_table(args.gwas_all) if args.gwas_all else pd.DataFrame()
    candidate = read_table(args.candidate_genes) if args.candidate_genes else pd.DataFrame()
    representative = read_table(args.representative_genes) if args.representative_genes else pd.DataFrame()

    if qtl.empty:
        raise ValueError("QTL regions table is empty.")
    if sig.empty:
        raise ValueError("GWAS significant table is empty.")

    rename_map = {}
    if "snp_id" not in qtl.columns and "lead_snp" in qtl.columns:
        rename_map["lead_snp"] = "lead_snp_id"
    if "snp_id" in qtl.columns and "lead_snp_id" not in qtl.columns:
        rename_map["snp_id"] = "lead_snp_id"
    if rename_map:
        qtl = qtl.rename(columns=rename_map)
    if "lead_bp" not in qtl.columns and "pos" in qtl.columns:
        qtl["lead_bp"] = qtl["pos"]
    if "panel_start" not in qtl.columns:
        qtl["panel_start"] = qtl["qtl_start"]
    if "panel_end" not in qtl.columns:
        qtl["panel_end"] = qtl["qtl_end"]
    if "qtl_type" not in qtl.columns:
        qtl["qtl_type"] = "qtl"
    if "locus_id" not in qtl.columns and "qtl_id" in qtl.columns:
        qtl["locus_id"] = qtl["qtl_id"]

    if gwas_all.empty:
        global_df = sig.loc[:, ["snp_id", "chrom", "pos", "p"]].drop_duplicates().rename(columns={"pos": "bp", "p": "p_value"})
    else:
        global_df = gwas_all.copy()
        if "bp" not in global_df.columns and "pos" in global_df.columns:
            global_df["bp"] = global_df["pos"]
        if "p_value" not in global_df.columns and "p" in global_df.columns:
            global_df["p_value"] = global_df["p"]
        global_df = global_df.loc[:, ["snp_id", "chrom", "bp", "p_value"]].drop_duplicates()
    global_df["mlog10_p"] = -np.log10(global_df["p_value"].clip(lower=1e-300))

    chrom_sizes = load_chrom_sizes(args, qtl, global_df)
    local_df = build_local(qtl, sig, args.ld_dir)

    if not candidate.empty:
        gene_models = candidate.loc[:, [col for col in ["gene_id", "chrom", "start", "end"] if col in candidate.columns]].drop_duplicates().copy()
        if "start" in gene_models.columns:
            gene_models = gene_models.rename(columns={"start": "start_bp", "end": "end_bp"})
        if "strand" not in gene_models.columns:
            gene_models["strand"] = "+"
    else:
        gene_models = pd.DataFrame(columns=["gene_id", "chrom", "start_bp", "end_bp", "strand"])

    if not representative.empty:
        highlight = representative.copy()
        if "representative_gene_id" in highlight.columns and "gene_id" not in highlight.columns:
            highlight = highlight.rename(columns={"representative_gene_id": "gene_id"})
        if "label" not in highlight.columns:
            highlight["label"] = highlight["gene_id"]
        if "rank" not in highlight.columns:
            highlight["rank"] = 1
        highlight = highlight.loc[:, ["qtl_id", "gene_id", "rank", "label"]]
        highlight = highlight.rename(columns={"qtl_id": "locus_id"})
    else:
        highlight = pd.DataFrame(columns=["locus_id", "gene_id", "rank", "label"])

    chrom_sizes.to_csv(outdir / "chrom_sizes.tsv", sep="\t", index=False)
    global_df.to_csv(outdir / "global_manhattan.tsv", sep="\t", index=False)
    qtl.loc[:, ["locus_id", "chrom", "lead_snp_id", "lead_bp", "lead_p", "qtl_start", "qtl_end", "panel_start", "panel_end", "qtl_type"]].to_csv(outdir / "qtl_regions.tsv", sep="\t", index=False)
    local_df.to_csv(outdir / "local_manhattan.tsv", sep="\t", index=False)
    gene_models.to_csv(outdir / "gene_models.tsv", sep="\t", index=False)
    highlight.to_csv(outdir / "highlight_genes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
