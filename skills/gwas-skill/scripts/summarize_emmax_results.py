#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import math
import multiprocessing as mp
from functools import partial
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2


plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.dpi"] = 150

GLOBAL = {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize EMMAX GWAS outputs and create figures.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument("--trait-meta", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def init_worker(reference_dir: str) -> None:
    ref = Path(reference_dir)
    GLOBAL["chrom"] = np.load(ref / "chrom.npy", mmap_mode="r")
    GLOBAL["pos"] = np.load(ref / "pos.npy", mmap_mode="r")
    GLOBAL["cum_pos"] = np.load(ref / "cum_pos.npy", mmap_mode="r")
    GLOBAL["bin_index"] = np.load(ref / "bin_index.npy", mmap_mode="r")
    GLOBAL["ticks"] = pd.read_csv(ref / "chrom_ticks.tsv", sep="\t")
    GLOBAL["bins"] = pd.read_csv(ref / "manhattan_bins.tsv", sep="\t")
    summary = pd.read_csv(ref / "reference_summary.tsv", sep="\t")
    GLOBAL["n_snp"] = int(summary.loc[summary["metric"] == "snp_count", "value"].iloc[0])
    GLOBAL["bonf"] = float(summary.loc[summary["metric"] == "bonferroni_threshold", "value"].iloc[0])
    GLOBAL["suggestive"] = float(summary.loc[summary["metric"] == "suggestive_threshold", "value"].iloc[0])


def load_reml_value(results_dir: Path, trait_id: str) -> float | None:
    reml_file = results_dir / f"{trait_id}.reml"
    if not reml_file.exists():
        return None
    values = [float(line.strip()) for line in reml_file.read_text().splitlines() if line.strip()]
    if not values:
        return None
    return values[-1]


def make_manhattan(trait_id: str, chrom: np.ndarray, cum_pos: np.ndarray, p: np.ndarray, out_png: Path) -> None:
    bins = GLOBAL["bins"].copy()
    bin_index = np.asarray(GLOBAL["bin_index"])
    min_p = np.full(bins.shape[0], np.inf, dtype=float)
    np.minimum.at(min_p, bin_index, p)
    bins["min_p"] = min_p
    bins["neglog10p"] = -np.log10(np.clip(bins["min_p"], 1e-300, 1.0))

    sig_mask = p < GLOBAL["bonf"]
    sig_x = cum_pos[sig_mask]
    sig_y = -np.log10(np.clip(p[sig_mask], 1e-300, 1.0))
    sig_chr = chrom[sig_mask]

    cmap = plt.get_cmap("RdYlGn")
    colors = [cmap(0.15), cmap(0.85)]

    fig, ax = plt.subplots(figsize=(14, 4.8))
    for color_group, sub in bins.groupby("color_group"):
        ax.scatter(sub["bin_mid"], sub["neglog10p"], s=3, color=colors[int(color_group)], alpha=0.75, rasterized=True)

    if sig_mask.any():
        for chrom_value in np.unique(sig_chr):
            sub_mask = sig_chr == chrom_value
            ax.scatter(sig_x[sub_mask], sig_y[sub_mask], s=6, color="black", alpha=0.9, rasterized=True)

    ax.axhline(-np.log10(GLOBAL["bonf"]), color="red", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(GLOBAL["suggestive"]), color="grey", linestyle=":", linewidth=0.8)
    ax.set_xticks(GLOBAL["ticks"]["midpoint"])
    ax.set_xticklabels(GLOBAL["ticks"]["label"], rotation=0)
    ax.set_ylabel("-log10(P)")
    ax.set_xlabel("Chromosome")
    ax.set_title(trait_id)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def make_qq(trait_id: str, p: np.ndarray, lambda_gc: float, out_png: Path) -> None:
    p = np.clip(p, 1e-300, 1.0)
    step = max(1, len(p) // 200000)
    sampled = np.sort(p[::step])
    expected = -np.log10((np.arange(1, len(sampled) + 1) - 0.5) / len(sampled))
    observed = -np.log10(sampled)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(expected, observed, s=4, color=plt.get_cmap("RdYlGn")(0.85), alpha=0.7, rasterized=True)
    max_val = max(expected.max(), observed.max())
    ax.plot([0, max_val], [0, max_val], linestyle="--", color="black", linewidth=0.8)
    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    ax.set_title(f"{trait_id}\nλGC={lambda_gc:.3f}")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def process_trait(row: dict[str, str], results_dir: str, outdir: str) -> dict[str, object]:
    trait_id = row["trait_id"]
    results_dir_path = Path(results_dir)
    outdir_path = Path(outdir)
    summary_dir = outdir_path / "per_trait_summaries"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary_file = summary_dir / f"{trait_id}.summary.tsv"
    if summary_file.exists():
        return pd.read_csv(summary_file, sep="\t").iloc[0].to_dict()

    result_file = results_dir_path / f"{trait_id}.ps.gz"
    if not result_file.exists():
        alt = results_dir_path / f"{trait_id}.ps"
        if alt.exists():
            result_file = alt
        else:
            raise FileNotFoundError(f"Missing result file for {trait_id}")

    df = pd.read_csv(result_file, sep="\t", header=None, names=["snp_id", "beta", "se", "p"], compression="infer")
    p = df["p"].to_numpy(dtype=float)
    p = np.clip(p, 1e-300, 1.0)
    chrom = np.asarray(GLOBAL["chrom"])
    pos = np.asarray(GLOBAL["pos"])
    cum_pos = np.asarray(GLOBAL["cum_pos"])

    sig_mask = p < GLOBAL["bonf"]
    suggestive_mask = p < GLOBAL["suggestive"]
    neglogp = -np.log10(p)
    top_idx = int(np.argmin(p))
    chisq = chi2.isf(p, 1)
    lambda_gc = float(np.nanmedian(chisq) / 0.454936423119572)
    pseudo_h2 = load_reml_value(results_dir_path, trait_id)

    sig_dir = outdir_path / "significant"
    top_dir = outdir_path / "top_hits"
    manhattan_dir = outdir_path / "figures" / "manhattan"
    qq_dir = outdir_path / "figures" / "qq"
    for directory in [sig_dir, top_dir, manhattan_dir, qq_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    sig_file = sig_dir / f"{trait_id}.significant.tsv.gz"
    if not sig_file.exists():
        sig_df = df.loc[sig_mask, ["snp_id", "beta", "se", "p"]].copy()
        if not sig_df.empty:
            sig_df["chrom"] = chrom[sig_mask]
            sig_df["pos"] = pos[sig_mask]
            sig_df["neglog10p"] = neglogp[sig_mask]
            sig_df["trait_id"] = trait_id
            sig_df["trait_level"] = row.get("trait_level", "")
            sig_df.to_csv(sig_file, sep="\t", index=False)

    top_file = top_dir / f"{trait_id}.top1000.tsv.gz"
    if not top_file.exists():
        top_n = min(1000, len(df))
        top_indices = np.argpartition(p, top_n - 1)[:top_n]
        top_df = df.iloc[top_indices].copy()
        top_df["chrom"] = chrom[top_indices]
        top_df["pos"] = pos[top_indices]
        top_df["neglog10p"] = neglogp[top_indices]
        top_df["trait_id"] = trait_id
        top_df = top_df.sort_values("p")
        top_df.to_csv(top_file, sep="\t", index=False)

    manhattan_file = manhattan_dir / f"{trait_id}.manhattan.png"
    if not manhattan_file.exists():
        make_manhattan(trait_id, chrom, cum_pos, p, manhattan_file)

    qq_file = qq_dir / f"{trait_id}.qq.png"
    if not qq_file.exists():
        make_qq(trait_id, p, lambda_gc, qq_file)

    summary_row = {
        "trait_id": trait_id,
        "trait_level": row.get("trait_level", ""),
        "display_name": row.get("display_name", ""),
        "tested_snp_count": len(df),
        "bonferroni_threshold": GLOBAL["bonf"],
        "suggestive_threshold": GLOBAL["suggestive"],
        "significant_snp_count": int(sig_mask.sum()),
        "suggestive_snp_count": int(suggestive_mask.sum()),
        "top_snp_id": df.loc[top_idx, "snp_id"],
        "top_chrom": str(chrom[top_idx]),
        "top_pos": int(pos[top_idx]),
        "top_beta": float(df.loc[top_idx, "beta"]),
        "top_se": float(df.loc[top_idx, "se"]),
        "top_p": float(df.loc[top_idx, "p"]),
        "lambda_gc": lambda_gc,
        "pseudo_heritability": pseudo_h2,
    }
    pd.DataFrame([summary_row]).to_csv(summary_file, sep="\t", index=False)
    return summary_row


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest, sep="\t")
    trait_meta = pd.read_csv(args.trait_meta, sep="\t")
    merged = manifest.merge(trait_meta[["trait_id", "trait_level", "display_name"]], on="trait_id", how="left")
    rows = merged.to_dict(orient="records")

    with mp.Pool(processes=args.workers, initializer=init_worker, initargs=(args.reference_dir,)) as pool:
        summaries = pool.map(partial(process_trait, results_dir=args.results_dir, outdir=args.outdir), rows)

    summary_df = pd.DataFrame(summaries).sort_values(["trait_level", "trait_id"])
    summary_df = summary_df.merge(trait_meta[["trait_id", "trait_level", "display_name"]], on="trait_id", how="left", suffixes=("", "_meta"))
    for col in ["trait_level", "display_name"]:
        summary_df[col] = summary_df[col].fillna(summary_df[f"{col}_meta"])
    drop_cols = [col for col in summary_df.columns if col.endswith("_meta")]
    if drop_cols:
        summary_df = summary_df.drop(columns=drop_cols)
    summary_df = summary_df.sort_values(["trait_level", "trait_id"])
    summary_df.to_csv(outdir / "trait_summaries.tsv", sep="\t", index=False)

    sig_files = sorted((outdir / "significant").glob("*.significant.tsv.gz"))
    if sig_files:
        sig_tables = [pd.read_csv(path, sep="\t") for path in sig_files]
        master_sig = pd.concat(sig_tables, ignore_index=True)
    else:
        master_sig = pd.DataFrame(columns=["snp_id", "beta", "se", "p", "chrom", "pos", "neglog10p", "trait_id", "trait_level"])
    master_sig.to_csv(outdir / "master_significant_snps.tsv.gz", sep="\t", index=False)

    if not master_sig.empty:
        chrom_dist = master_sig.groupby("chrom").size().reset_index(name="significant_snp_count")
        chrom_dist.to_csv(outdir / "chromosome_distribution.tsv", sep="\t", index=False)

        shared = (
            master_sig.groupby("snp_id")
            .agg(
                trait_count=("trait_id", "nunique"),
                trait_levels=("trait_level", lambda values: ",".join(sorted(set(map(str, values))))),
                chrom=("chrom", "first"),
                pos=("pos", "first"),
                min_p=("p", "min"),
            )
            .reset_index()
            .sort_values(["trait_count", "min_p"], ascending=[False, True])
        )
        shared.to_csv(outdir / "shared_snp_stats.tsv", sep="\t", index=False)
        shared.loc[shared["trait_count"] >= 2].to_csv(outdir / "pleiotropic_snps.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(columns=["chrom", "significant_snp_count"]).to_csv(outdir / "chromosome_distribution.tsv", sep="\t", index=False)
        pd.DataFrame(columns=["snp_id", "trait_count", "trait_levels", "chrom", "pos", "min_p"]).to_csv(outdir / "shared_snp_stats.tsv", sep="\t", index=False)
        pd.DataFrame(columns=["snp_id", "trait_count", "trait_levels", "chrom", "pos", "min_p"]).to_csv(outdir / "pleiotropic_snps.tsv", sep="\t", index=False)

    check = {
        "trait_count": int(summary_df.shape[0]),
        "result_directory": str(Path(args.results_dir)),
        "note": "Manhattan plots use binned genome-wide background plus full significant-point overlay to keep the workflow tractable for 1,534 traits.",
    }
    (outdir / "summary_check.json").write_text(json.dumps(check, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
