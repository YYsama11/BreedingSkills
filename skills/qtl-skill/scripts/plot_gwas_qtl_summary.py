#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec


@dataclass
class GeneModel:
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot layered GWAS-QTL summary figures.")
    parser.add_argument("--trait-label", default="QTL summary")
    parser.add_argument("--global-manhattan", required=True)
    parser.add_argument("--qtl-regions", required=True)
    parser.add_argument("--local-manhattan", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--chrom-sizes", default="")
    parser.add_argument("--gene-models", default="")
    parser.add_argument("--highlight-genes", default="")
    parser.add_argument("--y-min", type=float, default=3.0)
    parser.add_argument("--max-loci", type=int, default=0)
    parser.add_argument("--local-flank-bp", type=int, default=200000)
    return parser.parse_args()


def read_tsv(path: str) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t")


def normalize_inputs(args: argparse.Namespace):
    global_df = read_tsv(args.global_manhattan)
    qtl_df = read_tsv(args.qtl_regions)
    local_df = read_tsv(args.local_manhattan)
    gene_df = read_tsv(args.gene_models)
    highlight_df = read_tsv(args.highlight_genes)
    chrom_sizes = read_tsv(args.chrom_sizes)

    required_global = {"snp_id", "chrom", "bp"}
    required_qtl = {"locus_id", "chrom", "lead_snp_id", "lead_bp", "qtl_start", "qtl_end"}
    required_local = {"locus_id", "snp_id", "chrom", "bp"}

    if not required_global.issubset(global_df.columns):
        raise ValueError("global_manhattan must contain snp_id, chrom, bp")
    if not required_qtl.issubset(qtl_df.columns):
        raise ValueError("qtl_regions must contain locus_id, chrom, lead_snp_id, lead_bp, qtl_start, qtl_end")
    if not required_local.issubset(local_df.columns):
        raise ValueError("local_manhattan must contain locus_id, snp_id, chrom, bp")

    if "mlog10_p" not in global_df.columns:
        if "p_value" not in global_df.columns:
            raise ValueError("global_manhattan must contain either mlog10_p or p_value")
        global_df["mlog10_p"] = -np.log10(np.clip(global_df["p_value"].astype(float), 1e-300, 1.0))
    if "mlog10_p" not in local_df.columns:
        if "p_value" not in local_df.columns:
            raise ValueError("local_manhattan must contain either mlog10_p or p_value")
        local_df["mlog10_p"] = -np.log10(np.clip(local_df["p_value"].astype(float), 1e-300, 1.0))
    if "r2" not in local_df.columns:
        local_df["r2"] = np.nan
    if "is_lead" not in local_df.columns:
        local_df["is_lead"] = False
    if "lead_p" not in qtl_df.columns and "lead_mlog10_p" in qtl_df.columns:
        qtl_df["lead_p"] = np.power(10.0, -qtl_df["lead_mlog10_p"].astype(float))
    if "lead_p" not in qtl_df.columns:
        qtl_df["lead_p"] = np.nan
    if "panel_start" not in qtl_df.columns:
        qtl_df["panel_start"] = np.maximum(1, qtl_df["qtl_start"].astype(int) - args.local_flank_bp)
    if "panel_end" not in qtl_df.columns:
        qtl_df["panel_end"] = qtl_df["qtl_end"].astype(int) + args.local_flank_bp
    if "qtl_type" not in qtl_df.columns:
        qtl_df["qtl_type"] = "qtl"

    global_df["chrom"] = global_df["chrom"].astype(str)
    local_df["chrom"] = local_df["chrom"].astype(str)
    qtl_df["chrom"] = qtl_df["chrom"].astype(str)

    if chrom_sizes.empty:
        chrom_sizes = pd.concat(
            [
                global_df.groupby("chrom", as_index=False)["bp"].max().rename(columns={"bp": "length_bp"}),
                qtl_df.groupby("chrom", as_index=False)["panel_end"].max().rename(columns={"panel_end": "length_bp"}),
            ],
            ignore_index=True,
        ).groupby("chrom", as_index=False)["length_bp"].max()
    chrom_sizes["chrom"] = chrom_sizes["chrom"].astype(str)
    chrom_sizes = chrom_sizes.sort_values("chrom").reset_index(drop=True)

    if not gene_df.empty:
        expected = {"gene_id", "chrom", "start_bp", "end_bp", "strand"}
        if not expected.issubset(gene_df.columns):
            raise ValueError("gene_models must contain gene_id, chrom, start_bp, end_bp, strand")
        gene_df["chrom"] = gene_df["chrom"].astype(str)

    return chrom_sizes, global_df, qtl_df, local_df, gene_df, highlight_df


def chrom_offsets(chrom_sizes: pd.DataFrame):
    offsets = {}
    centers = {}
    cumulative = 0
    for _, row in chrom_sizes.iterrows():
        chrom = str(row["chrom"])
        offsets[chrom] = cumulative
        centers[chrom] = cumulative + float(row["length_bp"]) / 2
        cumulative += float(row["length_bp"])
    return offsets, centers


def lanes(genes: list[GeneModel], pad: int = 15000) -> dict[str, int]:
    if not genes:
        return {}
    genes = sorted(genes, key=lambda g: (g.start, g.end, g.gene_id))
    ends: list[int] = []
    out: dict[str, int] = {}
    for gene in genes:
        lane = 0
        while lane < len(ends) and gene.start <= ends[lane] + pad:
            lane += 1
        if lane == len(ends):
            ends.append(gene.end)
        else:
            ends[lane] = gene.end
        out[gene.gene_id] = lane
    return out


def r2_bin(r2: float) -> str:
    if pd.isna(r2):
        return "#bdbdbd"
    if r2 >= 0.8:
        return "#b2182b"
    if r2 >= 0.6:
        return "#ef8a62"
    if r2 >= 0.4:
        return "#fddbc7"
    if r2 >= 0.2:
        return "#67a9cf"
    return "#d1d1d1"


def alpha_col(color: str, alpha: float) -> tuple[float, float, float, float]:
    from matplotlib.colors import to_rgba

    r, g, b, _ = to_rgba(color)
    return (r, g, b, alpha)


def draw_global_panel(ax, trait_label, global_df, qtl_df, chrom_sizes, offsets, centers, y_min):
    base = global_df.copy()
    base["x"] = base["bp"] + base["chrom"].map(offsets)
    base["c"] = base["chrom"].map(lambda c: "#4C78A8" if (list(chrom_sizes["chrom"]).index(c) % 2 == 0) else "#F58518")
    ymax = max(base["mlog10_p"].max(), np.nanmax(-np.log10(np.clip(qtl_df["lead_p"].astype(float), 1e-300, 1.0))) if "lead_p" in qtl_df.columns else y_min)
    ax.scatter(base["x"], base["mlog10_p"], s=6, c=base["c"], alpha=0.35, linewidths=0, rasterized=True)
    locus_cols = ["#D73027", "#4575B4", "#1A9850", "#984EA3", "#FF7F00"]
    for idx, (_, locus) in enumerate(qtl_df.iterrows()):
        color = locus_cols[idx % len(locus_cols)]
        qtl_df.loc[qtl_df.index[idx], "plot_color"] = color
        x1 = offsets[locus["chrom"]] + locus["qtl_start"]
        x2 = offsets[locus["chrom"]] + locus["qtl_end"]
        ax.axvspan(x1, x2, color=color, alpha=0.12)
        y = -np.log10(max(float(locus["lead_p"]), 1e-300)) if pd.notna(locus["lead_p"]) else y_min
        ax.scatter([offsets[locus["chrom"]] + locus["lead_bp"]], [y], s=38, c=color, edgecolors="black", linewidths=0.4, zorder=5)
        ax.text(offsets[locus["chrom"]] + locus["lead_bp"], y + 0.03 * max(ymax - y_min, 1), locus["locus_id"], fontsize=8, ha="center")
    ax.set_xticks([centers[c] for c in chrom_sizes["chrom"]])
    ax.set_xticklabels(chrom_sizes["chrom"], fontsize=8)
    ax.set_ylabel("-log10(P)")
    ax.set_ylim(y_min, max(y_min + 0.1, ymax * 1.08))
    ax.set_title(trait_label, fontsize=11)


def genes_in_region(gene_df: pd.DataFrame, chrom: str, start: int, end: int) -> list[GeneModel]:
    if gene_df.empty:
        return []
    sub = gene_df[(gene_df["chrom"] == chrom) & (gene_df["start_bp"] <= end) & (gene_df["end_bp"] >= start)]
    return [GeneModel(row["gene_id"], str(row["chrom"]), int(row["start_bp"]), int(row["end_bp"]), str(row["strand"])) for _, row in sub.iterrows()]


def pick_highlight_ids(highlight_df: pd.DataFrame, locus_id: str) -> dict[str, str]:
    if highlight_df.empty:
        return {}
    sub = highlight_df[highlight_df["locus_id"] == locus_id]
    if sub.empty:
        return {}
    labels = {}
    for _, row in sub.iterrows():
        labels[str(row["gene_id"])] = str(row["label"]) if "label" in sub.columns and pd.notna(row.get("label", None)) else str(row["gene_id"])
    return labels


def draw_local_panel(sax, gax, locus, local_df, gene_df, highlight_df, y_min):
    panel_start = max(1, int(locus["panel_start"]))
    panel_end = int(locus["panel_end"])
    local = local_df[(local_df["locus_id"] == locus["locus_id"]) & (local_df["bp"] >= panel_start) & (local_df["bp"] <= panel_end)].copy()
    if local.empty:
        local = pd.DataFrame(
            {
                "bp": [int(locus["lead_bp"])],
                "mlog10_p": [-np.log10(max(float(locus["lead_p"]), 1e-300)) if pd.notna(locus["lead_p"]) else y_min],
                "r2": [1.0],
            }
        )
    local["c"] = local["r2"].apply(r2_bin)
    ymax = max(local["mlog10_p"].max(), y_min + 0.1)
    sax.axvspan(int(locus["qtl_start"]), int(locus["qtl_end"]), color=locus["plot_color"], alpha=0.12)
    sax.scatter(local["bp"], local["mlog10_p"], s=14, c=local["c"], alpha=0.55, linewidths=0, rasterized=True)
    sax.scatter([int(locus["lead_bp"])], [-np.log10(max(float(locus["lead_p"]), 1e-300)) if pd.notna(locus["lead_p"]) else y_min], s=55, c="#111111", marker="D", edgecolors="white", linewidths=0.4, zorder=6)
    sax.set_xlim(panel_start, panel_end)
    sax.set_ylim(y_min, ymax * 1.1)
    sax.set_ylabel("-log10(P)")
    sax.set_title(
        f"{locus['locus_id']} | {locus['chrom']}:{int(locus['qtl_start']):,}-{int(locus['qtl_end']):,} | lead {locus['lead_snp_id']} ({locus['qtl_type']})",
        fontsize=10,
        loc="left",
    )

    panel_genes = genes_in_region(gene_df, str(locus["chrom"]), panel_start, panel_end)
    lane_map = lanes(panel_genes)
    max_lane = max(lane_map.values(), default=0)
    gax.axvspan(int(locus["qtl_start"]), int(locus["qtl_end"]), color=locus["plot_color"], alpha=0.12)
    high = pick_highlight_ids(highlight_df, str(locus["locus_id"]))
    for gene in panel_genes:
        y = max_lane - lane_map[gene.gene_id]
        is_highlight = gene.gene_id in high
        col = locus["plot_color"] if is_highlight else "#666666"
        lw = 2.2 if is_highlight else 1.2
        gax.plot([gene.start, gene.end], [y, y], color=col, linewidth=lw, solid_capstyle="round")
        if is_highlight:
            anchor, ha = (gene.start, "left") if gene.strand == "+" else (gene.end, "right")
            gax.text(anchor, y + 0.18, high[gene.gene_id], fontsize=7, ha=ha, va="bottom", color=col)
    if not panel_genes:
        gax.text((panel_start + panel_end) / 2, 0.2, "No gene models provided in this interval", fontsize=8, ha="center", color="#666666")
    gax.set_ylim(-0.8, max_lane + 1.2)
    gax.set_yticks([])
    gax.set_xlabel(f"{locus['chrom']} position (bp)")
    for side in ["left", "right", "top"]:
        gax.spines[side].set_visible(False)
    plt.setp(sax.get_xticklabels(), visible=False)


def save_plot(trait_label, global_df, qtl_df, local_df, chrom_sizes, gene_df, highlight_df, out_path, y_min):
    if qtl_df.empty:
        raise ValueError("No QTL regions available for plotting.")
    offsets, centers = chrom_offsets(chrom_sizes)
    heights = [3.2] + [2.0, 0.85] * len(qtl_df)
    fig = plt.figure(figsize=(14, 4.8 + len(qtl_df) * 2.8))
    gs = GridSpec(1 + len(qtl_df) * 2, 1, figure=fig, height_ratios=heights, hspace=0.35)
    ax = fig.add_subplot(gs[0, 0])
    draw_global_panel(ax, trait_label, global_df, qtl_df, chrom_sizes, offsets, centers, y_min)
    for idx, (_, locus) in enumerate(qtl_df.iterrows()):
        sax = fig.add_subplot(gs[1 + idx * 2, 0])
        gax = fig.add_subplot(gs[2 + idx * 2, 0], sharex=sax)
        draw_local_panel(sax, gax, locus, local_df, gene_df, highlight_df, y_min)
    fig.subplots_adjust(top=0.97, bottom=0.05, left=0.06, right=0.99, hspace=0.35)
    if str(out_path).lower().endswith(".pdf"):
        with PdfPages(out_path) as pdf:
            pdf.savefig(fig)
    else:
        fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main():
    args = parse_args()
    chrom_sizes, global_df, qtl_df, local_df, gene_df, highlight_df = normalize_inputs(args)
    if args.max_loci > 0 and len(qtl_df) > args.max_loci:
        qtl_df = qtl_df.iloc[: args.max_loci].copy()
    save_plot(args.trait_label, global_df, qtl_df, local_df, chrom_sizes, gene_df, highlight_df, args.out, args.y_min)


if __name__ == "__main__":
    main()
