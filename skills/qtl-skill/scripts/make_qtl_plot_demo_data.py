#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create demo plotting tables for qtl-skill.")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    basic = outdir / "basic"
    rich = outdir / "rich"
    basic.mkdir(parents=True, exist_ok=True)
    rich.mkdir(parents=True, exist_ok=True)

    chrom_sizes = pd.DataFrame(
        [
            {"chrom": "1", "length_bp": 2_000_000},
            {"chrom": "2", "length_bp": 1_800_000},
        ]
    )

    global_rows = []
    for bp in range(100_000, 1_900_001, 100_000):
        p = 1e-8 if bp in {800_000, 850_000} else (5e-5 if bp in {1_500_000} else 0.1)
        global_rows.append({"snp_id": f"1:{bp}", "chrom": "1", "bp": bp, "p_value": p, "mlog10_p": -math.log10(p)})
    for bp in range(100_000, 1_700_001, 100_000):
        p = 2e-7 if bp in {500_000, 550_000} else 0.08
        global_rows.append({"snp_id": f"2:{bp}", "chrom": "2", "bp": bp, "p_value": p, "mlog10_p": -math.log10(p)})
    global_df = pd.DataFrame(global_rows)

    qtl_basic = pd.DataFrame(
        [
            {
                "locus_id": "L001",
                "chrom": "1",
                "lead_snp_id": "1:800000",
                "lead_bp": 800000,
                "lead_p": 1e-8,
                "qtl_start": 760000,
                "qtl_end": 880000,
                "panel_start": 650000,
                "panel_end": 950000,
                "qtl_type": "genomewide",
            }
        ]
    )

    local_basic = pd.DataFrame(
        [
            {"locus_id": "L001", "snp_id": "1:760000", "chrom": "1", "bp": 760000, "mlog10_p": 4.2, "r2": 0.3, "is_lead": False},
            {"locus_id": "L001", "snp_id": "1:800000", "chrom": "1", "bp": 800000, "mlog10_p": 8.0, "r2": 1.0, "is_lead": True},
            {"locus_id": "L001", "snp_id": "1:850000", "chrom": "1", "bp": 850000, "mlog10_p": 7.2, "r2": 0.8, "is_lead": False},
            {"locus_id": "L001", "snp_id": "1:900000", "chrom": "1", "bp": 900000, "mlog10_p": 3.8, "r2": 0.2, "is_lead": False},
        ]
    )

    gene_basic = pd.DataFrame(
        [
            {"gene_id": "GeneA", "chrom": "1", "start_bp": 700000, "end_bp": 780000, "strand": "+"},
            {"gene_id": "GeneB", "chrom": "1", "start_bp": 795000, "end_bp": 860000, "strand": "-"},
            {"gene_id": "GeneC", "chrom": "1", "start_bp": 880000, "end_bp": 930000, "strand": "+"},
        ]
    )

    qtl_rich = pd.DataFrame(
        [
            {
                "locus_id": "L001",
                "chrom": "1",
                "lead_snp_id": "1:800000",
                "lead_bp": 800000,
                "lead_p": 1e-8,
                "qtl_start": 760000,
                "qtl_end": 880000,
                "panel_start": 650000,
                "panel_end": 950000,
                "qtl_type": "genomewide",
            },
            {
                "locus_id": "L002",
                "chrom": "2",
                "lead_snp_id": "2:500000",
                "lead_bp": 500000,
                "lead_p": 2e-7,
                "qtl_start": 450000,
                "qtl_end": 600000,
                "panel_start": 350000,
                "panel_end": 700000,
                "qtl_type": "suggestive",
            },
        ]
    )

    local_rich = pd.concat(
        [
            local_basic,
            pd.DataFrame(
                [
                    {"locus_id": "L002", "snp_id": "2:420000", "chrom": "2", "bp": 420000, "mlog10_p": 3.5, "r2": 0.2, "is_lead": False},
                    {"locus_id": "L002", "snp_id": "2:500000", "chrom": "2", "bp": 500000, "mlog10_p": 6.7, "r2": 1.0, "is_lead": True},
                    {"locus_id": "L002", "snp_id": "2:550000", "chrom": "2", "bp": 550000, "mlog10_p": 6.1, "r2": 0.75, "is_lead": False},
                    {"locus_id": "L002", "snp_id": "2:620000", "chrom": "2", "bp": 620000, "mlog10_p": 3.1, "r2": 0.1, "is_lead": False},
                ]
            ),
        ],
        ignore_index=True,
    )

    gene_rich = pd.concat(
        [
            gene_basic,
            pd.DataFrame(
                [
                    {"gene_id": "GeneD", "chrom": "2", "start_bp": 430000, "end_bp": 470000, "strand": "+"},
                    {"gene_id": "GeneE", "chrom": "2", "start_bp": 495000, "end_bp": 545000, "strand": "-"},
                    {"gene_id": "GeneF", "chrom": "2", "start_bp": 560000, "end_bp": 650000, "strand": "+"},
                ]
            ),
        ],
        ignore_index=True,
    )

    highlight_rich = pd.DataFrame(
        [
            {"locus_id": "L001", "gene_id": "GeneB", "rank": 1, "label": "GeneB"},
            {"locus_id": "L002", "gene_id": "GeneE", "rank": 1, "label": "GeneE"},
            {"locus_id": "L002", "gene_id": "GeneF", "rank": 2, "label": "GeneF"},
        ]
    )

    for target in [basic, rich]:
        chrom_sizes.to_csv(target / "chrom_sizes.tsv", sep="\t", index=False)
        global_df.to_csv(target / "global_manhattan.tsv", sep="\t", index=False)

    qtl_basic.to_csv(basic / "qtl_regions.tsv", sep="\t", index=False)
    local_basic.to_csv(basic / "local_manhattan.tsv", sep="\t", index=False)
    gene_basic.to_csv(basic / "gene_models.tsv", sep="\t", index=False)

    qtl_rich.to_csv(rich / "qtl_regions.tsv", sep="\t", index=False)
    local_rich.to_csv(rich / "local_manhattan.tsv", sep="\t", index=False)
    gene_rich.to_csv(rich / "gene_models.tsv", sep="\t", index=False)
    highlight_rich.to_csv(rich / "highlight_genes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
