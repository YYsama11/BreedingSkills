#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


def main() -> None:
    outdir = Path(__file__).resolve().parent
    rows = []
    for chrom in ["1", "2"]:
        max_bp = 2_000_000 if chrom == "1" else 1_800_000
        for bp in range(100_000, max_bp + 1, 100_000):
            p = 0.1
            if chrom == "1" and bp == 800_000:
                p = 1e-8
            elif chrom == "2" and bp == 500_000:
                p = 2e-7
            elif bp % 300_000 == 0:
                p = 5e-5
            rows.append({"snp_id": f"{chrom}:{bp}", "chrom": chrom, "bp": bp, "p_value": p})
    pd.DataFrame(rows).to_csv(outdir / "demo_gwas.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"chrom": "1", "length_bp": 2_000_000}, {"chrom": "2", "length_bp": 1_800_000}]
    ).to_csv(outdir / "demo_chrom_sizes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
