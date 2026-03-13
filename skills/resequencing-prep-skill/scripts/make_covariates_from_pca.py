#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build an EMMAX-style covariate table from PLINK PCA output.")
    parser.add_argument("--eigenvec", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--num-pcs", type=int, default=3)
    parser.add_argument("--include-intercept", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    eigenvec = pd.read_csv(args.eigenvec, sep=r"\s+", header=None)
    if eigenvec.shape[1] < 2 + args.num_pcs:
        raise ValueError("Not enough principal components available in eigenvec file.")
    eigenvec.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, eigenvec.shape[1] - 1)]
    cols = ["FID", "IID"]
    if args.include_intercept:
        eigenvec["Intercept"] = 1
        cols.append("Intercept")
    cols.extend([f"PC{i}" for i in range(1, args.num_pcs + 1)])
    eigenvec.loc[:, cols].to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
