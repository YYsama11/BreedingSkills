#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run EMMAX association scans in parallel.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--tped-prefix", required=True)
    parser.add_argument("--kinship", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--parallel", type=int, default=10)
    parser.add_argument("--workspace", required=True)
    return parser.parse_args()


def run_one(row: pd.Series, tped_prefix: str, kinship: str, covariates: str, workspace: str) -> tuple[str, int]:
    trait_id = row["trait_id"]
    pheno_file = row["phenotype_file"]
    script_dir = Path(__file__).resolve().parent
    wrapper = script_dir / "run_with_top.sh"
    runner = script_dir / "run_emmax_trait.sh"
    cmd = [
        str(wrapper),
        f"emmax_{trait_id}",
        str(runner),
        trait_id,
        pheno_file,
        tped_prefix,
        kinship,
        covariates,
    ]
    proc = subprocess.run(cmd, cwd=workspace)
    return trait_id, proc.returncode


def main() -> None:
    args = parse_args()
    workspace = Path(args.workspace)
    log_dir = workspace / "analysis" / "logs" / "emmax"
    log_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest, sep="\t")
    status_records: list[dict[str, object]] = []

    with ThreadPoolExecutor(max_workers=args.parallel) as executor:
        futures = {
            executor.submit(
                run_one,
                row,
                args.tped_prefix,
                args.kinship,
                args.covariates,
                str(workspace),
            ): row["trait_id"]
            for _, row in manifest.iterrows()
        }
        for future in as_completed(futures):
            trait_id = futures[future]
            try:
                trait, code = future.result()
            except Exception as exc:
                status_records.append({"trait_id": trait_id, "status": "failed", "return_code": -1, "message": str(exc)})
            else:
                status_records.append(
                    {
                        "trait_id": trait,
                        "status": "ok" if code == 0 else "failed",
                        "return_code": code,
                        "message": "",
                    }
                )

    status_df = pd.DataFrame(status_records).sort_values(["status", "trait_id"])
    status_df.to_csv(log_dir / "emmax_batch_status.tsv", sep="\t", index=False)

    failed = status_df[status_df["status"] != "ok"]
    if not failed.empty:
        raise SystemExit(f"EMMAX failed for {failed.shape[0]} traits")


if __name__ == "__main__":
    main()
