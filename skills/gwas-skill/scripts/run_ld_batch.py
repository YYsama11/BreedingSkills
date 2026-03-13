#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run PLINK LD jobs for lead SNPs in parallel.")
    parser.add_argument("--lead-file", required=True)
    parser.add_argument("--bfile-prefix", required=True)
    parser.add_argument("--workspace", required=True)
    parser.add_argument("--parallel", type=int, default=4)
    return parser.parse_args()


def run_one(lead_snp: str, bfile_prefix: str, workspace: str) -> tuple[str, int]:
    workspace_path = Path(workspace)
    wrapper = workspace_path / "scripts" / "run_with_top.sh"
    runner = workspace_path / "scripts" / "run_ld_for_lead.sh"
    outdir = workspace_path / "analysis" / "qtl" / "ld"
    outdir.mkdir(parents=True, exist_ok=True)
    safe_name = lead_snp.replace(":", "_")
    out_prefix = outdir / safe_name
    cmd = [
        str(wrapper),
        f"ld_{safe_name}",
        str(runner),
        lead_snp,
        bfile_prefix,
        str(out_prefix),
    ]
    proc = subprocess.run(cmd, cwd=workspace)
    return lead_snp, proc.returncode


def main() -> None:
    args = parse_args()
    workspace = Path(args.workspace)
    leads = [line.strip() for line in Path(args.lead_file).read_text().splitlines() if line.strip()]
    status = []
    with ThreadPoolExecutor(max_workers=args.parallel) as executor:
        futures = {executor.submit(run_one, lead, args.bfile_prefix, str(workspace)): lead for lead in leads}
        for future in as_completed(futures):
            lead = futures[future]
            try:
                _, code = future.result()
            except Exception as exc:
                status.append({"lead_snp": lead, "status": "failed", "return_code": -1, "message": str(exc)})
            else:
                status.append({"lead_snp": lead, "status": "ok" if code == 0 else "failed", "return_code": code, "message": ""})

    status_df = pd.DataFrame(status).sort_values(["status", "lead_snp"])
    status_df.to_csv(workspace / "analysis" / "qtl" / "ld_batch_status.tsv", sep="\t", index=False)
    if (status_df["status"] != "ok").any():
        raise SystemExit("Some LD jobs failed")


if __name__ == "__main__":
    main()
