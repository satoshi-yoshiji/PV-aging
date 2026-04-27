"""Top-level entry point: runs Section 1 then Section 2."""
from __future__ import annotations
import argparse, logging, subprocess, sys
from pathlib import Path

HERE = Path(__file__).parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="/humgen/diabetes2/users/satoshi/misc/01.pv/curated_stats.tsv.gz")
    parser.add_argument("--outdir", default=str(HERE / "results"))
    parser.add_argument("--ukb-extract", default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-section1", action="store_true")
    parser.add_argument("--skip-section2", action="store_true")
    args = parser.parse_args()

    py = sys.executable
    common = ["--input", args.input, "--outdir", args.outdir, "--seed", str(args.seed)]
    if not args.skip_section1:
        cmd = [py, str(HERE / "run_section1_baseline_aging.py")] + common
        if args.ukb_extract:
            cmd += ["--ukb-extract", args.ukb_extract]
        print(">>", " ".join(cmd))
        subprocess.check_call(cmd)
    if not args.skip_section2:
        cmd = [py, str(HERE / "run_section2_aging_prospective.py")] + common
        print(">>", " ".join(cmd))
        subprocess.check_call(cmd)

if __name__ == "__main__":
    main()
