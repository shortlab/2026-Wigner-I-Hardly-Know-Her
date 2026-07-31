from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io_parser import parser_smoke
from .workflows.indium_analysis import run_indium_analysis
from .workflows.interface_experiment import run_interface_experiment


def _build_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser and subcommands for all pipeline workflows."""
    p = argparse.ArgumentParser(description="Wigner MATLAB-to-Python pipeline CLI")
    sub = p.add_subparsers(dest="cmd", required=True)

    s1 = sub.add_parser("parser-smoke", help="Parse a text file and print segment summary")
    s1.add_argument("--file", required=True, type=Path)

    s2 = sub.add_parser("run-interface-experiment", help="Run baseline subtraction plotting workflow")
    s2.add_argument("--workspace", default=Path.cwd(), type=Path)
    s2.add_argument("--baseline-pattern", default="cle_ind_baseline")
    s2.add_argument("--output-dir", default="outputs/interface_experiment")

    s3 = sub.add_parser("run-indium-analysis", help="Run parameterized indium analysis workflow")
    s3.add_argument("--file", required=True, type=Path)
    s3.add_argument("--output-dir", default="outputs/indium_analysis")
    s3.add_argument("--n-mc", type=int, default=1000)
    s3.add_argument("--seed", type=int, default=42)
    s3.add_argument("--fixed-lr-center", type=float, default=None)
    s3.add_argument("--fixed-rl-center", type=float, default=None)
    s3.add_argument("--fixed-outer-pad-c", type=float, default=15.0)
    s3.add_argument("--baseline-fit-offset-mw", type=float, default=1.0)
    return p


def main(argv: list[str] | None = None) -> int:
    """Dispatch CLI command, execute workflow, and print JSON result payload."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "parser-smoke":
        res = parser_smoke(args.file)
    elif args.cmd == "run-interface-experiment":
        res = run_interface_experiment(args.workspace, baseline_pattern=args.baseline_pattern, output_dir=args.output_dir)
    elif args.cmd == "run-indium-analysis":
        res = run_indium_analysis(
            args.file,
            output_dir=args.output_dir,
            n_mc=args.n_mc,
            seed=args.seed,
            fixed_lr_center=args.fixed_lr_center,
            fixed_rl_center=args.fixed_rl_center,
            fixed_outer_pad_c=args.fixed_outer_pad_c,
            baseline_fit_offset_mw=args.baseline_fit_offset_mw,
        )
    else:
        parser.error("Unknown command")
        return 2

    print(json.dumps(res, indent=2))
    return 0
