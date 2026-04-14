#!/usr/bin/env python3.10
"""
Run the smallest GAMMA-focused follow-up for the residual polarized-NLO deficit.

By default this wrapper launches two diagnostics under one tag:
  * a broad-window polarized-NLO branch study with the paired NNPDF setup
  * a GAMMA PDF/extractor object audit with the same paired NNPDF setup

The goal is not to replace the full validation workflow. It is to answer the
next discriminator question with one command and one shared campaign tag.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence


DEFAULT_BASE_DIR = Path(__file__).resolve().parent
DEFAULT_TAG = "plain56nn-gamma-residual-discriminator"
DEFAULT_PDF_PROFILE = "nnpdf_paired"
DEFAULT_BRANCH_JOBS = 4
DEFAULT_BRANCH_SHARDS = 4
DEFAULT_BRANCH_EVENTS_PER_SHARD = 2_000_000
DEFAULT_BRANCH_SEED_BASE = 740_000
DEFAULT_AUDIT_JOBS = 2
DEFAULT_AUDIT_SHARDS = 1
DEFAULT_AUDIT_EVENTS_PER_SHARD = 1_000_000
DEFAULT_AUDIT_SEED_BASE = 750_000


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        nargs="?",
        choices=("all", "branch", "pdf-audit"),
        default="all",
        help="Which discriminator(s) to run.",
    )
    parser.add_argument("--tag", default=DEFAULT_TAG, help="Shared campaign tag for both diagnostics.")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Directory containing the DISPOL helpers and card templates.",
    )
    parser.add_argument(
        "--pdf-profile",
        choices=("hybrid", "nnpdf_paired"),
        default=DEFAULT_PDF_PROFILE,
        help="PDF profile to materialize into both diagnostics.",
    )
    parser.add_argument(
        "--branch-windows",
        default="broad",
        help="Comma-separated raw-delta windows to run. Default keeps the follow-up minimal.",
    )
    parser.add_argument(
        "--branch-modes",
        default="eff,raw",
        help="Comma-separated polarized-NLO modes for the branch study (legacy,eff,raw).",
    )
    parser.add_argument(
        "--with-poldis",
        action="store_true",
        help="Also build missing POLDIS references for the raw-delta study.",
    )
    parser.add_argument(
        "--poldis-totals-json",
        type=Path,
        help="Optional broad-window GAMMA override from an existing campaign-level poldis-totals.json.",
    )
    parser.add_argument(
        "--skip-herwig",
        action="store_true",
        help="Reuse existing Herwig outputs. For the PDF audit this becomes collect-only.",
    )
    parser.add_argument(
        "--collect-only",
        action="store_true",
        help="Skip all execution and rebuild summaries from existing outputs only.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print the subcommands without executing them.")

    parser.add_argument("--branch-jobs", type=int, default=DEFAULT_BRANCH_JOBS, help="Concurrent Herwig jobs for the raw-delta study.")
    parser.add_argument("--branch-shards", type=int, default=DEFAULT_BRANCH_SHARDS, help="Shards per logical run for the raw-delta study.")
    parser.add_argument(
        "--branch-events-per-shard",
        type=int,
        default=DEFAULT_BRANCH_EVENTS_PER_SHARD,
        help="Events per raw-delta shard.",
    )
    parser.add_argument("--branch-seed-base", type=int, default=DEFAULT_BRANCH_SEED_BASE, help="Base RNG seed for the raw-delta study.")
    parser.add_argument("--poldis-jobs", type=int, default=2, help="Concurrent POLDIS helper jobs when --with-poldis is enabled.")
    parser.add_argument(
        "--poldis-variant-jobs",
        type=int,
        default=1,
        help="Concurrent polarized/unpolarized jobs inside each POLDIS helper.",
    )
    parser.add_argument("--poldis-events", type=int, default=50_000_000, help="Events per POLDIS reference run when enabled.")

    parser.add_argument("--audit-jobs", type=int, default=DEFAULT_AUDIT_JOBS, help="Concurrent Herwig jobs for the PDF audit.")
    parser.add_argument("--audit-shards", type=int, default=DEFAULT_AUDIT_SHARDS, help="Shards per logical run for the PDF audit.")
    parser.add_argument(
        "--audit-events-per-shard",
        type=int,
        default=DEFAULT_AUDIT_EVENTS_PER_SHARD,
        help="Events per PDF-audit shard.",
    )
    parser.add_argument("--audit-seed-base", type=int, default=DEFAULT_AUDIT_SEED_BASE, help="Base RNG seed for the PDF audit.")
    return parser


def build_branch_command(args: argparse.Namespace) -> List[str]:
    script = args.base_dir / "run_gamma_rawdelta_branch_study.py"
    cmd = [
        sys.executable,
        str(script),
        "--tag",
        args.tag,
        "--base-dir",
        str(args.base_dir),
        "--pdf-profiles",
        args.pdf_profile,
        "--windows",
        args.branch_windows,
        "--branch-modes",
        args.branch_modes,
        "--jobs",
        str(args.branch_jobs),
        "--shards",
        str(args.branch_shards),
        "--events-per-shard",
        str(args.branch_events_per_shard),
        "--seed-base",
        str(args.branch_seed_base),
    ]
    if args.with_poldis:
        cmd.extend(
            [
                "--poldis-jobs",
                str(args.poldis_jobs),
                "--poldis-variant-jobs",
                str(args.poldis_variant_jobs),
                "--poldis-events",
                str(args.poldis_events),
            ]
        )
    else:
        cmd.append("--skip-poldis")
    if args.poldis_totals_json is not None:
        cmd.extend(["--poldis-totals-json", str(args.poldis_totals_json)])
    if args.skip_herwig:
        cmd.append("--skip-herwig")
    if args.collect_only:
        cmd.append("--collect-only")
    if args.dry_run:
        cmd.append("--dry-run")
    return cmd


def build_audit_command(args: argparse.Namespace) -> List[str]:
    script = args.base_dir / "run_gamma_pdf_object_audit.py"
    cmd = [
        sys.executable,
        str(script),
        "--tag",
        args.tag,
        "--base-dir",
        str(args.base_dir),
        "--pdf-profile",
        args.pdf_profile,
        "--jobs",
        str(args.audit_jobs),
        "--shards",
        str(args.audit_shards),
        "--events-per-shard",
        str(args.audit_events_per_shard),
        "--seed-base",
        str(args.audit_seed_base),
    ]
    if args.collect_only or args.skip_herwig:
        cmd.append("--collect-only")
    if args.dry_run:
        cmd.append("--dry-run")
    return cmd


def selected_commands(args: argparse.Namespace) -> List[List[str]]:
    commands: List[List[str]] = []
    if args.mode in {"all", "branch"}:
        commands.append(build_branch_command(args))
    if args.mode in {"all", "pdf-audit"}:
        commands.append(build_audit_command(args))
    return commands


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    args.base_dir = args.base_dir.resolve()
    if not args.base_dir.exists():
        raise SystemExit(f"Base directory does not exist: {args.base_dir}")

    commands = selected_commands(args)
    if args.dry_run:
        for index, cmd in enumerate(commands, start=1):
            print(f"[plan:gamma-residual-discriminator] step={index}/{len(commands)}")
            print(shlex.join(cmd))
        return 0

    for index, cmd in enumerate(commands, start=1):
        print(f"[run:gamma-residual-discriminator] step={index}/{len(commands)}: {shlex.join(cmd)}", flush=True)
        proc = subprocess.run(cmd, cwd=args.base_dir, text=True)
        if proc.returncode != 0:
            return proc.returncode

    print("")
    print("Expected outputs:")
    if args.mode in {"all", "branch"}:
        print(f"  raw-delta study: {args.base_dir / 'campaigns' / args.tag / 'gamma-rawdelta-branch-study' / 'results.txt'}")
    if args.mode in {"all", "pdf-audit"}:
        print(f"  pdf audit:       {args.base_dir / 'campaigns' / args.tag / 'pdf-audit' / 'results.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
