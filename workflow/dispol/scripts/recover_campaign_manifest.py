#!/usr/bin/env python3.10
"""Rebuild a missing DISPOL campaign manifest from shard artifacts and launcher logs.

This helper is intended for recovery after a campaign crashes before
`campaigns/<tag>/manifest.json` is written. It reconstructs the expected shard
layout from the same campaign parameters used by `run_validation_campaign.py`,
marks shards with `.out` or `.yoda` artifacts as successful, marks launched
shards with only launcher/log artifacts as failed, and leaves missing shards
absent so the built-in `--rerun-failed-random-seed` flow can continue cleanly.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = SCRIPT_DIR.parent
DEFAULT_BASE_DIR = WORKFLOW_DIR
GAMMA_HELICITIES = ("00", "PP", "PM")
FULL_HELICITIES = ("00", "PP", "PM", "MP", "MM")
SETUP_ORDER = ("GAMMA", "Z", "ALL")
PS_SETUP_ORDER = ("SPINVAL", "SPINCOMP", "SPINHAD")
SUPPORTED_SETUPS = SETUP_ORDER + ("CC",) + PS_SETUP_ORDER
PS_SETUP_FAMILIES = {
    "SPINVAL": ("RIVETPS-SPIN",),
    "SPINCOMP": ("RIVETPS-SPIN", "RIVETPS-NOSPIN", "RIVETPS-NOSPIN-UNPOL"),
    "SPINHAD": ("RIVETPS-SPIN", "RIVETPS-NOSPIN", "RIVETPS-NOSPIN-UNPOL"),
}
SCALE_VARIATION_FACTORS = {
    "nominal": 1.0,
    "ScaleFactorDown": 0.5,
    "ScaleFactorUp": 2.0,
}
SCALE_VARIATION_ORDER = tuple(SCALE_VARIATION_FACTORS)
SCALE_VARIATION_STEM_SUFFIXES = {
    "nominal": "",
    "ScaleFactorDown": "-SCALEDOWN",
    "ScaleFactorUp": "-SCALEUP",
}


def normalize_campaign_setups(setups: Sequence[str]) -> List[str]:
    if not setups:
        return list(SETUP_ORDER)
    ordering = {name: index for index, name in enumerate(SUPPORTED_SETUPS)}
    seen = set()
    ordered: List[str] = []
    for setup in sorted((str(item).upper() for item in setups), key=lambda item: ordering.get(item, len(ordering))):
        if setup in seen:
            continue
        ordered.append(setup)
        seen.add(setup)
    return ordered


def is_ps_setup(setup: str) -> bool:
    return str(setup).upper() in PS_SETUP_ORDER


def contains_ps_setups(setups: Sequence[str]) -> bool:
    return any(is_ps_setup(setup) for setup in setups)


def contains_legacy_setups(setups: Sequence[str]) -> bool:
    return any(not is_ps_setup(setup) for setup in setups)


def validate_setup_family(setups: Sequence[str]) -> None:
    if contains_ps_setups(setups) and contains_legacy_setups(setups):
        raise ValueError("Do not mix SPINVAL/SPINCOMP/SPINHAD with GAMMA/Z/ALL/CC in the same invocation.")


def analysis_variant_from_args(args: argparse.Namespace) -> str:
    selected = int(bool(args.rivet)) + int(bool(args.rivetfo))
    if selected > 1:
        raise ValueError("Use at most one of --rivet and --rivetfo.")
    if args.rivetfo:
        return "RIVETFO"
    if args.rivet:
        return "RIVET"
    return ""


def use_lo_runtime(analysis_variant: str) -> bool:
    return analysis_variant not in {"RIVETFO", "RIVETPS"} and not analysis_variant.startswith("RIVETPS-")


def resolve_include_lo_requested(requested: Optional[bool], analysis_variant: str) -> bool:
    if requested is not None:
        return bool(requested)
    return use_lo_runtime(analysis_variant)


def selected_scale_variations(enabled: bool) -> Sequence[str]:
    return SCALE_VARIATION_ORDER if enabled else ("nominal",)


def scale_variation_stem_suffix(scale_variation: str) -> str:
    try:
        return SCALE_VARIATION_STEM_SUFFIXES[scale_variation]
    except KeyError as exc:
        raise ValueError(f"Unsupported scale variation {scale_variation!r}") from exc


def build_logical_run_stem(
    setup: str,
    order: str,
    helicity: str,
    analysis_variant: str = "",
    scale_variation: str = "nominal",
) -> str:
    analysis_suffix = f"-{analysis_variant}" if analysis_variant else ""
    if order == "LO":
        base = f"DIS-POL-LO_{helicity}-{setup}"
    else:
        base = f"DIS-POL-POWHEG_{helicity}-{order}-{setup}"
    return f"{base}{analysis_suffix}{scale_variation_stem_suffix(scale_variation)}"


def generated_scale_variation_in_file(tag: str, stem: str) -> str:
    return str(Path("campaigns") / tag / "generated-inputs" / f"{stem}.in")


def build_jobs(
    lo_events: int,
    posnlo_events: int,
    negnlo_events: int,
    analysis_variant: str = "",
    analysis_variants_by_setup: Optional[Dict[str, Sequence[str]]] = None,
    setups: Sequence[str] = (),
    scale_variations: bool = False,
    campaign_tag: str = "",
    include_lo: Optional[bool] = None,
) -> List[Dict[str, object]]:
    jobs: List[Dict[str, object]] = []
    include_lo = resolve_include_lo_requested(include_lo, analysis_variant)
    for setup in normalize_campaign_setups(setups):
        helicities = GAMMA_HELICITIES if setup == "GAMMA" else FULL_HELICITIES
        variants = (
            tuple(analysis_variants_by_setup.get(setup, (analysis_variant,)))
            if analysis_variants_by_setup is not None
            else (analysis_variant,)
        )
        if include_lo:
            for variant in variants:
                for hel in helicities:
                    nominal_stem = build_logical_run_stem(setup, "LO", hel, variant, "nominal")
                    for scale_variation in selected_scale_variations(scale_variations):
                        stem = build_logical_run_stem(setup, "LO", hel, variant, scale_variation)
                        jobs.append(
                            {
                                "setup": setup,
                                "order": "LO",
                                "helicity": hel,
                                "stem": stem,
                                "run_file": f"{stem}.run",
                                "in_file": (
                                    generated_scale_variation_in_file(campaign_tag, stem)
                                    if scale_variation != "nominal"
                                    else f"{stem}.in"
                                ),
                                "events": lo_events,
                                "analysis_variant": variant,
                                "scale_variation": scale_variation,
                                "scale_factor": SCALE_VARIATION_FACTORS[scale_variation],
                                "source_in_file": f"{nominal_stem}.in",
                            }
                        )
        for piece, events in (("POSNLO", posnlo_events), ("NEGNLO", negnlo_events)):
            for variant in variants:
                for hel in helicities:
                    nominal_stem = build_logical_run_stem(setup, piece, hel, variant, "nominal")
                    for scale_variation in selected_scale_variations(scale_variations):
                        stem = build_logical_run_stem(setup, piece, hel, variant, scale_variation)
                        jobs.append(
                            {
                                "setup": setup,
                                "order": piece,
                                "helicity": hel,
                                "stem": stem,
                                "run_file": f"{stem}.run",
                                "in_file": (
                                    generated_scale_variation_in_file(campaign_tag, stem)
                                    if scale_variation != "nominal"
                                    else f"{stem}.in"
                                ),
                                "events": events,
                                "analysis_variant": variant,
                                "scale_variation": scale_variation,
                                "scale_factor": SCALE_VARIATION_FACTORS[scale_variation],
                                "source_in_file": f"{nominal_stem}.in",
                            }
                        )
    return jobs


def choose_shards(logical_jobs: int, max_jobs: int, requested_shards: int) -> int:
    if requested_shards > 0:
        return requested_shards
    if logical_jobs <= 0:
        return 1
    return max(1, (max_jobs + logical_jobs - 1) // logical_jobs)


def split_events(total_events: int, shards: int) -> List[int]:
    if shards <= 1:
        return [total_events]
    base = total_events // shards
    remainder = total_events % shards
    parts = [base + (1 if index < remainder else 0) for index in range(shards)]
    return [part for part in parts if part > 0]


def build_shards(
    jobs: Sequence[Dict[str, object]],
    base_tag: str,
    requested_shards: int,
    max_jobs: int,
    seed_base: int,
) -> List[Dict[str, object]]:
    shard_target = choose_shards(len(jobs), max_jobs, requested_shards)
    shards: List[Dict[str, object]] = []
    seed = seed_base
    for job in jobs:
        event_parts = split_events(int(job["events"]), shard_target)
        shard_count = len(event_parts)
        for index, events in enumerate(event_parts, start=1):
            tag = base_tag if shard_count == 1 else f"{base_tag}-s{index:03d}"
            shards.append(
                {
                    "job": job,
                    "shard_index": index,
                    "shard_count": shard_count,
                    "tag": tag,
                    "seed": seed,
                    "events": events,
                    "rerun_parent_tag": "",
                }
            )
            seed += 1
    return shards


def build_artifact_index(base_dir: Path, expected_prefixes: Iterable[str]) -> Dict[str, Dict[str, List[str]]]:
    expected = set(expected_prefixes)
    index = {
        prefix: {"out_files": [], "log_files": [], "yoda_files": []}
        for prefix in expected
    }
    if not expected:
        return index

    for path in base_dir.iterdir():
        if not path.is_file():
            continue
        name = path.name
        prefix: Optional[str] = None
        bucket: Optional[str] = None
        if name.endswith(".out"):
            prefix = name[:-4]
            bucket = "out_files"
        elif name.endswith(".log"):
            prefix = name[:-4]
            bucket = "log_files"
        else:
            yoda_pos = name.find(".yoda")
            if yoda_pos != -1:
                prefix = name[:yoda_pos]
                bucket = "yoda_files"
        if prefix is None or bucket is None or prefix not in index:
            continue
        index[prefix][bucket].append(str(path))

    for artifact_lists in index.values():
        for values in artifact_lists.values():
            values.sort()
    return index


def build_launcher_log_index(launcher_dir: Path) -> set[str]:
    if not launcher_dir.exists():
        return set()
    return {path.name for path in launcher_dir.iterdir() if path.is_file()}


def resolve_analysis_configuration(args: argparse.Namespace) -> tuple[str, bool, bool, Optional[Dict[str, Sequence[str]]]]:
    requested_setups = normalize_campaign_setups(args.setup)
    validate_setup_family(requested_setups)
    analysis_variant = analysis_variant_from_args(args)
    include_lo = resolve_include_lo_requested(args.include_lo, analysis_variant)
    diagnostics_enabled = bool(args.diagnostics)
    analysis_variants_by_setup: Optional[Dict[str, Sequence[str]]] = None
    if contains_ps_setups(requested_setups):
        if analysis_variant == "RIVETFO":
            raise ValueError("SPINVAL/SPINCOMP/SPINHAD are not supported with --rivetfo.")
        if args.raw_powheg:
            raise ValueError("--raw-powheg is not supported for SPINVAL/SPINCOMP/SPINHAD.")
        if args.include_lo:
            raise ValueError("SPINVAL/SPINCOMP/SPINHAD are NLO-only workflows; do not use --include-lo.")
        analysis_variant = "RIVETPS"
        include_lo = False
        diagnostics_enabled = False
        analysis_variants_by_setup = {
            setup: PS_SETUP_FAMILIES[setup]
            for setup in requested_setups
            if is_ps_setup(setup)
        }
    return analysis_variant, include_lo, diagnostics_enabled, analysis_variants_by_setup


def manifest_finished_entry(
    spec: Dict[str, object],
    command: List[str],
    returncode: int,
    launcher_log: Optional[str],
    out_files: List[str],
    log_files: List[str],
    yoda_files: List[str],
) -> Dict[str, object]:
    return {
        "spec": spec,
        "command": command,
        "returncode": returncode,
        "started_at": None,
        "ended_at": None,
        "duration_s": None,
        "launcher_log": launcher_log,
        "out_files": out_files,
        "log_files": log_files,
        "yoda_files": yoda_files,
    }


def manifest_prepared_entry(job: Dict[str, object]) -> Dict[str, object]:
    return {
        "spec": {
            "job": job,
            "shard_index": 1,
            "shard_count": 1,
            "tag": "",
            "seed": 0,
            "events": int(job["events"]),
            "rerun_parent_tag": "",
        },
        "prepare_command": ["Herwig", "read", str(job["in_file"])],
        "prepare_returncode": 0,
    }


def build_manifest(args: argparse.Namespace) -> tuple[Dict[str, object], Dict[str, int]]:
    base_dir = args.base_dir.resolve()
    campaign_dir = base_dir / "campaigns" / args.tag
    analysis_variant, include_lo, diagnostics_enabled, analysis_variants_by_setup = resolve_analysis_configuration(args)
    requested_setups = normalize_campaign_setups(args.setup)

    jobs = build_jobs(
        args.lo_events,
        args.posnlo_events,
        args.negnlo_events,
        analysis_variant=analysis_variant,
        analysis_variants_by_setup=analysis_variants_by_setup,
        setups=requested_setups,
        scale_variations=bool(args.scale_variations),
        campaign_tag=args.tag,
        include_lo=include_lo,
    )
    shards = build_shards(jobs, args.tag, args.shards, max(1, args.jobs), args.seed_base)
    prepared = [
        manifest_prepared_entry(job)
        for job in jobs
        if (base_dir / str(job["run_file"])).exists()
    ]

    expected_prefixes = [f"{spec['job']['stem']}-S{spec['seed']}-{spec['tag']}" for spec in shards]
    artifact_index = build_artifact_index(base_dir, expected_prefixes)
    launcher_dir = campaign_dir / "launcher-logs"
    launcher_log_names = build_launcher_log_index(launcher_dir)

    finished: List[Dict[str, object]] = []
    failed: List[Dict[str, object]] = []
    finished_slots = 0
    failed_slots = 0

    for spec in shards:
        job = spec["job"]
        prefix = f"{job['stem']}-S{spec['seed']}-{spec['tag']}"
        launcher_name = f"{job['stem']}-{spec['tag']}.launcher.log"
        launcher_log = str(launcher_dir / launcher_name) if launcher_name in launcher_log_names else None
        artifacts = artifact_index.get(prefix, {"out_files": [], "log_files": [], "yoda_files": []})
        out_files = artifacts["out_files"]
        log_files = artifacts["log_files"]
        yoda_files = artifacts["yoda_files"]
        command = [
            "Herwig",
            "run",
            str(job["run_file"]),
            "-N",
            str(spec["events"]),
            "-t",
            str(spec["tag"]),
            "-s",
            str(spec["seed"]),
        ]
        has_success_artifacts = bool(out_files or yoda_files)
        has_failure_artifacts = bool(launcher_log or log_files)
        if has_success_artifacts:
            finished.append(
                manifest_finished_entry(spec, command, 0, launcher_log, out_files, log_files, yoda_files)
            )
            finished_slots += 1
        elif has_failure_artifacts:
            failed.append(
                manifest_finished_entry(spec, command, 1, launcher_log, out_files, log_files, yoda_files)
            )
            failed_slots += 1

    analysis_variants = sorted(
        {
            str(job["analysis_variant"])
            for job in jobs
            if str(job["analysis_variant"])
        }
    )
    manifest = {
        "tag": args.tag,
        "base_dir": str(base_dir),
        "jobs_limit": args.jobs,
        "shards_per_logical_run": args.shards,
        "seed_base": args.seed_base,
        "merge_yoda": not args.no_merge_yoda,
        "force_prepare": bool(args.force_prepare),
        "raw_powheg": bool(args.raw_powheg),
        "include_lo": include_lo,
        "extract_diagnostics": diagnostics_enabled,
        "scale_variations": bool(args.scale_variations),
        "scale_variation_factors": dict(SCALE_VARIATION_FACTORS),
        "yoda_merge_tool": args.yoda_merge_tool,
        "rivet": bool(analysis_variant),
        "analysis_variant": analysis_variant,
        "analysis_variants": analysis_variants,
        "ps_mode": bool(analysis_variants and all(item.startswith("RIVETPS-") for item in analysis_variants)),
        "campaign_setups": requested_setups,
        "event_counts": {
            "LO": args.lo_events,
            "POSNLO": args.posnlo_events,
            "NEGNLO": args.negnlo_events,
        },
        "prepared": prepared,
        "finished": finished,
        "failed": failed,
        "finished_count": len(finished),
        "failed_count": len(failed),
        "extract_status": {},
        "all_yoda_files": sorted({y for item in finished for y in item["yoda_files"]}),
        "merged_yoda": [],
        "nlo_yoda": [],
        "all_merged_yoda_files": [],
        "all_nlo_yoda_files": [],
    }
    summary = {
        "logical_jobs": len(jobs),
        "total_shards": len(shards),
        "prepared_jobs": len(prepared),
        "finished_shards": finished_slots,
        "failed_shards": failed_slots,
        "missing_shards": len(shards) - finished_slots - failed_slots,
    }
    return manifest, summary


def resume_command(args: argparse.Namespace) -> str:
    cmd = [
        "python3.10",
        str(SCRIPT_DIR / "run_validation_campaign.py"),
        "full",
        "--base-dir",
        str(args.base_dir.resolve()),
        "-t",
        args.tag,
    ]
    for setup in normalize_campaign_setups(args.setup):
        cmd.extend(["--setup", setup])
    cmd.extend(["--jobs", str(args.jobs)])
    cmd.extend(["--shards", str(args.shards)])
    cmd.extend(["--seed-base", str(args.seed_base)])
    cmd.extend(["--posnlo-events", str(args.posnlo_events)])
    cmd.extend(["--negnlo-events", str(args.negnlo_events)])
    if args.include_lo:
        cmd.append("--include-lo")
        cmd.extend(["--lo-events", str(args.lo_events)])
    elif not contains_ps_setups(args.setup) and not args.rivetfo:
        cmd.extend(["--lo-events", str(args.lo_events)])
    if args.rivet:
        cmd.append("--rivet")
    if args.rivetfo:
        cmd.append("--rivetfo")
    if args.scale_variations:
        cmd.append("--scale-variations")
    if args.raw_powheg:
        cmd.append("--raw-powheg")
    if args.diagnostics:
        cmd.append("--diagnostics")
    if args.keep_going:
        cmd.append("--keep-going")
    if args.force_prepare:
        cmd.append("--force-prepare")
    cmd.append("--rerun-failed-random-seed")
    if args.yoda_merge_tool:
        cmd.extend(["--yoda-merge-tool", args.yoda_merge_tool])
    return shlex.join(cmd)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Directory containing the DIS workflow cards, outputs, and campaigns/.",
    )
    parser.add_argument("-t", "--tag", required=True, help="Campaign tag to recover.")
    parser.add_argument(
        "--setup",
        action="append",
        choices=SUPPORTED_SETUPS,
        default=[],
        help="Setup(s) used for the campaign. May be given multiple times.",
    )
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1, help="Maximum concurrent Herwig jobs used for the campaign.")
    parser.add_argument("--shards", type=int, default=0, help="Shards per logical run used for the campaign.")
    parser.add_argument("--seed-base", type=int, default=100000, help="Base RNG seed used for the campaign.")
    parser.add_argument("--lo-events", type=int, default=1_000_000_000, help="Events per LO run.")
    parser.add_argument("--posnlo-events", type=int, default=1_000_000_000, help="Events per POSNLO run.")
    parser.add_argument("--negnlo-events", type=int, default=10_000_000, help="Events per NEGNLO run.")
    parser.add_argument("--rivet", action="store_true", help="Recover a --rivet campaign.")
    parser.add_argument("--rivetfo", action="store_true", help="Recover a --rivetfo campaign.")
    parser.add_argument("--raw-powheg", action="store_true", help="Campaign used --raw-powheg.")
    parser.add_argument("--include-lo", action="store_true", default=None, help="Campaign used --include-lo.")
    parser.add_argument("--diagnostics", action="store_true", default=None, help="Campaign used --diagnostics.")
    parser.add_argument("--scale-variations", action="store_true", default=None, help="Campaign used --scale-variations.")
    parser.add_argument("--force-prepare", action="store_true", help="Include --force-prepare in the printed resume command.")
    parser.add_argument("--no-merge-yoda", action="store_true", help="Campaign used --no-merge-yoda.")
    parser.add_argument("--keep-going", action="store_true", help="Include --keep-going in the printed resume command.")
    parser.add_argument("--yoda-merge-tool", default="auto", help="YODA merge tool to record and reuse in the printed resume command.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite an existing manifest.json.")
    parser.add_argument("--dry-run", action="store_true", help="Print the recovered summary without writing manifest.json.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    campaign_dir = args.base_dir.resolve() / "campaigns" / args.tag
    manifest_path = campaign_dir / "manifest.json"
    if manifest_path.exists() and not (args.overwrite or args.dry_run):
        raise FileExistsError(f"{manifest_path} already exists. Use --overwrite to replace it.")

    manifest, summary = build_manifest(args)
    if args.dry_run:
        print(f"[dry-run] would write {manifest_path}")
    else:
        campaign_dir.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
        print(f"Wrote {manifest_path}")

    print(
        "Recovered campaign summary: "
        f"logical_jobs={summary['logical_jobs']} "
        f"total_shards={summary['total_shards']} "
        f"prepared_jobs={summary['prepared_jobs']} "
        f"finished_shards={summary['finished_shards']} "
        f"failed_shards={summary['failed_shards']} "
        f"missing_shards={summary['missing_shards']}"
    )
    print()
    print("Resume command:")
    print(resume_command(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
