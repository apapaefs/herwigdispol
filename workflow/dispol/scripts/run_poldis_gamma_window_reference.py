#!/usr/bin/env python3.10
"""
Build focused POLDIS GAMMA total-cross-section references for a chosen DIS window.

The helper:
  * templates user_dijet_rivetplots.f into a campaign-local work area
  * compiles a local poldis.x in separate unpolarized/polarized subdirectories
  * runs the two jobs
  * parses the printed LO/NLO/NNLO inclusive totals
  * writes reference.txt and reference.json

The default target is the interior GAMMA window used by the focused Herwig
follow-up study:
  Q^2 in [100,1000] GeV^2, y in [0.3,0.5]
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import re
import shlex
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence

from extract_dis_out_results import Measurement, POLDIS_POL_REFS, POLDIS_UNPOL_REFS


DEFAULT_BASE_DIR = Path(__file__).resolve().parent
DEFAULT_ROOT_DIR = DEFAULT_BASE_DIR.parent.parent
DEFAULT_POLDIS_DIR = DEFAULT_ROOT_DIR / "POLDIS" / "POLDIS-public"
DEFAULT_TAG = "plain50-gamma-poldis-ref"
DEFAULT_EVENTS = 200_000_000
DEFAULT_PDF_PROFILE = "hybrid"
DEFAULT_JOBS = 1

PDF_PROFILES: Dict[str, Dict[str, str]] = {
    "hybrid": {
        "unpolarized": "PDF4LHC15_nnlo_100_pdfas",
        "polarized_diff": "BDSSV24-NNLO",
    },
    "nnpdf_paired": {
        "unpolarized": "NNPDF40_nlo_pch_as_01180",
        "polarized_diff": "NNPDFpol20_nlo_as_01180",
    },
}


@dataclass(frozen=True)
class CutWindow:
    label: str
    q2_min: float
    q2_max: float
    y_min: float
    y_max: float


BROAD_WINDOW = CutWindow(label="broad", q2_min=49.0, q2_max=2500.0, y_min=0.2, y_max=0.6)
INTERIOR_WINDOW = CutWindow(label="interior", q2_min=100.0, q2_max=1000.0, y_min=0.3, y_max=0.5)
WINDOWS = {"broad": BROAD_WINDOW, "interior": INTERIOR_WINDOW}

TOTAL_RE = re.compile(
    r"^(?P<order>NNLO|NLO|LO)\s*=\s*(?P<value>[+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*\+\-\s*"
    r"(?P<error>[+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*$"
)
PROGRESS_RE = re.compile(r"^\s*(?P<events>\d+),\s*ISEED=")

RUNTIME_FILES = ("poldis.f", "gbook.f", "jetalg.f")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default=DEFAULT_TAG, help="Campaign tag used for the reference work area.")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Directory containing the DISPOL helpers.",
    )
    parser.add_argument(
        "--poldis-dir",
        type=Path,
        default=DEFAULT_POLDIS_DIR,
        help="Directory containing the POLDIS source tree and compile helper.",
    )
    parser.add_argument(
        "--window",
        choices=tuple(WINDOWS),
        default="interior",
        help="DIS window to evaluate.",
    )
    parser.add_argument(
        "--pdf-profile",
        choices=tuple(PDF_PROFILES),
        default=DEFAULT_PDF_PROFILE,
        help="Named PDF profile used for the unpolarized and polarized POLDIS runs.",
    )
    parser.add_argument(
        "--unpolarized-pdf",
        help="Optional override for the unpolarized POLDIS PDF set name.",
    )
    parser.add_argument(
        "--polarized-diff-pdf",
        help="Optional override for the polarized-difference POLDIS PDF set name.",
    )
    parser.add_argument("--events", type=int, default=DEFAULT_EVENTS, help="Number of POLDIS events per run.")
    parser.add_argument(
        "--jobs",
        type=int,
        default=DEFAULT_JOBS,
        help="Maximum concurrent unpolarized/polarized POLDIS jobs inside this reference build.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print planned compile/run commands only.")
    parser.add_argument(
        "--collect-only",
        action="store_true",
        help="Skip compile/run and rebuild the parsed reference from existing logs.",
    )
    return parser


def work_dir(base_dir: Path, tag: str, window: CutWindow) -> Path:
    suffix = "poldis-gamma-interior" if window.label == "interior" else "poldis-gamma-broad"
    return base_dir / "campaigns" / tag / suffix


def variant_dir(base_dir: Path, tag: str, window: CutWindow, label: str) -> Path:
    return work_dir(base_dir, tag, window) / label


def reference_txt_path(base_dir: Path, tag: str, window: CutWindow) -> Path:
    return work_dir(base_dir, tag, window) / "reference.txt"


def reference_json_path(base_dir: Path, tag: str, window: CutWindow) -> Path:
    return work_dir(base_dir, tag, window) / "reference.json"


def monitor_dir(base_dir: Path, tag: str, window: CutWindow) -> Path:
    return work_dir(base_dir, tag, window) / "monitor"


def status_txt_path(base_dir: Path, tag: str, window: CutWindow) -> Path:
    return monitor_dir(base_dir, tag, window) / "status.txt"


def status_json_path(base_dir: Path, tag: str, window: CutWindow) -> Path:
    return monitor_dir(base_dir, tag, window) / "status.json"


def measurement_payload(measurement: Measurement) -> dict:
    return {"value_pb": measurement.value_pb, "error_pb": measurement.error_pb}


def measurement_from_payload(payload: Mapping[str, float]) -> Measurement:
    return Measurement(float(payload["value_pb"]), float(payload["error_pb"]))


def fmt_measurement(payload: Mapping[str, float]) -> str:
    return f"{payload['value_pb']:.10f} +- {payload['error_pb']:.10f} pb"


def fmt_seconds(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "n/a"
    total = int(round(value))
    hours, rem = divmod(total, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{hours:d}:{minutes:02d}:{seconds:02d}"


def sub_measurement(a: Measurement, b: Measurement) -> Measurement:
    return Measurement(a.value_pb - b.value_pb, math.sqrt(a.error_pb ** 2 + b.error_pb ** 2))


def parse_poldis_totals(text: str, context: str) -> Dict[str, Measurement]:
    totals: Dict[str, Measurement] = {}
    for line in text.splitlines():
        match = TOTAL_RE.match(line.strip())
        if not match:
            continue
        totals[match.group("order")] = Measurement(
            value_pb=float(match.group("value")),
            error_pb=float(match.group("error")),
        )
    missing = [order for order in ("LO", "NLO", "NNLO") if order not in totals]
    if missing:
        raise RuntimeError(f"Failed to parse {context}: missing {', '.join(missing)} totals")
    return totals


def compile_command(poldis_dir: Path) -> list[str]:
    compile_script = poldis_dir / "compile_dijet_rivetplots"
    if not compile_script.exists():
        raise FileNotFoundError(f"Missing compile helper {compile_script}")
    command = compile_script.read_text().strip()
    if not command:
        raise RuntimeError(f"Compile helper {compile_script} is empty")
    return shlex.split(command)


def resolve_pdf_choice(profile: str, unpolarized_override: str | None, polarized_override: str | None) -> Dict[str, str]:
    if profile not in PDF_PROFILES:
        raise KeyError(f"Unknown PDF profile {profile!r}")
    resolved = dict(PDF_PROFILES[profile])
    if unpolarized_override:
        resolved["unpolarized"] = unpolarized_override
    if polarized_override:
        resolved["polarized_diff"] = polarized_override
    return resolved


def is_builtin_hybrid_pdf_choice(pdfs: Mapping[str, str]) -> bool:
    return (
        pdfs.get("unpolarized") == PDF_PROFILES["hybrid"]["unpolarized"]
        and pdfs.get("polarized_diff") == PDF_PROFILES["hybrid"]["polarized_diff"]
    )


def patch_user_source(
    source_text: str,
    ipol: int,
    window: CutWindow,
    events: int,
    unpolarized_pdf: str,
    polarized_diff_pdf: str,
) -> str:
    replacements = {
        "IPOL=": f"      IPOL={ipol}",
        "IBOSON=": "      IBOSON=0",
        "LEPCH=": "      LEPCH=-1D0",
        "IFRAME=": "      IFRAME=1",
        "ICH=": "      ICH=0",
        "INNLO=": "      INNLO=1",
        "NEV=": f"      NEV={events}",
        "YMIN=": f"      YMIN={window.y_min}",
        "YMAX=": f"      YMAX={window.y_max}",
        "Q2MIN=": f"      Q2MIN={window.q2_min:.0f}",
        "Q2MAX=": f"      Q2MAX={window.q2_max:.0f}",
    }
    seen = {key: False for key in replacements}
    seen_unpolarized_pdf = False
    seen_polarized_pdf = False
    patched = []
    pdf_branch: str | None = None
    for line in source_text.splitlines():
        stripped = line.strip()
        if stripped.startswith("IF (IPOL.EQ.0) THEN"):
            pdf_branch = "unpolarized"
            patched.append(line)
            continue
        if stripped.startswith("ELSEIF (IPOL.EQ.1) THEN"):
            pdf_branch = "polarized_diff"
            patched.append(line)
            continue
        if stripped.startswith("ENDIF"):
            pdf_branch = None
            patched.append(line)
            continue
        if stripped.startswith("call InitPDFsetByName("):
            if pdf_branch == "unpolarized":
                patched.append(f'         call InitPDFsetByName("{unpolarized_pdf}")')
                seen_unpolarized_pdf = True
                continue
            if pdf_branch == "polarized_diff":
                patched.append(f'         call InitPDFsetByName("{polarized_diff_pdf}")')
                seen_polarized_pdf = True
                continue
        replaced = False
        for prefix, replacement in replacements.items():
            if stripped.startswith(prefix):
                patched.append(replacement)
                seen[prefix] = True
                replaced = True
                break
        if not replaced:
            patched.append(line)
    missing = [prefix for prefix, ok in seen.items() if not ok]
    if missing:
        raise RuntimeError(f"Failed to patch POLDIS source lines: {missing}")
    if not seen_unpolarized_pdf or not seen_polarized_pdf:
        raise RuntimeError("Failed to patch POLDIS PDF set lines in user_dijet_rivetplots.f")
    rendered = "\n".join(patched)
    if source_text.endswith("\n"):
        rendered += "\n"
    return rendered


def expected_top_name(ipol: int) -> str:
    return "dijets_unpol_GAM.top" if ipol == 0 else "dijets_pol_GAM.top"


def render_monitor_text(payload: Mapping[str, object]) -> str:
    lines = [
        f"Tag: {payload['tag']}",
        f"Window: {payload['window']}",
        f"Phase: {payload['phase']}",
        f"Variant: {payload['variant']}",
        f"Elapsed: {fmt_seconds(float(payload['elapsed_s']))}",
    ]
    if payload.get("events_total") is not None:
        lines.extend(
            [
                f"Events done: {payload['events_done']}/{payload['events_total']}",
                f"Variant progress: {float(payload['variant_percent']):.2f}%",
                f"Overall progress: {float(payload['overall_percent']):.2f}%",
            ]
        )
    if payload.get("last_progress_line"):
        lines.append(f"Last progress line: {payload['last_progress_line']}")
    if payload.get("log"):
        lines.append(f"Log: {payload['log']}")
    return "\n".join(lines).rstrip() + "\n"


def write_monitor_files(base_dir: Path, tag: str, window: CutWindow, payload: Mapping[str, object]) -> None:
    mon_dir = monitor_dir(base_dir, tag, window)
    mon_dir.mkdir(parents=True, exist_ok=True)
    status_json_path(base_dir, tag, window).write_text(json.dumps(payload, indent=2, sort_keys=True))
    status_txt_path(base_dir, tag, window).write_text(render_monitor_text(payload))


def build_monitor_payload(
    tag: str,
    window: CutWindow,
    phase: str,
    variant: str,
    started_at: float,
    run_index: int,
    run_count: int,
    events_total: int | None = None,
    events_done: int = 0,
    last_progress_line: str | None = None,
    log_path: Path | None = None,
) -> dict:
    variant_fraction = 0.0
    if events_total and events_total > 0:
        variant_fraction = min(max(events_done / events_total, 0.0), 1.0)
    overall_fraction = (run_index + variant_fraction) / max(run_count, 1)
    return {
        "tag": tag,
        "window": window.label,
        "phase": phase,
        "variant": variant,
        "elapsed_s": time.time() - started_at,
        "events_total": events_total,
        "events_done": events_done,
        "variant_percent": 100.0 * variant_fraction,
        "overall_percent": 100.0 * overall_fraction,
        "last_progress_line": last_progress_line,
        "log": str(log_path) if log_path is not None else None,
    }


def materialize_variant_sources(
    poldis_dir: Path,
    target_dir: Path,
    ipol: int,
    window: CutWindow,
    events: int,
    unpolarized_pdf: str,
    polarized_diff_pdf: str,
) -> None:
    target_dir.mkdir(parents=True, exist_ok=True)
    for filename in RUNTIME_FILES:
        source_path = poldis_dir / filename
        if not source_path.exists():
            raise FileNotFoundError(f"Missing POLDIS runtime source {source_path}")
        shutil.copy2(source_path, target_dir / filename)

    template_path = poldis_dir / "user_dijet_rivetplots.f"
    if not template_path.exists():
        raise FileNotFoundError(f"Missing POLDIS user template {template_path}")
    rendered = patch_user_source(
        template_path.read_text(),
        ipol=ipol,
        window=window,
        events=events,
        unpolarized_pdf=unpolarized_pdf,
        polarized_diff_pdf=polarized_diff_pdf,
    )
    rendered_path = target_dir / "user_dijet_rivetplots.f"
    if not rendered_path.exists() or rendered_path.read_text() != rendered:
        rendered_path.write_text(rendered)


def run_logged(cmd: Sequence[str], cwd: Path, log_path: Path, dry_run: bool) -> None:
    if dry_run:
        print(shlex.join(list(cmd)))
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as handle:
        handle.write(f"$ {' '.join(cmd)}\n")
        handle.flush()
        proc = subprocess.run(cmd, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed in {cwd}: {shlex.join(list(cmd))}\nSee {log_path}")


def run_poldis_with_progress(
    base_dir: Path,
    tag: str,
    window: CutWindow,
    run_dir: Path,
    variant_label: str,
    events: int,
    run_index: int,
    run_count: int,
    started_at: float,
    dry_run: bool,
) -> None:
    cmd = ["./poldis.x"]
    log_path = run_dir / "run.log"
    if dry_run:
        print(shlex.join(cmd))
        return

    log_path.parent.mkdir(parents=True, exist_ok=True)
    payload = build_monitor_payload(
        tag=tag,
        window=window,
        phase=f"running-{variant_label}",
        variant=variant_label,
        started_at=started_at,
        run_index=run_index,
        run_count=run_count,
        events_total=events,
        events_done=0,
        log_path=log_path,
    )
    write_monitor_files(base_dir, tag, window, payload)

    with log_path.open("w") as handle:
        handle.write(f"$ {' '.join(cmd)}\n")
        handle.flush()
        proc = subprocess.Popen(
            cmd,
            cwd=run_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        last_events = 0
        last_line = None
        for line in proc.stdout:
            handle.write(line)
            handle.flush()
            stripped = line.rstrip()
            match = PROGRESS_RE.match(stripped)
            if match:
                last_events = int(match.group("events"))
                last_line = stripped
                payload = build_monitor_payload(
                    tag=tag,
                    window=window,
                    phase=f"running-{variant_label}",
                    variant=variant_label,
                    started_at=started_at,
                    run_index=run_index,
                    run_count=run_count,
                    events_total=events,
                    events_done=last_events,
                    last_progress_line=last_line,
                    log_path=log_path,
                )
                write_monitor_files(base_dir, tag, window, payload)
        returncode = proc.wait()

    if returncode != 0:
        raise RuntimeError(f"Command failed in {run_dir}: {shlex.join(cmd)}\nSee {log_path}")

    payload = build_monitor_payload(
        tag=tag,
        window=window,
        phase=f"completed-{variant_label}",
        variant=variant_label,
        started_at=started_at,
        run_index=run_index + 1,
        run_count=run_count,
        events_total=events,
        events_done=events,
        last_progress_line=last_line,
        log_path=log_path,
    )
    write_monitor_files(base_dir, tag, window, payload)


def ensure_compiled(poldis_dir: Path, run_dir: Path, dry_run: bool) -> None:
    cmd = compile_command(poldis_dir)
    binary = run_dir / "poldis.x"
    if not dry_run and binary.exists():
        source_mtime = max((run_dir / filename).stat().st_mtime for filename in (*RUNTIME_FILES, "user_dijet_rivetplots.f"))
        if binary.stat().st_mtime >= source_mtime:
            return
    run_logged(cmd, cwd=run_dir, log_path=run_dir / "compile.log", dry_run=dry_run)


def ensure_ran(run_dir: Path, dry_run: bool) -> None:
    run_logged(["./poldis.x"], cwd=run_dir, log_path=run_dir / "run.log", dry_run=dry_run)


def process_variant(
    *,
    base_dir: Path,
    tag: str,
    window: CutWindow,
    poldis_dir: Path,
    label: str,
    ipol: int,
    events: int,
    pdfs: Mapping[str, str],
    started_at: float,
    dry_run: bool,
    parallel_mode: bool,
) -> str:
    run_dir = variant_dir(base_dir, tag, window, label)
    materialize_variant_sources(
        poldis_dir,
        run_dir,
        ipol=ipol,
        window=window,
        events=events,
        unpolarized_pdf=pdfs["unpolarized"],
        polarized_diff_pdf=pdfs["polarized_diff"],
    )
    if not dry_run and not parallel_mode:
        write_monitor_files(
            base_dir,
            tag,
            window,
            build_monitor_payload(
                tag=tag,
                window=window,
                phase=f"compiling-{label}",
                variant=label,
                started_at=started_at,
                run_index=0 if label == "unpolarized" else 1,
                run_count=2,
                log_path=run_dir / "compile.log",
            ),
        )
    ensure_compiled(poldis_dir, run_dir, dry_run=dry_run)
    if parallel_mode:
        ensure_ran(run_dir, dry_run=dry_run)
    else:
        run_poldis_with_progress(
            base_dir=base_dir,
            tag=tag,
            window=window,
            run_dir=run_dir,
            variant_label=label,
            events=events,
            run_index=0 if label == "unpolarized" else 1,
            run_count=2,
            started_at=started_at,
            dry_run=dry_run,
        )
    return label


def load_variant_payload(run_dir: Path, ipol: int) -> dict:
    log_path = run_dir / "run.log"
    if not log_path.exists():
        raise FileNotFoundError(f"Missing POLDIS run log {log_path}")
    totals = parse_poldis_totals(log_path.read_text(), context=str(log_path))
    top_path = run_dir / expected_top_name(ipol)
    payload = {
        "log": str(log_path),
        "top": str(top_path) if top_path.exists() else None,
        "totals": {order: measurement_payload(totals[order]) for order in ("LO", "NLO", "NNLO")},
    }
    return payload


def build_broad_comparisons(payload: dict) -> dict:
    comparisons: dict = {"unpolarized": {}, "polarized": {}}
    for order in ("LO", "NLO", "NNLO"):
        unpol_meas = measurement_from_payload(payload["unpolarized"][order])
        pol_meas = measurement_from_payload(payload["polarized"][order])
        ref_unpol = POLDIS_UNPOL_REFS["GAMMA"][order]
        ref_pol = POLDIS_POL_REFS["GAMMA"][order]
        comparisons["unpolarized"][order] = measurement_payload(sub_measurement(unpol_meas, ref_unpol))
        comparisons["polarized"][order] = measurement_payload(sub_measurement(pol_meas, ref_pol))
    return comparisons


def build_reference_payload(
    base_dir: Path,
    tag: str,
    window: CutWindow,
    events: int,
    pdf_profile: str,
    pdfs: Mapping[str, str],
) -> dict:
    unpol_dir = variant_dir(base_dir, tag, window, "unpolarized")
    pol_dir = variant_dir(base_dir, tag, window, "polarized")
    unpol_payload = load_variant_payload(unpol_dir, ipol=0)
    pol_payload = load_variant_payload(pol_dir, ipol=1)
    payload = {
        "tag": tag,
        "mode": "poldis-gamma-window-reference",
        "window": asdict(window),
        "events": events,
        "pdf_profile": pdf_profile,
        "pdfs": {
            "unpolarized": pdfs["unpolarized"],
            "polarized_diff": pdfs["polarized_diff"],
        },
        "work_dir": str(work_dir(base_dir, tag, window)),
        "runs": {
            "unpolarized": unpol_payload,
            "polarized": pol_payload,
        },
        "unpolarized": unpol_payload["totals"],
        "polarized": pol_payload["totals"],
    }
    if window.label == "broad" and is_builtin_hybrid_pdf_choice(pdfs):
        payload["builtin_broad_comparisons"] = build_broad_comparisons(payload)
    return payload


def render_reference_text(payload: dict) -> str:
    window = payload["window"]
    lines = [
        "POLDIS GAMMA window reference",
        "=============================",
        f"Tag: {payload['tag']}",
        f"Window: {window['label']}",
        f"Q^2 in [{window['q2_min']:.0f}, {window['q2_max']:.0f}] GeV^2, y in [{window['y_min']:.1f}, {window['y_max']:.1f}]",
        f"Events per run: {payload['events']}",
        f"PDF profile: {payload['pdf_profile']}",
        f"Unpolarized PDF: {payload['pdfs']['unpolarized']}",
        f"Polarized diff PDF: {payload['pdfs']['polarized_diff']}",
        "",
        "Unpolarized totals",
        "------------------",
        f"LO:   {fmt_measurement(payload['unpolarized']['LO'])}",
        f"NLO:  {fmt_measurement(payload['unpolarized']['NLO'])}",
        f"NNLO: {fmt_measurement(payload['unpolarized']['NNLO'])}",
        f"log:  {payload['runs']['unpolarized']['log']}",
        f"top:  {payload['runs']['unpolarized']['top'] or 'missing'}",
        "",
        "Polarized totals",
        "----------------",
        f"LO:   {fmt_measurement(payload['polarized']['LO'])}",
        f"NLO:  {fmt_measurement(payload['polarized']['NLO'])}",
        f"NNLO: {fmt_measurement(payload['polarized']['NNLO'])}",
        f"log:  {payload['runs']['polarized']['log']}",
        f"top:  {payload['runs']['polarized']['top'] or 'missing'}",
    ]
    if "builtin_broad_comparisons" in payload:
        comparisons = payload["builtin_broad_comparisons"]
        lines.extend(
            [
                "",
                "Broad sanity check vs built-in DISPOL references",
                "-----------------------------------------------",
            ]
        )
        for order in ("LO", "NLO", "NNLO"):
            lines.append(
                f"Unpolarized {order} - stored ref: {fmt_measurement(comparisons['unpolarized'][order])}"
            )
            lines.append(
                f"Polarized {order} - stored ref:   {fmt_measurement(comparisons['polarized'][order])}"
            )
    return "\n".join(lines).rstrip() + "\n"


def write_results(base_dir: Path, tag: str, window: CutWindow, payload: dict) -> None:
    output_dir = work_dir(base_dir, tag, window)
    output_dir.mkdir(parents=True, exist_ok=True)
    text = render_reference_text(payload)
    reference_txt_path(base_dir, tag, window).write_text(text)
    reference_json_path(base_dir, tag, window).write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(text, end="")
    print(f"Wrote text reference: {reference_txt_path(base_dir, tag, window)}")
    print(f"Wrote JSON reference: {reference_json_path(base_dir, tag, window)}")


def run(args: argparse.Namespace) -> None:
    base_dir = args.base_dir.resolve()
    poldis_dir = args.poldis_dir.resolve()
    window = WINDOWS[args.window]
    pdfs = resolve_pdf_choice(args.pdf_profile, args.unpolarized_pdf, args.polarized_diff_pdf)
    started_at = time.time()
    if not args.dry_run:
        write_monitor_files(
            base_dir,
            args.tag,
            window,
            build_monitor_payload(
                tag=args.tag,
                window=window,
                phase="initializing",
                variant="none",
                started_at=started_at,
                run_index=0,
                run_count=2,
            ),
        )

    if not args.collect_only:
        variants = (("unpolarized", 0), ("polarized", 1))
        parallel_variants = max(1, min(args.jobs, len(variants))) > 1

        if parallel_variants and not args.dry_run:
            write_monitor_files(
                base_dir,
                args.tag,
                window,
                build_monitor_payload(
                    tag=args.tag,
                    window=window,
                    phase="running-parallel",
                    variant="parallel",
                    started_at=started_at,
                    run_index=0,
                    run_count=2,
                ),
            )

        if parallel_variants:
            completed = 0
            with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, min(args.jobs, len(variants)))) as executor:
                future_map = {
                    executor.submit(
                        process_variant,
                        base_dir=base_dir,
                        tag=args.tag,
                        window=window,
                        poldis_dir=poldis_dir,
                        label=label,
                        ipol=ipol,
                        events=args.events,
                        pdfs=pdfs,
                        started_at=started_at,
                        dry_run=args.dry_run,
                        parallel_mode=True,
                    ): label
                    for label, ipol in variants
                }
                for future in concurrent.futures.as_completed(future_map):
                    label = future.result()
                    completed += 1
                    if not args.dry_run:
                        write_monitor_files(
                            base_dir,
                            args.tag,
                            window,
                            build_monitor_payload(
                                tag=args.tag,
                                window=window,
                                phase="running-parallel",
                                variant=label,
                                started_at=started_at,
                                run_index=completed,
                                run_count=2,
                                log_path=variant_dir(base_dir, args.tag, window, label) / "run.log",
                            ),
                        )
                    print(f"[progress:poldis-gamma-window-reference] completed {label} ({completed}/2)", flush=True)
        else:
            for run_index, (label, ipol) in enumerate(variants):
                process_variant(
                    base_dir=base_dir,
                    tag=args.tag,
                    window=window,
                    poldis_dir=poldis_dir,
                    label=label,
                    ipol=ipol,
                    events=args.events,
                    pdfs=pdfs,
                    started_at=started_at,
                    dry_run=args.dry_run,
                    parallel_mode=False,
                )
        if args.dry_run:
            return

    payload = build_reference_payload(
        base_dir,
        args.tag,
        window,
        events=args.events,
        pdf_profile=args.pdf_profile,
        pdfs=pdfs,
    )
    write_results(base_dir, args.tag, window, payload)
    write_monitor_files(
        base_dir,
        args.tag,
        window,
        build_monitor_payload(
            tag=args.tag,
            window=window,
            phase="complete",
            variant="all",
            started_at=started_at,
            run_index=2,
            run_count=2,
            events_total=args.events,
            events_done=args.events,
            log_path=reference_txt_path(base_dir, args.tag, window),
        ),
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.base_dir.resolve().exists():
        raise SystemExit(f"Base directory does not exist: {args.base_dir}")
    if not args.poldis_dir.resolve().exists():
        raise SystemExit(f"POLDIS directory does not exist: {args.poldis_dir}")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
