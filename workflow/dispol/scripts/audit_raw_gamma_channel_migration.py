#!/usr/bin/env python3.10
"""Audit raw GAMMA channel migration between winner-level and selected Q2 bins.

This helper is a focused follow-up to the raw POWHEG channel studies. It is
designed to answer one specific question:

  Is the low-Q2 polarized channel distortion already present in the raw winner
  sample, or does it appear only after the selected MC_DIS_BREIT-like cuts?

For one GAMMA channel and one order branch, it prints coarse Q2-bin tables for

  - winner-level channel weights from the summary CSV (DIS-window only)
  - selected channel weights from the raw channel YODAs
  - per-helicity efficiencies epsilon_h = selected_h / winner_h
  - the polarized combination D = 0.5 * (PP - PM) before and after selection
  - the channel asymmetry A = D / 00 before and after selection

The helper is intentionally self-contained and does not import ``yoda`` or
``run_validation_campaign.py`` so it can be run with a plain Python + gzip
environment.
"""

from __future__ import annotations

import argparse
import ast
import csv
import gzip
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence


HELICITIES = ("00", "PP", "PM")


@dataclass(frozen=True)
class HistogramMeasurement:
    edges: tuple[float, ...]
    values: tuple[float, ...]
    errors: tuple[float, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag", help="Campaign tag, e.g. rivetfo_raw13.")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="DISPOL base directory containing campaigns/ (default: this script's directory).",
    )
    parser.add_argument(
        "--order",
        choices=("POSNLO", "NEGNLO"),
        default="POSNLO",
        help="Order branch to inspect (default: POSNLO).",
    )
    parser.add_argument(
        "--channel",
        choices=("Compton", "BGF"),
        default="Compton",
        help="Winner channel to inspect (default: Compton).",
    )
    parser.add_argument(
        "--analysis",
        default="MC_DIS_BREIT",
        help="Analysis prefix in the raw channel YODA paths (default: MC_DIS_BREIT).",
    )
    parser.add_argument("--q2-min", type=float, default=100.0, help="Only print bins with upper edge above this value.")
    parser.add_argument("--q2-max", type=float, default=1000.0, help="Only print bins with lower edge below this value.")
    parser.add_argument(
        "--precision",
        type=int,
        default=6,
        help="Digits after the decimal point for floating-point output (default: 6).",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=0,
        help="If > 0, print a progress line every N processed raw combined entries.",
    )
    return parser.parse_args()


def resolve_campaign_dir(base_dir: Path, tag: str) -> Path:
    campaign_dir = base_dir / "campaigns" / tag
    if campaign_dir.exists():
        return campaign_dir.resolve()
    direct = Path(tag).expanduser()
    if direct.exists():
        return direct.resolve()
    raise FileNotFoundError(f"Could not resolve campaign directory for tag '{tag}' under {base_dir}")


def load_manifest(campaign_dir: Path) -> dict:
    manifest_path = campaign_dir / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing manifest: {manifest_path}")
    return json.loads(manifest_path.read_text())


def resolve_artifact_path(path: Path, base_dir: Path, campaign_dir: Path) -> Path:
    path = path.expanduser()
    if path.exists():
        return path.resolve()

    candidates = (
        base_dir / path.name,
        campaign_dir / path.name,
        campaign_dir / "raw-powheg" / path.name,
        campaign_dir / "analysis" / path.name,
        campaign_dir / "launcher-logs" / path.name,
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(
        f"Could not resolve raw POWHEG artifact '{path}'. Tried: "
        + ", ".join(str(candidate) for candidate in candidates)
    )


def iter_manifest_entries(manifest: dict, setup: str, order: str, helicity: str) -> Iterable[dict]:
    prefix = f"DIS-POL-POWHEG_{helicity}-{order}-{setup}"
    for entry in manifest.get("raw_powheg_yoda", []):
        if not isinstance(entry, dict):
            continue
        if entry.get("channel") not in (None, ""):
            continue
        if entry.get("returncode") != 0:
            continue
        logical_run = str(entry.get("logical_run", ""))
        if logical_run.startswith(prefix):
            yield entry


def stdout_path_for_yoda(yoda_path: Path) -> Path:
    stdout_path = yoda_path.with_suffix("").with_suffix(".yoda.stdout")
    if not stdout_path.exists():
        raise FileNotFoundError(f"Missing raw POWHEG stdout sidecar: {stdout_path}")
    return stdout_path


def parse_parsed_event_count(stdout_path: Path) -> int:
    prefix = "[raw-powheg] parsed events:"
    for line in stdout_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(prefix):
            return int(line.split(":", 1)[1].strip())
    raise RuntimeError(f"Could not find parsed event count in {stdout_path}")


def parse_cross_section_pb(command: Sequence[object]) -> float:
    for index, token in enumerate(command):
        if token == "--cross-section-pb" and index + 1 < len(command):
            return float(command[index + 1])
    raise ValueError(f"Could not find --cross-section-pb in command: {command}")


def channel_output_path_for_entry(entry: dict, channel: str, base_dir: Path, campaign_dir: Path) -> Path:
    output_file = entry.get("output_file")
    if not isinstance(output_file, str) or not output_file:
        raise FileNotFoundError(f"Raw POWHEG entry is missing output_file: {entry}")
    combined_path = resolve_artifact_path(Path(output_file), base_dir, campaign_dir)
    suffix = ".rawpowheg.yoda.gz"
    if not combined_path.name.endswith(suffix):
        raise FileNotFoundError(f"Unexpected combined raw POWHEG filename: {combined_path}")
    channel_name = combined_path.name[: -len(suffix)] + f".rawpowheg.{channel}.yoda.gz"
    channel_path = combined_path.with_name(channel_name)
    if not channel_path.exists():
        raise FileNotFoundError(f"Missing channel raw POWHEG YODA: {channel_path}")
    return channel_path


def channel_normalization_scale(yoda_path: Path, channel: Optional[str]) -> float:
    if channel is None:
        return 1.0

    suffix = f".rawpowheg.{channel}.yoda.gz"
    if not yoda_path.name.endswith(suffix):
        raise FileNotFoundError(f"Unexpected raw POWHEG channel filename: {yoda_path}")
    combined_name = yoda_path.name[: -len(suffix)] + ".rawpowheg.yoda.gz"
    combined_path = yoda_path.with_name(combined_name)
    if not combined_path.exists():
        raise FileNotFoundError(f"Missing combined raw POWHEG YODA partner: {combined_path}")

    combined_parsed = parse_parsed_event_count(stdout_path_for_yoda(combined_path))
    channel_parsed = parse_parsed_event_count(stdout_path_for_yoda(yoda_path))
    if combined_parsed <= 0 or channel_parsed <= 0:
        return 1.0
    return channel_parsed / combined_parsed


def parse_error_token(token: str) -> float:
    lowered = token.lower()
    if lowered in {"---", "nan"}:
        return 0.0
    value = float(token)
    return 0.0 if not math.isfinite(value) else abs(value)


def read_estimate1d(yoda_path: Path, object_path: str, scale: float = 1.0) -> HistogramMeasurement:
    opener = gzip.open if yoda_path.suffix == ".gz" else open
    edges: list[float] | None = None
    values: list[float] = []
    errors: list[float] = []
    in_target = False

    with opener(yoda_path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("BEGIN ") and line.endswith(object_path):
                in_target = True
                continue
            if not in_target:
                continue
            if line.startswith("END "):
                break
            if line.startswith("Edges(A1):"):
                edges = [float(item) for item in ast.literal_eval(line.split(":", 1)[1].strip())]
                continue
            if not line or line.startswith("#") or line == "---":
                continue
            if ":" in line:
                continue
            parts = line.split()
            if not parts:
                continue
            value_token = parts[0]
            down_token = parts[1] if len(parts) > 1 else "0"
            up_token = parts[2] if len(parts) > 2 else down_token
            value = float("nan") if value_token.lower() == "nan" else float(value_token)
            error = 0.5 * (parse_error_token(down_token) + parse_error_token(up_token))
            values.append(value)
            errors.append(error)

    if edges is None:
        raise KeyError(f"Could not find {object_path} in {yoda_path}")
    if len(values) == len(edges) + 1:
        values = values[1:-1]
        errors = errors[1:-1]
    if len(values) != len(edges) - 1:
        raise RuntimeError(
            f"Parsed {len(values)} values for {object_path} but found {len(edges) - 1} bins in {yoda_path}"
        )

    cleaned_values = []
    cleaned_errors = []
    for value, error in zip(values, errors):
        cleaned_values.append(scale * (0.0 if not math.isfinite(value) else value))
        cleaned_errors.append(scale * (0.0 if not math.isfinite(error) else error))
    return HistogramMeasurement(
        edges=tuple(edges),
        values=tuple(cleaned_values),
        errors=tuple(cleaned_errors),
    )


def combine_shard_measurements(paths: Sequence[Path], object_path: str, channel: Optional[str]) -> HistogramMeasurement:
    if not paths:
        raise ValueError(f"No raw YODA paths were provided for {object_path}")
    shards = [
        read_estimate1d(path, object_path, scale=channel_normalization_scale(path, channel))
        for path in paths
    ]
    edges = shards[0].edges
    for shard in shards[1:]:
        if shard.edges != edges:
            raise RuntimeError(f"Mismatched binning for {object_path}")

    nbins = len(edges) - 1
    values: list[float] = []
    errors: list[float] = []
    nshards = len(shards)
    for index in range(nbins):
        mean_value = math.fsum(shard.values[index] for shard in shards) / nshards
        mean_error = math.sqrt(math.fsum(shard.errors[index] ** 2 for shard in shards)) / nshards
        values.append(mean_value)
        errors.append(mean_error)
    return HistogramMeasurement(edges=edges, values=tuple(values), errors=tuple(errors))


def load_selected_measurements(
    *,
    base_dir: Path,
    tag: str,
    order: str,
    channel: str,
    analysis: str,
    progress_every: int = 0,
) -> Mapping[str, HistogramMeasurement]:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    manifest = load_manifest(campaign_dir)
    object_path = f"/{analysis}/Q2"
    out: dict[str, HistogramMeasurement] = {}
    for helicity in HELICITIES:
        paths = [
            channel_output_path_for_entry(entry, channel, base_dir, campaign_dir)
            for entry in iter_manifest_entries(manifest, "GAMMA", order, helicity)
        ]
        out[helicity] = combine_shard_measurements(paths, object_path, channel=channel)
        if progress_every > 0:
            print(
                f"[progress] selected measurements: loaded {helicity} {order} {channel} "
                f"from {len(paths)} shard groups",
                flush=True,
            )
    return out


def new_bin_array(edges: Sequence[float]) -> list[float]:
    return [0.0] * (len(edges) - 1)


def find_bin(edges: Sequence[float], value: float) -> int | None:
    for index in range(len(edges) - 1):
        low = float(edges[index])
        high = float(edges[index + 1])
        if low <= value < high:
            return index
    if value == float(edges[-1]):
        return len(edges) - 2
    return None


def measurement_integrals(measurement: HistogramMeasurement) -> list[float]:
    integrals: list[float] = []
    for index, value in enumerate(measurement.values):
        width = measurement.edges[index + 1] - measurement.edges[index]
        integrals.append(float(value) * float(width))
    return integrals


def load_winner_counts(
    *,
    base_dir: Path,
    tag: str,
    order: str,
    channel: str,
    edges: Sequence[float],
    progress_every: int = 0,
) -> Mapping[str, list[float]]:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    manifest = load_manifest(campaign_dir)
    out = {helicity: new_bin_array(edges) for helicity in HELICITIES}

    processed = 0
    for helicity in HELICITIES:
        for entry in iter_manifest_entries(manifest, "GAMMA", order, helicity):
            processed += 1
            command = entry.get("command")
            if not isinstance(command, list):
                raise ValueError(f"Raw POWHEG entry is missing a command list: {entry}")
            cross_section_pb = parse_cross_section_pb(command)

            summary_csv = entry.get("summary_csv")
            if not isinstance(summary_csv, str) or not summary_csv:
                raise FileNotFoundError(f"Raw POWHEG entry is missing summary_csv: {entry}")
            summary_path = resolve_artifact_path(Path(summary_csv), base_dir, campaign_dir)

            output_file = entry.get("output_file")
            if not isinstance(output_file, str) or not output_file:
                raise FileNotFoundError(f"Raw POWHEG entry is missing output_file: {entry}")
            combined_yoda = resolve_artifact_path(Path(output_file), base_dir, campaign_dir)
            parsed_events = parse_parsed_event_count(stdout_path_for_yoda(combined_yoda))
            if parsed_events <= 0:
                continue

            scale = cross_section_pb / float(parsed_events)
            with summary_path.open("rt", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle):
                    if row.get("winner", "") != channel:
                        continue
                    q2 = float(row.get("Q2", "0") or 0.0)
                    y_value = float(row.get("y", "0") or 0.0)
                    if not (0.2 <= y_value <= 0.6):
                        continue
                    q2_index = find_bin(edges, q2)
                    if q2_index is None:
                        continue
                    out[helicity][q2_index] += scale

            if progress_every > 0 and processed % progress_every == 0:
                print(
                    f"[progress] winner counts: processed {processed} combined entries "
                    f"(latest {helicity} {order})",
                    flush=True,
                )

    return out


def fmt(value: float | None, precision: int) -> str:
    if value is None or not math.isfinite(value):
        return "nan"
    return f"{value:.{precision}f}"


def safe_ratio(num: float, den: float) -> float | None:
    if den == 0.0 or not math.isfinite(num) or not math.isfinite(den):
        return None
    return num / den


def window_sum(edges: Sequence[float], values: Sequence[float], q2_min: float, q2_max: float) -> float:
    total = 0.0
    for low, high, value in zip(edges[:-1], edges[1:], values):
        if low >= q2_min and high <= q2_max:
            total += value
    return total


def print_table(
    *,
    edges: Sequence[float],
    winner: Mapping[str, Sequence[float]],
    selected: Mapping[str, Sequence[float]],
    q2_min: float | None,
    q2_max: float | None,
    precision: int,
) -> None:
    print(
        "bin  Q2-range"
        "                 win00      sel00      eps00"
        "       winPP      selPP      epsPP"
        "       winPM      selPM      epsPM"
        "       D_win      D_sel      A_win      A_sel"
    )
    for index in range(len(edges) - 1):
        low = float(edges[index])
        high = float(edges[index + 1])
        if q2_min is not None and high <= q2_min:
            continue
        if q2_max is not None and low >= q2_max:
            continue

        win00 = float(winner["00"][index])
        winpp = float(winner["PP"][index])
        winpm = float(winner["PM"][index])
        sel00 = float(selected["00"][index])
        selpp = float(selected["PP"][index])
        selpm = float(selected["PM"][index])

        d_win = 0.5 * (winpp - winpm)
        d_sel = 0.5 * (selpp - selpm)
        a_win = safe_ratio(d_win, win00)
        a_sel = safe_ratio(d_sel, sel00)

        print(
            f"{index + 1:3d}  [{low:7.2f}, {high:7.2f})"
            f"  {fmt(win00, precision):>10}"
            f"  {fmt(sel00, precision):>10}"
            f"  {fmt(safe_ratio(sel00, win00), precision):>10}"
            f"  {fmt(winpp, precision):>10}"
            f"  {fmt(selpp, precision):>10}"
            f"  {fmt(safe_ratio(selpp, winpp), precision):>10}"
            f"  {fmt(winpm, precision):>10}"
            f"  {fmt(selpm, precision):>10}"
            f"  {fmt(safe_ratio(selpm, winpm), precision):>10}"
            f"  {fmt(d_win, precision):>10}"
            f"  {fmt(d_sel, precision):>10}"
            f"  {fmt(a_win, precision):>10}"
            f"  {fmt(a_sel, precision):>10}"
        )


def main() -> int:
    args = parse_args()
    base_dir = args.base_dir.resolve()
    selected_measurements = load_selected_measurements(
        base_dir=base_dir,
        tag=args.tag,
        order=args.order,
        channel=args.channel,
        analysis=args.analysis,
        progress_every=args.progress_every,
    )
    edges = selected_measurements["00"].edges
    winner = load_winner_counts(
        base_dir=base_dir,
        tag=args.tag,
        order=args.order,
        channel=args.channel,
        edges=edges,
        progress_every=args.progress_every,
    )
    # The raw channel YODAs store differential bin heights; convert them back to
    # per-bin integrals before comparing against winner-level weighted counts.
    selected = {helicity: measurement_integrals(selected_measurements[helicity]) for helicity in HELICITIES}

    print(f"Campaign: {resolve_campaign_dir(base_dir, args.tag)}")
    print("Setup: GAMMA")
    print(f"Order: {args.order}")
    print(f"Channel: {args.channel}")
    print("Definitions:")
    print("  winner_h  = raw winner-channel weight after the DIS window (from summary_csv)")
    print("  sel_h     = selected channel weight in /MC_DIS_BREIT/Q2 (from raw channel YODA)")
    print("  eps_h     = sel_h / winner_h")
    print("  D         = 0.5 * (PP - PM)")
    print("  A         = D / 00")
    print()

    print_table(
        edges=edges,
        winner=winner,
        selected=selected,
        q2_min=args.q2_min,
        q2_max=args.q2_max,
        precision=args.precision,
    )

    winner_d = [0.5 * (winner["PP"][i] - winner["PM"][i]) for i in range(len(edges) - 1)]
    selected_d = [0.5 * (selected["PP"][i] - selected["PM"][i]) for i in range(len(edges) - 1)]

    print()
    print(f"Window sums in [{args.q2_min:.0f}, {args.q2_max:.0f}] GeV^2:")
    print(f"  winner_00 = {fmt(window_sum(edges, winner['00'], args.q2_min, args.q2_max), args.precision)}")
    print(f"  sel_00    = {fmt(window_sum(edges, selected['00'], args.q2_min, args.q2_max), args.precision)}")
    print(f"  winner_D  = {fmt(window_sum(edges, winner_d, args.q2_min, args.q2_max), args.precision)}")
    print(f"  sel_D     = {fmt(window_sum(edges, selected_d, args.q2_min, args.q2_max), args.precision)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
