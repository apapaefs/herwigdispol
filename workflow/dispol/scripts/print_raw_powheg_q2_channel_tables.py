#!/usr/bin/env python3.10
"""Print raw-POWHEG helicity/channel Q2 tables for winner fractions and efficiencies.

Definitions used here:

  f_c^h(Q2) = N_winner,c^h(Q2) / sum_c' N_winner,c'^h(Q2)
  epsilon_c^h(Q2) = N_selected,c^h(Q2) / N_winner,c^h(Q2)

where:

- c runs over the raw POWHEG winner channels {Compton, BGF}
- h runs over the beam-helicity configurations
- N_winner,c^h counts the raw winning channel after the DIS-window cut only
- N_selected,c^h counts the same raw winner after the full MC_DIS_BREIT-like
  raw-pair selection in ``powheg_raw_momenta_to_yoda.py``

This intentionally excludes Born fallback / winner=None events from f_c^h.
The goal is to isolate channel competition and acceptance inside the emitted
real-emission sample.
"""

from __future__ import annotations

import argparse
import ast
import csv
import gzip
import json
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

from powheg_raw_momenta_to_yoda import BINNING, RAW_CHANNELS, find_bin


DISPLAY_HELICITIES = {
    "ALL": ("PP", "PM", "MP", "MM"),
    "Z": ("PP", "PM", "MP", "MM"),
    "CC": ("PP", "PM", "MP", "MM"),
    "GAMMA": ("PP", "PM"),
}
ORDERS = ("POSNLO", "NEGNLO")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag", help="Campaign tag, e.g. rivetfo_raw08.")
    parser.add_argument(
        "--base-dir",
        default=str(Path(__file__).resolve().parent),
        help="DISPOL base directory containing campaigns/ and the raw logs.",
    )
    parser.add_argument(
        "--setup",
        default="ALL",
        choices=tuple(sorted(DISPLAY_HELICITIES)),
        help="Physics setup to inspect (default: ALL).",
    )
    parser.add_argument(
        "--order",
        default="POSNLO",
        choices=("POSNLO", "NEGNLO", "NLO"),
        help=(
            "Which contribution to inspect. "
            "Use POSNLO/NEGNLO for the physical raw real-emission samples, "
            "or NLO for the signed POSNLO-NEGNLO diagnostic combination "
            "(default: POSNLO)."
        ),
    )
    parser.add_argument(
        "--q2-min",
        type=float,
        default=None,
        help="Only print bins with upper edge above this value.",
    )
    parser.add_argument(
        "--q2-max",
        type=float,
        default=None,
        help="Only print bins with lower edge below this value.",
    )
    parser.add_argument(
        "--all-bins",
        action="store_true",
        help="Print every Q2 bin in the raw-sidecar grid, even if it is empty.",
    )
    parser.add_argument(
        "--precision",
        type=int,
        default=6,
        help="Digits after the decimal point for table entries (default: 6).",
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        default=None,
        help="Optional cap on the number of raw combined entries to process.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=0,
        help="If > 0, print a progress line every N processed entries.",
    )
    return parser.parse_args()


def resolve_campaign_dir(base_dir: Path, tag: str) -> Path:
    campaign_dir = base_dir / "campaigns" / tag
    if campaign_dir.exists():
        return campaign_dir
    direct = Path(tag).expanduser()
    if direct.exists():
        return direct.resolve()
    raise SystemExit(f"Could not resolve campaign directory for tag '{tag}' under {base_dir}.")


def load_manifest(campaign_dir: Path) -> dict:
    manifest_path = campaign_dir / "manifest.json"
    if not manifest_path.exists():
        raise SystemExit(f"Missing manifest: {manifest_path}")
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def parse_cross_section_pb(command: Sequence[str]) -> float:
    for index, token in enumerate(command):
        if token == "--cross-section-pb" and index + 1 < len(command):
            return float(command[index + 1])
    raise ValueError(f"Could not find --cross-section-pb in command: {command}")


def resolve_manifest_path(raw_path: str, base_dir: Path, campaign_dir: Path) -> Path:
    path = Path(raw_path).expanduser()
    if path.exists():
        return path.resolve()

    candidates = (
        base_dir / path.name,
        campaign_dir / path.name,
        campaign_dir / "raw-powheg" / path.name,
        campaign_dir / "launcher-logs" / path.name,
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    raise FileNotFoundError(
        f"Could not resolve artifact '{raw_path}'. Tried: "
        + ", ".join(str(candidate) for candidate in candidates)
    )


def selected_helicities(setup: str) -> Sequence[str]:
    return DISPLAY_HELICITIES[setup]


def new_q2_bin_array() -> List[float]:
    return [0.0] * (len(BINNING["Q2"]) - 1)


def new_nested_arrays(setup: str) -> Dict[str, Dict[str, List[float]]]:
    return {
        helicity: {channel: new_q2_bin_array() for channel in RAW_CHANNELS}
        for helicity in selected_helicities(setup)
    }


def accumulate(target: List[float], source: Sequence[float], factor: float = 1.0) -> None:
    for index, value in enumerate(source):
        target[index] += factor * float(value)


def iter_relevant_manifest_entries(
    manifest: dict,
    setup: str,
) -> Iterable[tuple[str, str, dict]]:
    helicities = set(selected_helicities(setup))
    entries = manifest.get("raw_powheg_yoda", [])
    if not isinstance(entries, list):
        return

    for entry in entries:
        if not isinstance(entry, dict):
            continue
        if entry.get("channel") not in (None, ""):
            continue
        if entry.get("returncode") != 0:
            continue

        logical_run = str(entry.get("logical_run", ""))
        for helicity in helicities:
            for order in ORDERS:
                prefix = f"DIS-POL-POWHEG_{helicity}-{order}-{setup}"
                if logical_run.startswith(prefix):
                    yield helicity, order, entry
                    break


def read_summary_rows(summary_csv: Path) -> List[dict[str, str]]:
    with summary_csv.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def channel_output_path_for_entry(entry: dict, channel: str, base_dir: Path, campaign_dir: Path) -> Path:
    output_file = entry.get("output_file")
    if not isinstance(output_file, str) or not output_file:
        raise FileNotFoundError(f"Raw POWHEG entry is missing output_file: {entry}")
    combined_path = resolve_manifest_path(output_file, base_dir, campaign_dir)
    suffix = ".rawpowheg.yoda.gz"
    if not combined_path.name.endswith(suffix):
        raise FileNotFoundError(f"Unexpected combined raw POWHEG filename: {combined_path}")
    channel_name = combined_path.name[: -len(suffix)] + f".rawpowheg.{channel}.yoda.gz"
    channel_path = combined_path.with_name(channel_name)
    if not channel_path.exists():
        raise FileNotFoundError(f"Missing channel raw POWHEG YODA: {channel_path}")
    return channel_path


def stdout_path_for_entry(entry: dict, base_dir: Path, campaign_dir: Path) -> Path:
    output_file = entry.get("output_file")
    if not isinstance(output_file, str) or not output_file:
        raise FileNotFoundError(f"Raw POWHEG entry is missing output_file: {entry}")
    combined_path = resolve_manifest_path(output_file, base_dir, campaign_dir)
    stdout_path = combined_path.with_suffix("").with_suffix(".yoda.stdout")
    if not stdout_path.exists():
        raise FileNotFoundError(f"Missing raw POWHEG stdout sidecar: {stdout_path}")
    return stdout_path


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


def read_q2_bin_integrals(yoda_path: Path, analysis: str = "MC_DIS_BREIT") -> List[float]:
    object_path = f"/{analysis}/Q2"
    opener = gzip.open if yoda_path.suffix == ".gz" else open
    edges: List[float] | None = None
    values: List[float] = []
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
            token = line.split()[0]
            values.append(0.0 if token.lower() == "nan" else float(token))

    if edges is None:
        raise KeyError(f"Could not find {object_path} in {yoda_path}")
    # The raw Estimate1D text dump includes explicit underflow/overflow rows.
    if len(values) == len(edges) + 1:
        values = values[1:-1]
    if len(values) != len(edges) - 1:
        raise RuntimeError(
            f"Parsed {len(values)} Q2 bin values but found {len(edges) - 1} edges in {yoda_path}"
        )

    integrals: List[float] = []
    for index, value in enumerate(values):
        width = edges[index + 1] - edges[index]
        integrals.append(value * width)
    return integrals


def collect_weighted_counts(
    base_dir: Path,
    campaign_dir: Path,
    manifest: dict,
    setup: str,
    max_entries: int | None = None,
    progress_every: int = 0,
) -> tuple[
    Dict[str, Dict[str, Dict[str, List[float]]]],
    Dict[str, Dict[str, Dict[str, List[float]]]],
    Dict[str, Dict[str, Dict[str, float]]],
]:
    winner_counts = {order: new_nested_arrays(setup) for order in ORDERS}
    selected_counts = {order: new_nested_arrays(setup) for order in ORDERS}
    stats = {
        order: {
            helicity: {"entries": 0.0, "parsed_events": 0.0}
            for helicity in selected_helicities(setup)
        }
        for order in ORDERS
    }

    seen_any = False
    seen_nonzero = False
    processed = 0
    for helicity, order, entry in iter_relevant_manifest_entries(manifest, setup):
        if max_entries is not None and processed >= max_entries:
            break
        seen_any = True
        processed += 1
        stats[order][helicity]["entries"] += 1.0
        command = entry.get("command")
        if not isinstance(command, list):
            raise SystemExit(f"Raw POWHEG entry is missing a command list: {entry}")
        cross_section_pb = parse_cross_section_pb(command)
        summary_csv = entry.get("summary_csv")
        if not isinstance(summary_csv, str) or not summary_csv:
            raise FileNotFoundError(f"Raw POWHEG entry is missing summary_csv: {entry}")
        summary_rows = read_summary_rows(resolve_manifest_path(summary_csv, base_dir, campaign_dir))
        parsed_events = parse_parsed_event_count(stdout_path_for_entry(entry, base_dir, campaign_dir))
        if parsed_events <= 0:
            continue
        seen_nonzero = True
        stats[order][helicity]["parsed_events"] += float(parsed_events)

        scale = cross_section_pb / float(parsed_events)

        for channel in RAW_CHANNELS:
            # Pre-cut winner counts: only the DIS window, before the raw dijet selection.
            for row in summary_rows:
                if row.get("winner", "") != channel:
                    continue
                q2 = float(row.get("Q2", "0") or 0.0)
                y_value = float(row.get("y", "0") or 0.0)
                if not (0.2 <= y_value <= 0.6):
                    continue
                q2_index = find_bin(BINNING["Q2"], q2)
                if q2_index is None:
                    continue
                winner_counts[order][helicity][channel][q2_index] += scale

            # Selected counts come from the already-converted channel YODA.
            channel_yoda = channel_output_path_for_entry(entry, channel, base_dir, campaign_dir)
            selected_integrals = read_q2_bin_integrals(channel_yoda)
            channel_parsed_events = parse_parsed_event_count(stdout_path_for_yoda(channel_yoda))
            if channel_parsed_events <= 0:
                continue

            # The per-channel raw YODAs are produced after filtering to a
            # single winner channel, so they are normalized with that channel's
            # own parsed-event count. Rescale them back to the combined parsed
            # event count so epsilon_c = N_selected,c / N_winner,c is formed in
            # a consistent normalization.
            accumulate(
                selected_counts[order][helicity][channel],
                selected_integrals,
                factor=float(channel_parsed_events) / float(parsed_events),
            )

        if progress_every > 0 and processed % progress_every == 0:
            print(
                f"[progress] processed {processed} raw combined entries "
                f"(latest {helicity} {order})"
            )

    if not seen_any:
        raise SystemExit(
            f"No successful combined raw POWHEG entries found for setup {setup} in {campaign_dir}."
        )
    if not seen_nonzero:
        raise SystemExit(
            "Found matching raw POWHEG entries, but all had zero parsed raw events. "
            "That usually means the converter saw no POWHEG_RAW_MOM lines in the logs."
        )

    return winner_counts, selected_counts, stats


def combine_nlo(
    pos: Dict[str, Dict[str, List[float]]],
    neg: Dict[str, Dict[str, List[float]]],
) -> Dict[str, Dict[str, List[float]]]:
    out: Dict[str, Dict[str, List[float]]] = {}
    for helicity in pos:
        out[helicity] = {}
        for channel in pos[helicity]:
            out[helicity][channel] = [
                float(pval) - float(nval)
                for pval, nval in zip(pos[helicity][channel], neg[helicity][channel])
            ]
    return out


def should_print_bin(
    index: int,
    q2_edges: Sequence[float],
    winner_nlo: Dict[str, Dict[str, List[float]]],
    q2_min: float | None,
    q2_max: float | None,
    print_all_bins: bool,
) -> bool:
    xlow = float(q2_edges[index])
    xhigh = float(q2_edges[index + 1])
    if q2_min is not None and xhigh <= q2_min:
        return False
    if q2_max is not None and xlow >= q2_max:
        return False
    if print_all_bins:
        return True

    total = 0.0
    for helicity in winner_nlo:
        total += sum(abs(winner_nlo[helicity][channel][index]) for channel in RAW_CHANNELS)
    return total > 0.0


def safe_ratio(numerator: float, denominator: float) -> float | None:
    if abs(denominator) <= 1.0e-30:
        return None
    return numerator / denominator


def format_value(value: float | None, precision: int) -> str:
    if value is None:
        return "nan"
    return f"{value:.{precision}f}"


def print_table(
    title: str,
    q2_edges: Sequence[float],
    winner_view: Dict[str, Dict[str, List[float]]],
    selected_view: Dict[str, Dict[str, List[float]]],
    setup: str,
    precision: int,
    q2_min: float | None,
    q2_max: float | None,
    print_all_bins: bool,
    mode: str,
    channel: str,
) -> None:
    helicities = list(selected_helicities(setup))
    print(title)
    print("bin  Q2-range              " + "  ".join(f"{hel:>8s}" for hel in helicities))
    printed_any = False
    for index in range(len(q2_edges) - 1):
        if not should_print_bin(index, q2_edges, winner_view, q2_min, q2_max, print_all_bins):
            continue
        printed_any = True
        row: List[str] = []
        xlow = float(q2_edges[index])
        xhigh = float(q2_edges[index + 1])
        for helicity in helicities:
            winner = winner_view[helicity][channel][index]
            if mode == "epsilon":
                selected = selected_view[helicity][channel][index]
                value = safe_ratio(selected, winner)
            elif mode == "fraction":
                total = sum(winner_view[helicity][name][index] for name in RAW_CHANNELS)
                value = safe_ratio(winner, total)
            else:
                raise ValueError(f"Unknown table mode: {mode}")
            row.append(f"{format_value(value, precision):>8s}")
        print(f"{index + 1:>3d}  [{xlow:7.2f},{xhigh:7.2f})  " + "  ".join(row))
    if not printed_any:
        print("# no nonzero bins for the current order/selection")
    print()


def main() -> int:
    args = parse_args()
    base_dir = Path(args.base_dir).expanduser().resolve()
    campaign_dir = resolve_campaign_dir(base_dir, args.tag)
    manifest = load_manifest(campaign_dir)

    winner_counts, selected_counts, stats = collect_weighted_counts(
        base_dir,
        campaign_dir,
        manifest,
        args.setup,
        max_entries=args.max_entries,
        progress_every=args.progress_every,
    )
    if args.order == "NLO":
        winner_view = combine_nlo(winner_counts["POSNLO"], winner_counts["NEGNLO"])
        selected_view = combine_nlo(selected_counts["POSNLO"], selected_counts["NEGNLO"])
    else:
        winner_view = winner_counts[args.order]
        selected_view = selected_counts[args.order]

    q2_edges = BINNING["Q2"]

    print(f"Campaign: {campaign_dir}")
    print(f"Setup: {args.setup}")
    print(f"Order view: {args.order}")
    print("Definitions:")
    print("  epsilon_c^h(Q2) = N_selected,c^h / N_winner,c^h")
    print("  f_c^h(Q2)       = N_winner,c^h / (N_winner,Compton^h + N_winner,BGF^h)")
    print("  winner sample excludes Born fallback / winner=None events by construction.")
    print("  selected sample uses the current raw-sidecar MC_DIS_BREIT-like cuts.")
    if args.order == "NLO":
        print("  NLO view uses the signed POSNLO - NEGNLO combination before forming ratios.")
        print("  This is a cancellation diagnostic, not a positive-definite efficiency/fraction.")
    else:
        print("  Order view uses the raw positive sample for that contribution only.")
    print("Entry summary:")
    for order in ORDERS:
        for helicity in selected_helicities(args.setup):
            order_stats = stats[order][helicity]
            print(
                f"  {order:>6s} {helicity}:"
                f" entries={int(order_stats['entries'])}"
                f" parsed_events={int(order_stats['parsed_events'])}"
            )
    print()

    print_table(
        "epsilon_Compton^h(Q2)",
        q2_edges,
        winner_view,
        selected_view,
        args.setup,
        args.precision,
        args.q2_min,
        args.q2_max,
        args.all_bins,
        mode="epsilon",
        channel="Compton",
    )
    print_table(
        "epsilon_BGF^h(Q2)",
        q2_edges,
        winner_view,
        selected_view,
        args.setup,
        args.precision,
        args.q2_min,
        args.q2_max,
        args.all_bins,
        mode="epsilon",
        channel="BGF",
    )
    print_table(
        "f_Compton^h(Q2)",
        q2_edges,
        winner_view,
        selected_view,
        args.setup,
        args.precision,
        args.q2_min,
        args.q2_max,
        args.all_bins,
        mode="fraction",
        channel="Compton",
    )
    print_table(
        "f_BGF^h(Q2)",
        q2_edges,
        winner_view,
        selected_view,
        args.setup,
        args.precision,
        args.q2_min,
        args.q2_max,
        args.all_bins,
        mode="fraction",
        channel="BGF",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
