#!/usr/bin/env python3.10
"""Print raw-POWHEG selected histogram pieces split by order and channel.

This helper works directly from the raw-POWHEG shard YODAs written by
``run_validation_campaign.py ... --raw-powheg``. It avoids ``yoda.read()``
and instead parses the Estimate1D text blocks explicitly, which is more robust
for the raw sidecar files produced on the server.

For a requested histogram such as ``DQ2``, it prints the selected:

  - POSNLO Compton piece
  - POSNLO BGF piece
  - NEGNLO Compton piece
  - NEGNLO BGF piece
  - signed channel pieces POS - NEG
  - signed combined total

bin by bin, together with a closure check that the signed total is equal to the
sum of the signed channel pieces.
"""

from __future__ import annotations

import argparse
import ast
import gzip
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from run_validation_campaign import collect_raw_powheg_inputs, resolve_campaign_dir


RAW_CHANNELS = ("Compton", "BGF")


@dataclass(frozen=True)
class HistogramMeasurement:
    edges: Tuple[float, ...]
    values: Tuple[float, ...]
    errors: Tuple[float, ...]

    @property
    def nbins(self) -> int:
        return len(self.values)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag", help="Campaign tag, e.g. rivetfo_raw12.")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="DISPOL base directory containing campaigns/ (default: this script's directory).",
    )
    parser.add_argument("--setup", default="ALL", help="Physics setup to inspect (default: ALL).")
    parser.add_argument(
        "--analysis",
        default="MC_DIS_BREIT",
        help="Analysis prefix in the raw YODA paths (default: MC_DIS_BREIT).",
    )
    parser.add_argument(
        "--hist",
        default="DQ2",
        help=(
            "Analyzed histogram suffix to print, e.g. DQ2, Q2, DpT1. "
            "For D/ALL histograms this helper automatically maps back to the "
            "raw base histogram before combining helicities."
        ),
    )
    parser.add_argument("--xmin", type=float, default=None, help="Only print bins with xmax > xmin.")
    parser.add_argument("--xmax", type=float, default=None, help="Only print bins with xmin < xmax.")
    parser.add_argument(
        "--precision",
        type=int,
        default=6,
        help="Digits after the decimal point for floating-point output (default: 6).",
    )
    parser.add_argument(
        "--all-bins",
        action="store_true",
        help="Print all bins, even if every printed contribution is exactly zero.",
    )
    return parser.parse_args()


def resolve_hist_request(hist: str) -> tuple[str, str]:
    if hist.startswith("ALL") and len(hist) > 3:
        return "all", hist[3:]
    if hist.startswith("D") and len(hist) > 1:
        return "delta", hist[1:]
    return "sigma", hist


def resolve_artifact_path(path: Path, base_dir: Path, campaign_dir: Path) -> Path:
    path = path.expanduser()
    if path.exists():
        return path.resolve()

    candidates = (
        base_dir / path.name,
        campaign_dir / path.name,
        campaign_dir / "raw-powheg" / path.name,
        campaign_dir / "analysis" / path.name,
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    raise FileNotFoundError(
        f"Could not resolve raw POWHEG artifact '{path}'. Tried: "
        + ", ".join(str(candidate) for candidate in candidates)
    )


def normalize_input_paths(
    inputs: Dict[str, List[Path]],
    base_dir: Path,
    campaign_dir: Path,
) -> Dict[str, List[Path]]:
    return {
        helicity: [resolve_artifact_path(path, base_dir, campaign_dir) for path in paths]
        for helicity, paths in inputs.items()
    }


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
    edges: Optional[List[float]] = None
    values: List[float] = []
    errors: List[float] = []
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


def combine_shard_measurements(
    paths: Sequence[Path],
    object_path: str,
    channel: Optional[str] = None,
) -> HistogramMeasurement:
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
    values: List[float] = []
    errors: List[float] = []
    nshards = len(shards)
    for index in range(nbins):
        mean_value = math.fsum(shard.values[index] for shard in shards) / nshards
        mean_error = math.sqrt(math.fsum(shard.errors[index] ** 2 for shard in shards)) / nshards
        values.append(mean_value)
        errors.append(mean_error)
    return HistogramMeasurement(edges=edges, values=tuple(values), errors=tuple(errors))


def subtract_measurements(pos: HistogramMeasurement, neg: Optional[HistogramMeasurement]) -> HistogramMeasurement:
    if neg is None:
        return pos
    if pos.edges != neg.edges:
        raise RuntimeError("Cannot subtract measurements with different binning")
    values = []
    errors = []
    for pos_value, pos_error, neg_value, neg_error in zip(pos.values, pos.errors, neg.values, neg.errors):
        values.append(pos_value - neg_value)
        errors.append(math.hypot(pos_error, neg_error))
    return HistogramMeasurement(edges=pos.edges, values=tuple(values), errors=tuple(errors))


def ratio_measurement(numerator: HistogramMeasurement, denominator: HistogramMeasurement) -> HistogramMeasurement:
    if numerator.edges != denominator.edges:
        raise RuntimeError("Cannot form a ratio for measurements with different binning")
    values = []
    errors = []
    for num, err_num, den, err_den in zip(
        numerator.values, numerator.errors, denominator.values, denominator.errors
    ):
        if (not math.isfinite(num) or not math.isfinite(err_num) or
                not math.isfinite(den) or not math.isfinite(err_den) or den == 0.0):
            values.append(float("nan"))
            errors.append(0.0)
            continue
        value = num / den
        error = math.hypot(err_num / den, num * err_den / (den * den))
        values.append(value)
        errors.append(error)
    return HistogramMeasurement(edges=numerator.edges, values=tuple(values), errors=tuple(errors))


def required_helicities(setup: str) -> Tuple[str, ...]:
    if setup == "GAMMA":
        return ("00", "PP", "PM")
    return ("00", "PP", "PM", "MP", "MM")


def load_helicity_measurements(
    base_dir: Path,
    tag: str,
    setup: str,
    order: str,
    analysis: str,
    base_label: str,
    channel: Optional[str] = None,
) -> Dict[str, HistogramMeasurement]:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    inputs = collect_raw_powheg_inputs(campaign_dir, tag, setup, order, channel=channel)
    inputs = normalize_input_paths(inputs, base_dir, campaign_dir)
    object_path = f"/{analysis}/{base_label}"
    return {
        helicity: combine_shard_measurements(paths, object_path, channel=channel)
        for helicity, paths in inputs.items()
    }


def combine_setup_measurements(
    setup: str,
    pos_map: Dict[str, HistogramMeasurement],
    neg_map: Optional[Dict[str, HistogramMeasurement]] = None,
) -> tuple[HistogramMeasurement, HistogramMeasurement, HistogramMeasurement]:
    neg_map = neg_map or {}
    helicities = required_helicities(setup)
    edges = pos_map[helicities[0]].edges

    sigma_values: List[float] = []
    sigma_errors: List[float] = []
    delta_values: List[float] = []
    delta_errors: List[float] = []

    for bin_index in range(len(edges) - 1):
        def value_error(helicity: str) -> tuple[float, float]:
            pos = pos_map[helicity]
            neg = neg_map.get(helicity)
            value = pos.values[bin_index]
            error = pos.errors[bin_index]
            if neg is not None:
                value -= neg.values[bin_index]
                error = math.hypot(error, neg.errors[bin_index])
            return value, error

        n00, e00 = value_error("00")
        npp, epp = value_error("PP")
        npm, epm = value_error("PM")

        if setup == "GAMMA":
            sigma_value = n00
            sigma_error = e00
            delta_value = 0.5 * math.fsum((npp, -npm))
            delta_error = 0.5 * math.sqrt(math.fsum((epp**2, epm**2)))
        else:
            nmp, emp = value_error("MP")
            nmm, emm = value_error("MM")
            sigma_value = n00
            sigma_error = e00
            delta_value = 0.25 * math.fsum((npp, nmm, -npm, -nmp))
            delta_error = 0.25 * math.sqrt(math.fsum((epp**2, epm**2, emp**2, emm**2)))

        sigma_values.append(sigma_value)
        sigma_errors.append(sigma_error)
        delta_values.append(delta_value)
        delta_errors.append(delta_error)

    sigma = HistogramMeasurement(edges=edges, values=tuple(sigma_values), errors=tuple(sigma_errors))
    delta = HistogramMeasurement(edges=edges, values=tuple(delta_values), errors=tuple(delta_errors))
    ratio = ratio_measurement(delta, sigma)
    return sigma, delta, ratio


def select_histogram(
    mode: str,
    sigma: HistogramMeasurement,
    delta: HistogramMeasurement,
    ratio: HistogramMeasurement,
) -> HistogramMeasurement:
    if mode == "sigma":
        return sigma
    if mode == "delta":
        return delta
    if mode == "all":
        return ratio
    raise ValueError(f"Unsupported histogram mode: {mode}")


def format_float(value: float, precision: int) -> str:
    if math.isnan(value):
        return "nan"
    return f"{value:.{precision}f}"


def main() -> int:
    args = parse_args()
    base_dir = args.base_dir.resolve()
    setup = args.setup.upper()
    hist_mode, base_label = resolve_hist_request(args.hist)
    campaign_dir = resolve_campaign_dir(base_dir, args.tag)

    pos_total_map = load_helicity_measurements(base_dir, args.tag, setup, "POSNLO", args.analysis, base_label)
    neg_total_map = load_helicity_measurements(base_dir, args.tag, setup, "NEGNLO", args.analysis, base_label)
    pos_total_sigma, pos_total_delta, pos_total_ratio = combine_setup_measurements(setup, pos_total_map)
    neg_total_sigma, neg_total_delta, neg_total_ratio = combine_setup_measurements(setup, neg_total_map)
    nlo_total_sigma, nlo_total_delta, nlo_total_ratio = combine_setup_measurements(setup, pos_total_map, neg_total_map)

    pos_channel_hists = {}
    neg_channel_hists = {}
    nlo_channel_hists = {}
    for channel in RAW_CHANNELS:
        pos_map = load_helicity_measurements(
            base_dir, args.tag, setup, "POSNLO", args.analysis, base_label, channel=channel
        )
        neg_map = load_helicity_measurements(
            base_dir, args.tag, setup, "NEGNLO", args.analysis, base_label, channel=channel
        )
        pos_sigma, pos_delta, pos_ratio = combine_setup_measurements(setup, pos_map)
        neg_sigma, neg_delta, neg_ratio = combine_setup_measurements(setup, neg_map)
        nlo_sigma, nlo_delta, nlo_ratio = combine_setup_measurements(setup, pos_map, neg_map)
        pos_channel_hists[channel] = select_histogram(hist_mode, pos_sigma, pos_delta, pos_ratio)
        neg_channel_hists[channel] = select_histogram(hist_mode, neg_sigma, neg_delta, neg_ratio)
        nlo_channel_hists[channel] = select_histogram(hist_mode, nlo_sigma, nlo_delta, nlo_ratio)

    pos_total = select_histogram(hist_mode, pos_total_sigma, pos_total_delta, pos_total_ratio)
    neg_total = select_histogram(hist_mode, neg_total_sigma, neg_total_delta, neg_total_ratio)
    nlo_total = select_histogram(hist_mode, nlo_total_sigma, nlo_total_delta, nlo_total_ratio)

    print(f"Campaign: {campaign_dir}")
    print(f"Setup: {setup}")
    print(f"Histogram: /{args.analysis}/{args.hist}")
    print(f"Raw base histogram: /{args.analysis}/{base_label}")
    print("Columns:")
    print("  POS_C / POS_B = positive-order channel-resolved selected contributions")
    print("  NEG_C / NEG_B = negative-order channel-resolved selected contributions")
    print("  NLO_C / NLO_B = signed channel pieces (POS - NEG)")
    print("  NLO_T         = signed combined total (POS - NEG)")
    print("  closure       = NLO_T - (NLO_C + NLO_B)")
    print()
    print(
        "bin  range"
        "                          POS_C       POS_B       NEG_C       NEG_B"
        "       NLO_C       NLO_B       NLO_T      closure"
    )

    for index in range(pos_total.nbins):
        xlo = pos_total.edges[index]
        xhi = pos_total.edges[index + 1]
        if args.xmin is not None and xhi <= args.xmin:
            continue
        if args.xmax is not None and xlo >= args.xmax:
            continue

        pos_c = pos_channel_hists["Compton"].values[index]
        pos_b = pos_channel_hists["BGF"].values[index]
        neg_c = neg_channel_hists["Compton"].values[index]
        neg_b = neg_channel_hists["BGF"].values[index]
        nlo_c = nlo_channel_hists["Compton"].values[index]
        nlo_b = nlo_channel_hists["BGF"].values[index]
        nlo_t = nlo_total.values[index]
        closure = nlo_t - (nlo_c + nlo_b)

        values = (pos_c, pos_b, neg_c, neg_b, nlo_c, nlo_b, nlo_t, closure)
        if not args.all_bins and all(value == 0.0 for value in values):
            continue

        print(
            f"{index + 1:3d}  [{xlo:7.2f}, {xhi:7.2f})"
            f"  {format_float(pos_c, args.precision):>10}"
            f"  {format_float(pos_b, args.precision):>10}"
            f"  {format_float(neg_c, args.precision):>10}"
            f"  {format_float(neg_b, args.precision):>10}"
            f"  {format_float(nlo_c, args.precision):>10}"
            f"  {format_float(nlo_b, args.precision):>10}"
            f"  {format_float(nlo_t, args.precision):>10}"
            f"  {format_float(closure, args.precision):>10}"
        )

    pos_total_sum = math.fsum(pos_total.values)
    neg_total_sum = math.fsum(neg_total.values)
    nlo_total_sum = math.fsum(nlo_total.values)
    print()
    print("Totals (sum of printed bin heights, not width-weighted integrals):")
    print(f"  POS_total = {format_float(pos_total_sum, args.precision)}")
    print(f"  NEG_total = {format_float(neg_total_sum, args.precision)}")
    print(f"  NLO_total = {format_float(nlo_total_sum, args.precision)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
