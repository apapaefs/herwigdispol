#!/usr/bin/env python3.10
"""Compare the selected GAMMA fixed-order pattern between POLDIS and raw Herwig.

This helper is the direct follow-up to the focused GAMMA channel-separated
POLDIS runner.
It compares

  POLDIS:   selected / nlo_real_tree

against the closest raw-Herwig objects:

  Herwig:   setup=GAMMA, order=POSNLO, channel={BGF|Compton}, hist=DQ2

The raw-Herwig YODA stores differential bin heights, while the focused POLDIS
reference stores coarse-bin integrals in pb. The helper therefore rebins the
raw-Herwig selected channel contribution by integrating the piecewise-constant
histogram over the target coarse Q^2 bins before comparing sizes.

For context it also prints the rebinned raw-Herwig

  setup=GAMMA, order=NEGNLO, channel={BGF|Compton}, hist=DQ2

piece and the POS-NEG combination, which helps tell us whether the comparison
surface is dominated by the positive-real tree or whether subtraction-side
effects are numerically important.
"""

from __future__ import annotations

import argparse
import ast
import gzip
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence


DEFAULT_BASE_DIR = Path(__file__).resolve().parent
HERWIG_CHANNELS = ("BGF", "Compton")

_FLOAT_RE = r"[+-]?\d+(?:\.\d*)?(?:[DdEe][+-]?\d+)?"
_TXT_ROW_RE = re.compile(
    r"^\s*(?P<q2_low>" + _FLOAT_RE + r")\s+"
    r"(?P<q2_high>" + _FLOAT_RE + r")\s+"
    r"(?P<value>" + _FLOAT_RE + r")\s+"
    r"(?P<error>" + _FLOAT_RE + r")\s*$"
)


@dataclass(frozen=True)
class BinIntegral:
    q2_low: float
    q2_high: float
    value_pb: float
    error_pb: float


@dataclass(frozen=True)
class HistogramMeasurement:
    edges: tuple[float, ...]
    values: tuple[float, ...]
    errors: tuple[float, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("herwig_tag", help="Raw-Herwig campaign tag, e.g. rivetfo_raw13.")
    parser.add_argument(
        "--poldis-reference",
        type=Path,
        required=True,
        help="Focused channel-separated POLDIS reference (.json preferred, .txt supported).",
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="DISPOL base directory containing campaigns/ (default: this script's directory).",
    )
    parser.add_argument(
        "--herwig-channel",
        choices=HERWIG_CHANNELS,
        default="BGF",
        help="Raw-Herwig selected channel to compare against the POLDIS reference (default: BGF).",
    )
    parser.add_argument(
        "--analysis",
        default="MC_DIS_BREIT",
        help="Raw-Herwig analysis prefix (default: MC_DIS_BREIT).",
    )
    parser.add_argument(
        "--precision",
        type=int,
        default=6,
        help="Digits after the decimal point for printed floating-point values (default: 6).",
    )
    parser.add_argument(
        "--shape-window-min",
        type=float,
        default=100.0,
        help="Lower edge of the shape-only comparison window (default: 100).",
    )
    parser.add_argument(
        "--shape-window-max",
        type=float,
        default=1000.0,
        help="Upper edge of the shape-only comparison window (default: 1000).",
    )
    return parser.parse_args()


def parse_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def maybe_json_sibling(path: Path) -> Path | None:
    if path.suffix.lower() == ".json":
        return path
    sibling = path.with_suffix(".json")
    return sibling if sibling.exists() else None


def load_poldis_selected_real_tree(path: Path) -> list[BinIntegral]:
    json_path = maybe_json_sibling(path)
    if json_path is not None and json_path.exists():
        payload = json.loads(json_path.read_text())
        rows = (
            payload.get("q2_component_tables", {})
            .get("selected", {})
            .get("nlo_real_tree", [])
        )
        result = [
            BinIntegral(
                q2_low=float(row["q2_low"]),
                q2_high=float(row["q2_high"]),
                value_pb=float(row["value_pb"]),
                error_pb=float(row["error_pb"]),
            )
            for row in rows
        ]
        if result:
            return result

    text = path.read_text()
    section_header = "Selected Q2 tables"
    component_header = "nlo_real_tree"
    lines = text.splitlines()
    in_selected = False
    in_component = False
    rows: list[BinIntegral] = []
    for raw_line in lines:
        line = raw_line.rstrip()
        stripped = line.strip()
        if stripped == section_header:
            in_selected = True
            in_component = False
            continue
        if not in_selected:
            continue
        if stripped.endswith("Q2 tables") and stripped != section_header:
            break
        if stripped == component_header:
            in_component = True
            continue
        if in_component and stripped and not line.startswith(" "):
            # Next top-level component header or end of section.
            in_component = False
        if not in_component:
            continue
        match = _TXT_ROW_RE.match(stripped)
        if not match:
            continue
        rows.append(
            BinIntegral(
                q2_low=parse_float(match.group("q2_low")),
                q2_high=parse_float(match.group("q2_high")),
                value_pb=parse_float(match.group("value")),
                error_pb=parse_float(match.group("error")),
            )
        )
    if not rows:
        raise RuntimeError(f"Could not parse selected/nlo_real_tree rows from {path}")
    return rows


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


def load_gamma_helicity_measurements(
    *,
    base_dir: Path,
    tag: str,
    order: str,
    analysis: str,
    channel: str,
) -> Mapping[str, HistogramMeasurement]:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    manifest = load_manifest(campaign_dir)
    object_path = f"/{analysis}/Q2"
    out: dict[str, HistogramMeasurement] = {}
    for helicity in ("PP", "PM"):
        paths = [
            channel_output_path_for_entry(entry, channel, base_dir, campaign_dir)
            for entry in iter_manifest_entries(manifest, "GAMMA", order, helicity)
        ]
        out[helicity] = combine_shard_measurements(paths, object_path, channel=channel)
    return out


def rebin_density_to_integrals(
    measurement: HistogramMeasurement,
    target_edges: Sequence[float],
) -> HistogramMeasurement:
    values: list[float] = []
    errors: list[float] = []
    source_edges = measurement.edges

    for index in range(len(target_edges) - 1):
        lo = target_edges[index]
        hi = target_edges[index + 1]
        integral = 0.0
        err2 = 0.0

        for source_index, density in enumerate(measurement.values):
            source_lo = source_edges[source_index]
            source_hi = source_edges[source_index + 1]
            overlap = min(hi, source_hi) - max(lo, source_lo)
            if overlap <= 0.0:
                continue
            integral += density * overlap
            err2 += (measurement.errors[source_index] * overlap) ** 2

        values.append(integral)
        errors.append(math.sqrt(err2))

    return HistogramMeasurement(
        edges=tuple(target_edges),
        values=tuple(values),
        errors=tuple(errors),
    )


def load_raw_gamma_delta(
    *,
    base_dir: Path,
    tag: str,
    order: str,
    analysis: str,
    channel: str,
) -> HistogramMeasurement:
    measurements = load_gamma_helicity_measurements(
        base_dir=base_dir,
        tag=tag,
        order=order,
        analysis=analysis,
        channel=channel,
    )
    pp = measurements["PP"]
    pm = measurements["PM"]
    if pp.edges != pm.edges:
        raise RuntimeError("PP and PM histograms have different binning")
    values = []
    errors = []
    for pp_value, pp_error, pm_value, pm_error in zip(pp.values, pp.errors, pm.values, pm.errors):
        values.append(0.5 * (pp_value - pm_value))
        errors.append(0.5 * math.hypot(pp_error, pm_error))
    return HistogramMeasurement(edges=pp.edges, values=tuple(values), errors=tuple(errors))


def subtract_measurements(pos: HistogramMeasurement, neg: HistogramMeasurement) -> HistogramMeasurement:
    if pos.edges != neg.edges:
        raise RuntimeError("Cannot subtract measurements with different binning")
    values = []
    errors = []
    for pos_value, pos_error, neg_value, neg_error in zip(pos.values, pos.errors, neg.values, neg.errors):
        values.append(pos_value - neg_value)
        errors.append(math.hypot(pos_error, neg_error))
    return HistogramMeasurement(edges=pos.edges, values=tuple(values), errors=tuple(errors))


def fmt(value: float, precision: int) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.{precision}f}"


def ratio(num: float, den: float) -> float:
    if den == 0.0 or not math.isfinite(num) or not math.isfinite(den):
        return float("nan")
    return num / den


def raw_channel_label(channel: str) -> str:
    return "gluon/BGF" if channel == "BGF" else "quark/Compton"


def iter_shape_rows(
    edges: Sequence[float],
    values: Sequence[float],
    total: float,
    window_min: float,
    window_max: float,
) -> Iterable[tuple[float, float, float]]:
    for index, value in enumerate(values):
        lo = edges[index]
        hi = edges[index + 1]
        if lo < window_min or hi > window_max:
            continue
        yield lo, hi, (value / total if total != 0.0 else float("nan"))


def main() -> int:
    args = parse_args()
    base_dir = args.base_dir.resolve()
    poldis_rows = load_poldis_selected_real_tree(args.poldis_reference.resolve())
    target_edges = [poldis_rows[0].q2_low] + [row.q2_high for row in poldis_rows]
    herwig_channel = args.herwig_channel

    herwig_pos = load_raw_gamma_delta(
        base_dir=base_dir,
        tag=args.herwig_tag,
        order="POSNLO",
        analysis=args.analysis,
        channel=herwig_channel,
    )
    herwig_neg = load_raw_gamma_delta(
        base_dir=base_dir,
        tag=args.herwig_tag,
        order="NEGNLO",
        analysis=args.analysis,
        channel=herwig_channel,
    )

    herwig_pos_rebinned = rebin_density_to_integrals(herwig_pos, target_edges)
    herwig_neg_rebinned = rebin_density_to_integrals(herwig_neg, target_edges)
    herwig_nlo_rebinned = subtract_measurements(herwig_pos_rebinned, herwig_neg_rebinned)

    print(f"Herwig tag: {args.herwig_tag}")
    print(f"POLDIS reference: {args.poldis_reference.resolve()}")
    print("Comparison surface:")
    print(
        "  POLDIS selected/nlo_real_tree"
        f"  vs  raw Herwig GAMMA/POSNLO/{herwig_channel} ({raw_channel_label(herwig_channel)}) selected DQ2 rebinned as integrals"
    )
    print(f"  Also shown: raw Herwig GAMMA/NEGNLO/{herwig_channel} and POS-NEG on the same coarse bins")
    print()
    print(
        "q2 range"
        "             POLDIS_sel_real"
        "      Herwig_POS"
        "      Herwig_NEG"
        "      Herwig_NLO"
        "      POS/POLDIS"
        "         diff"
        "        pull"
    )

    for index, row in enumerate(poldis_rows):
        p_value = row.value_pb
        p_error = row.error_pb
        h_pos = herwig_pos_rebinned.values[index]
        h_neg = herwig_neg_rebinned.values[index]
        h_nlo = herwig_nlo_rebinned.values[index]
        h_error = herwig_pos_rebinned.errors[index]
        diff = h_pos - p_value
        denom = math.hypot(h_error, p_error)
        pull = diff / denom if denom > 0.0 else float("nan")
        print(
            f"[{row.q2_low:7.1f}, {row.q2_high:7.1f})"
            f"  {fmt(p_value, args.precision):>14}"
            f"  {fmt(h_pos, args.precision):>14}"
            f"  {fmt(h_neg, args.precision):>14}"
            f"  {fmt(h_nlo, args.precision):>14}"
            f"  {fmt(ratio(h_pos, p_value), args.precision):>14}"
            f"  {fmt(diff, args.precision):>12}"
            f"  {fmt(pull, args.precision):>11}"
        )

    poldis_total = math.fsum(row.value_pb for row in poldis_rows if row.q2_low >= args.shape_window_min and row.q2_high <= args.shape_window_max)
    herwig_pos_total = math.fsum(
        value
        for lo, hi, value in zip(herwig_pos_rebinned.edges[:-1], herwig_pos_rebinned.edges[1:], herwig_pos_rebinned.values)
        if lo >= args.shape_window_min and hi <= args.shape_window_max
    )
    herwig_nlo_total = math.fsum(
        value
        for lo, hi, value in zip(herwig_nlo_rebinned.edges[:-1], herwig_nlo_rebinned.edges[1:], herwig_nlo_rebinned.values)
        if lo >= args.shape_window_min and hi <= args.shape_window_max
    )

    print()
    print(
        f"Shape-only fractions in [{args.shape_window_min:.0f}, {args.shape_window_max:.0f}] GeV^2"
    )
    print("q2 range             POLDIS_frac     Herwig_POS_frac   Herwig_NLO_frac")
    poldis_shape = [(row.q2_low, row.q2_high, row.value_pb / poldis_total if poldis_total != 0.0 else float('nan')) for row in poldis_rows if row.q2_low >= args.shape_window_min and row.q2_high <= args.shape_window_max]
    herwig_pos_shape = list(
        iter_shape_rows(
            herwig_pos_rebinned.edges,
            herwig_pos_rebinned.values,
            herwig_pos_total,
            args.shape_window_min,
            args.shape_window_max,
        )
    )
    herwig_nlo_shape = list(
        iter_shape_rows(
            herwig_nlo_rebinned.edges,
            herwig_nlo_rebinned.values,
            herwig_nlo_total,
            args.shape_window_min,
            args.shape_window_max,
        )
    )
    for (lo, hi, pfrac), (_, _, hposfrac), (_, _, hnlofrac) in zip(poldis_shape, herwig_pos_shape, herwig_nlo_shape):
        print(
            f"[{lo:7.1f}, {hi:7.1f})"
            f"  {fmt(pfrac, args.precision):>14}"
            f"  {fmt(hposfrac, args.precision):>16}"
            f"  {fmt(hnlofrac, args.precision):>16}"
        )

    print()
    print("Window totals (integrals over the coarse bins):")
    print(f"  POLDIS selected/nlo_real_tree = {fmt(poldis_total, args.precision)} pb")
    print(f"  Herwig POSNLO/{herwig_channel:<8}       = {fmt(herwig_pos_total, args.precision)} pb")
    print(f"  Herwig POS-NEG/{herwig_channel:<8}      = {fmt(herwig_nlo_total, args.precision)} pb")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
