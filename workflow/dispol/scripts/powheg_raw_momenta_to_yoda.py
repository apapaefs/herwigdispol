#!/usr/bin/env python3.10
"""Convert raw POWHEG momentum diagnostics from Herwig logs into YODA."""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from poldis_top_to_yoda import new_binned_estimate, set_bin_val_err, write_yoda_gz


LINE_RE = re.compile(r"POWHEG_RAW_MOM\b(.*)")
SUMMARY_LINE_RE = re.compile(r"POWHEG_RAW_SUMMARY\b(.*)")

SUMMARY_FIELD_ORDER = [
    "event", "winner", "fallback", "Q2", "xB", "y", "s",
    "comp_status", "comp_trials", "comp_rejectXP", "comp_rejectVeto",
    "comp_weightNeg", "comp_weightHigh", "comp_xp", "comp_zp",
    "comp_xMapped", "comp_xT", "comp_xTMin", "comp_pT", "comp_phase",
    "comp_pdfRatio", "comp_alphaRatio", "comp_meAvg", "comp_wgt",
    "comp_pdfScale", "comp_alphaScale",
    "bgf_status", "bgf_trials", "bgf_rejectXP", "bgf_rejectVeto",
    "bgf_weightNeg", "bgf_weightHigh", "bgf_xp", "bgf_zp",
    "bgf_xMapped", "bgf_xT", "bgf_xTMin", "bgf_pT", "bgf_phase",
    "bgf_pdfRatio", "bgf_alphaRatio", "bgf_meAvg", "bgf_wgt",
    "bgf_pdfScale", "bgf_alphaScale",
]


@dataclass(frozen=True)
class FourVector:
    px: float
    py: float
    pz: float
    energy: float

    @property
    def pt(self) -> float:
        return math.hypot(self.px, self.py)

    @property
    def p(self) -> float:
        return math.sqrt(max(0.0, self.px * self.px + self.py * self.py + self.pz * self.pz))

    @property
    def mass2(self) -> float:
        return self.energy * self.energy - self.p * self.p

    @property
    def mass(self) -> float:
        return math.sqrt(max(0.0, self.mass2))

    @property
    def rapidity(self) -> float:
        denom = self.energy - self.pz
        num = self.energy + self.pz
        if denom <= 0.0 or num <= 0.0:
            return math.copysign(float("inf"), self.pz)
        return 0.5 * math.log(num / denom)

    @property
    def eta(self) -> float:
        momentum = self.p
        denom = momentum - self.pz
        num = momentum + self.pz
        if denom <= 0.0 or num <= 0.0:
            return math.copysign(float("inf"), self.pz)
        return 0.5 * math.log(num / denom)

    def __add__(self, other: "FourVector") -> "FourVector":
        return FourVector(
            self.px + other.px,
            self.py + other.py,
            self.pz + other.pz,
            self.energy + other.energy,
        )


@dataclass(frozen=True)
class RawPowhegEvent:
    event: Optional[int]
    proc: str
    contrib: int
    q2: float
    xb: float
    y: float
    ee: Optional[float]
    s: Optional[float]
    xp: Optional[float]
    zp: Optional[float]
    pt: float
    p1_breit: FourVector
    p2_breit: FourVector
    p1_lab: FourVector
    p2_lab: FourVector


RAW_CHANNELS = ("Compton", "BGF")


@dataclass
class CutflowCounter:
    parsed: int = 0
    pass_q2: int = 0
    pass_y: int = 0
    pass_pt: int = 0
    pass_rapidity: int = 0
    accepted: int = 0


def linear_edges(nbins: int, xmin: float, xmax: float) -> List[float]:
    width = (xmax - xmin) / nbins
    return [xmin + index * width for index in range(nbins + 1)]


def logspace_edges(nbins: int, xmin: float, xmax: float) -> List[float]:
    log_min = math.log10(xmin)
    log_max = math.log10(xmax)
    step = (log_max - log_min) / nbins
    return [10.0 ** (log_min + index * step) for index in range(nbins + 1)]


BINNING: Dict[str, List[float]] = {
    "Q2": linear_edges(100, 25.0, 2500.0),
    "Pt": linear_edges(15, 5.0, 30.0),
    "XBj": linear_edges(20, 0.0, 1.0),
    "Mjj": logspace_edges(15, 10.0, 100.0),
    "Eta": linear_edges(15, 0.0, 2.5),
    "Zeta": linear_edges(12, -1.75, -0.25),
    "pT1": linear_edges(15, 5.0, 30.0),
    "pT2": linear_edges(15, 5.0, 30.0),
    "pT2OverpT1": linear_edges(15, 0.0, 1.0),
    "pTAsym": linear_edges(15, 0.0, 1.0),
}

RAW_PT_MAP_EDGES = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 30.0]
RAW_Q2_MAP_EDGES = [25.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2500.0]
RAW_Y_MAP_EDGES = [0.2, 0.3, 0.4, 0.5, 0.6]
RAW_XP_MAP_EDGES = [0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0]


def parse_vector(fields: Dict[str, str], prefix: str) -> FourVector:
    return FourVector(
        float(fields[f"{prefix}x"]),
        float(fields[f"{prefix}y"]),
        float(fields[f"{prefix}z"]),
        float(fields[f"{prefix}E"]),
    )


def parse_key_value_fields(payload: str) -> Dict[str, str]:
    fields: Dict[str, str] = {}
    for token in payload.strip().split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        fields[key] = value
    return fields


def parse_raw_event(line: str) -> RawPowhegEvent | None:
    match = LINE_RE.search(line)
    if not match:
        return None
    fields = parse_key_value_fields(match.group(1))
    required = {
        "proc", "contrib", "Q2", "xB", "y", "pT",
        "p1Bx", "p1By", "p1Bz", "p1BE",
        "p2Bx", "p2By", "p2Bz", "p2BE",
        "p1Lx", "p1Ly", "p1Lz", "p1LE",
        "p2Lx", "p2Ly", "p2Lz", "p2LE",
    }
    if not required.issubset(fields):
        return None
    return RawPowhegEvent(
        event=int(fields["event"]) if "event" in fields else None,
        proc=fields["proc"],
        contrib=int(fields["contrib"]),
        q2=float(fields["Q2"]),
        xb=float(fields["xB"]),
        y=float(fields["y"]),
        ee=float(fields["Ee"]) if "Ee" in fields else None,
        s=float(fields["s"]) if "s" in fields else None,
        xp=float(fields["xp"]) if "xp" in fields else None,
        zp=float(fields["zp"]) if "zp" in fields else None,
        pt=float(fields["pT"]),
        p1_breit=parse_vector(fields, "p1B"),
        p2_breit=parse_vector(fields, "p2B"),
        p1_lab=parse_vector(fields, "p1L"),
        p2_lab=parse_vector(fields, "p2L"),
    )


def parse_raw_summary(line: str) -> Dict[str, str] | None:
    match = SUMMARY_LINE_RE.search(line)
    if not match:
        return None
    fields = parse_key_value_fields(match.group(1))
    if not set(SUMMARY_FIELD_ORDER).issubset(fields):
        return None
    return fields


def iter_raw_events(paths: Sequence[str]) -> Iterable[RawPowhegEvent]:
    for path in paths:
        with open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                event = parse_raw_event(line)
                if event is not None:
                    yield event


def iter_raw_summaries(paths: Sequence[str]) -> Iterable[Tuple[str, Dict[str, str]]]:
    for path in paths:
        with open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                summary = parse_raw_summary(line)
                if summary is not None:
                    yield str(Path(path)), summary


def write_summary_csv(output: str, rows: Sequence[Tuple[str, Dict[str, str]]]) -> None:
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*SUMMARY_FIELD_ORDER, "source_log"])
        writer.writeheader()
        for source_log, fields in rows:
            row = {name: fields.get(name, "") for name in SUMMARY_FIELD_ORDER}
            row["source_log"] = source_log
            writer.writerow(row)


def find_bin(edges: Sequence[float], value: float) -> int | None:
    if value < edges[0] or value >= edges[-1]:
        if value == edges[-1]:
            return len(edges) - 2
        return None
    for index in range(len(edges) - 1):
        if edges[index] <= value < edges[index + 1]:
            return index
    return None


def add_histogram_value(storage: Dict[str, List[float]], name: str, value: float, weight: float) -> None:
    edges = BINNING[name]
    index = find_bin(edges, value)
    if index is None:
        return
    storage[name][index] += weight


def new_cutflow_dict() -> Dict[str, int]:
    return {
        "parsed": 0,
        "pass_q2": 0,
        "pass_y": 0,
        "pass_pt": 0,
        "pass_rapidity": 0,
        "accepted": 0,
    }


def increment_cutflow(counter: Dict[str, int], stage: str) -> None:
    counter[stage] += 1


def new_pt_map() -> List[int]:
    return [0] * (len(RAW_PT_MAP_EDGES) - 1)


def new_xb_pt_map(xb_edges: Sequence[float]) -> List[List[int]]:
    return [[0] * (len(RAW_PT_MAP_EDGES) - 1) for _ in range(len(xb_edges) - 1)]


def new_xb_other_map(xb_edges: Sequence[float], other_edges: Sequence[float]) -> List[List[int]]:
    return [[0] * (len(other_edges) - 1) for _ in range(len(xb_edges) - 1)]


def format_edges(edges: Sequence[float], index: int) -> str:
    return f"[{edges[index]:.2f},{edges[index + 1]:.2f})"


def build_histograms(
    events: Sequence[RawPowhegEvent],
    scale: float,
) -> Tuple[
    Dict[str, List[float]],
    Dict[str, List[float]],
    int,
    CutflowCounter,
    Dict[str, CutflowCounter],
    List[Dict[str, int]],
    Dict[str, List[Dict[str, int]]],
    List[List[int]],
    Dict[str, List[List[int]]],
    List[List[int]],
    Dict[str, List[List[int]]],
    List[List[int]],
    Dict[str, List[List[int]]],
    List[List[int]],
    Dict[str, List[List[int]]],
    List[List[int]],
    Dict[str, List[List[int]]],
    Dict[str, float],
    Dict[str, Dict[str, float]],
    List[Dict[str, float]],
    Dict[str, List[Dict[str, float]]],
]:
    counts = {name: [0.0] * (len(edges) - 1) for name, edges in BINNING.items()}
    errs2 = {name: [0.0] * (len(edges) - 1) for name, edges in BINNING.items()}
    accepted = 0
    cutflow = CutflowCounter()
    channel_cutflows: Dict[str, CutflowCounter] = {}
    xb_edges = BINNING["XBj"]
    xb_cutflows = [new_cutflow_dict() for _ in range(len(xb_edges) - 1)]
    channel_xb_cutflows: Dict[str, List[Dict[str, int]]] = {
        channel: [new_cutflow_dict() for _ in range(len(xb_edges) - 1)]
        for channel in RAW_CHANNELS
    }
    xb_pt_counts = new_xb_pt_map(xb_edges)
    channel_xb_pt_counts: Dict[str, List[List[int]]] = {
        channel: new_xb_pt_map(xb_edges) for channel in RAW_CHANNELS
    }
    xb_q2_counts = new_xb_other_map(xb_edges, RAW_Q2_MAP_EDGES)
    channel_xb_q2_counts: Dict[str, List[List[int]]] = {
        channel: new_xb_other_map(xb_edges, RAW_Q2_MAP_EDGES) for channel in RAW_CHANNELS
    }
    xb_y_counts = new_xb_other_map(xb_edges, RAW_Y_MAP_EDGES)
    channel_xb_y_counts: Dict[str, List[List[int]]] = {
        channel: new_xb_other_map(xb_edges, RAW_Y_MAP_EDGES) for channel in RAW_CHANNELS
    }
    xb_xp_counts = new_xb_other_map(xb_edges, RAW_XP_MAP_EDGES)
    channel_xb_xp_counts: Dict[str, List[List[int]]] = {
        channel: new_xb_other_map(xb_edges, RAW_XP_MAP_EDGES) for channel in RAW_CHANNELS
    }
    xp_q2_counts = new_xb_other_map(RAW_XP_MAP_EDGES, RAW_Q2_MAP_EDGES)
    channel_xp_q2_counts: Dict[str, List[List[int]]] = {
        channel: new_xb_other_map(RAW_XP_MAP_EDGES, RAW_Q2_MAP_EDGES) for channel in RAW_CHANNELS
    }
    xb_consistency = {
        "count": 0.0,
        "sum_abs": 0.0,
        "sum_rel": 0.0,
        "max_abs": 0.0,
        "max_rel": 0.0,
    }
    channel_xb_consistency: Dict[str, Dict[str, float]] = {
        channel: {
            "count": 0.0,
            "sum_abs": 0.0,
            "sum_rel": 0.0,
            "max_abs": 0.0,
            "max_rel": 0.0,
        }
        for channel in RAW_CHANNELS
    }
    xb_means = [
        {
            "count": 0.0,
            "sum_pt": 0.0,
            "count_xp": 0.0,
            "sum_xp": 0.0,
            "count_zp": 0.0,
            "sum_zp": 0.0,
        }
        for _ in range(len(xb_edges) - 1)
    ]
    channel_xb_means: Dict[str, List[Dict[str, float]]] = {
        channel: [
            {
                "count": 0.0,
                "sum_pt": 0.0,
                "count_xp": 0.0,
                "sum_xp": 0.0,
                "count_zp": 0.0,
                "sum_zp": 0.0,
            }
            for _ in range(len(xb_edges) - 1)
        ]
        for channel in RAW_CHANNELS
    }

    for event in events:
        cutflow.parsed += 1
        chan = channel_cutflows.setdefault(event.proc, CutflowCounter())
        chan.parsed += 1
        xb_index = find_bin(xb_edges, event.xb)
        q2_index = find_bin(RAW_Q2_MAP_EDGES, event.q2)
        y_index = find_bin(RAW_Y_MAP_EDGES, event.y)
        xp_index = find_bin(RAW_XP_MAP_EDGES, event.xp) if event.xp is not None else None
        if xb_index is not None:
            increment_cutflow(xb_cutflows[xb_index], "parsed")
            if event.proc in channel_xb_cutflows:
                increment_cutflow(channel_xb_cutflows[event.proc][xb_index], "parsed")
            if q2_index is not None:
                xb_q2_counts[xb_index][q2_index] += 1
                if event.proc in channel_xb_q2_counts:
                    channel_xb_q2_counts[event.proc][xb_index][q2_index] += 1
            if y_index is not None:
                xb_y_counts[xb_index][y_index] += 1
                if event.proc in channel_xb_y_counts:
                    channel_xb_y_counts[event.proc][xb_index][y_index] += 1
            if xp_index is not None:
                xb_xp_counts[xb_index][xp_index] += 1
                if event.proc in channel_xb_xp_counts:
                    channel_xb_xp_counts[event.proc][xb_index][xp_index] += 1
            xb_means[xb_index]["count"] += 1.0
            xb_means[xb_index]["sum_pt"] += event.pt
            if event.xp is not None:
                xb_means[xb_index]["count_xp"] += 1.0
                xb_means[xb_index]["sum_xp"] += event.xp
            if event.zp is not None:
                xb_means[xb_index]["count_zp"] += 1.0
                xb_means[xb_index]["sum_zp"] += event.zp
            if event.proc in channel_xb_means:
                stats = channel_xb_means[event.proc][xb_index]
                stats["count"] += 1.0
                stats["sum_pt"] += event.pt
                if event.xp is not None:
                    stats["count_xp"] += 1.0
                    stats["sum_xp"] += event.xp
                if event.zp is not None:
                    stats["count_zp"] += 1.0
                    stats["sum_zp"] += event.zp
        if xp_index is not None and q2_index is not None:
            xp_q2_counts[xp_index][q2_index] += 1
            if event.proc in channel_xp_q2_counts:
                channel_xp_q2_counts[event.proc][xp_index][q2_index] += 1
        if event.s is not None and event.s > 0.0 and event.y > 0.0:
            xb_reco = event.q2 / (event.y * event.s)
            abs_diff = abs(xb_reco - event.xb)
            rel_diff = abs_diff / abs(event.xb) if event.xb != 0.0 else 0.0
            xb_consistency["count"] += 1.0
            xb_consistency["sum_abs"] += abs_diff
            xb_consistency["sum_rel"] += rel_diff
            xb_consistency["max_abs"] = max(xb_consistency["max_abs"], abs_diff)
            xb_consistency["max_rel"] = max(xb_consistency["max_rel"], rel_diff)
            if event.proc in channel_xb_consistency:
                stats = channel_xb_consistency[event.proc]
                stats["count"] += 1.0
                stats["sum_abs"] += abs_diff
                stats["sum_rel"] += rel_diff
                stats["max_abs"] = max(stats["max_abs"], abs_diff)
                stats["max_rel"] = max(stats["max_rel"], rel_diff)

        if not (25.0 <= event.q2 <= 2500.0):
            continue
        cutflow.pass_q2 += 1
        chan.pass_q2 += 1
        if xb_index is not None:
            increment_cutflow(xb_cutflows[xb_index], "pass_q2")
            if event.proc in channel_xb_cutflows:
                increment_cutflow(channel_xb_cutflows[event.proc][xb_index], "pass_q2")
        if not (0.2 <= event.y <= 0.6):
            continue
        cutflow.pass_y += 1
        chan.pass_y += 1
        if xb_index is not None:
            increment_cutflow(xb_cutflows[xb_index], "pass_y")
            if event.proc in channel_xb_cutflows:
                increment_cutflow(channel_xb_cutflows[event.proc][xb_index], "pass_y")

        hadrons = [
            (event.p1_breit, event.p1_lab),
            (event.p2_breit, event.p2_lab),
        ]
        hadrons.sort(key=lambda pair: pair[0].pt, reverse=True)
        (j1_breit, j1_lab), (j2_breit, j2_lab) = hadrons

        if j1_breit.pt < 5.0 or j2_breit.pt < 4.0:
            continue
        cutflow.pass_pt += 1
        chan.pass_pt += 1
        if xb_index is not None:
            increment_cutflow(xb_cutflows[xb_index], "pass_pt")
            pt_index = find_bin(RAW_PT_MAP_EDGES, event.pt)
            if pt_index is not None:
                xb_pt_counts[xb_index][pt_index] += 1
            if event.proc in channel_xb_cutflows:
                increment_cutflow(channel_xb_cutflows[event.proc][xb_index], "pass_pt")
                if pt_index is not None and event.proc in channel_xb_pt_counts:
                    channel_xb_pt_counts[event.proc][xb_index][pt_index] += 1
        if not (-3.5 <= j1_lab.rapidity <= 3.5):
            continue
        if not (-3.5 <= j2_lab.rapidity <= 3.5):
            continue
        cutflow.pass_rapidity += 1
        chan.pass_rapidity += 1
        if xb_index is not None:
            increment_cutflow(xb_cutflows[xb_index], "pass_rapidity")
            if event.proc in channel_xb_cutflows:
                increment_cutflow(channel_xb_cutflows[event.proc][xb_index], "pass_rapidity")

        accepted += 1
        cutflow.accepted += 1
        chan.accepted += 1
        if xb_index is not None:
            increment_cutflow(xb_cutflows[xb_index], "accepted")
            if event.proc in channel_xb_cutflows:
                increment_cutflow(channel_xb_cutflows[event.proc][xb_index], "accepted")
        weight = scale
        weight2 = scale * scale

        dijet_pt = 0.5 * (j1_breit.pt + j2_breit.pt)
        mjj = (j1_breit + j2_breit).mass
        eta_star = 0.5 * abs(j1_breit.eta - j2_breit.eta)
        zeta = math.log10(event.xb * (1.0 + (mjj * mjj) / event.q2))
        ratio = j2_breit.pt / j1_breit.pt
        asym = (j1_breit.pt - j2_breit.pt) / (j1_breit.pt + j2_breit.pt)

        values = {
            "Q2": event.q2,
            "Pt": dijet_pt,
            "XBj": event.xb,
            "Mjj": mjj,
            "Eta": eta_star,
            "Zeta": zeta,
            "pT1": j1_breit.pt,
            "pT2": j2_breit.pt,
            "pT2OverpT1": ratio,
            "pTAsym": asym,
        }

        for name, value in values.items():
            index = find_bin(BINNING[name], value)
            if index is None:
                continue
            counts[name][index] += weight
            errs2[name][index] += weight2

    return (
        counts,
        errs2,
        accepted,
        cutflow,
        channel_cutflows,
        xb_cutflows,
        channel_xb_cutflows,
        xb_pt_counts,
        channel_xb_pt_counts,
        xb_q2_counts,
        channel_xb_q2_counts,
        xb_y_counts,
        channel_xb_y_counts,
        xb_xp_counts,
        channel_xb_xp_counts,
        xp_q2_counts,
        channel_xp_q2_counts,
        xb_consistency,
        channel_xb_consistency,
        xb_means,
        channel_xb_means,
    )


def write_yoda(
    output: str,
    counts: Dict[str, List[float]],
    errs2: Dict[str, List[float]],
    analysis: str,
    legend: str,
) -> None:
    objects: Dict[str, object] = {}
    for name, edges in BINNING.items():
        path = f"/{analysis}/{name}"
        hist = new_binned_estimate(edges)
        hist.setPath(path)
        try:
            hist.setAnnotation("Legend", legend)
        except Exception:
            pass

        widths = [edges[index + 1] - edges[index] for index in range(len(edges) - 1)]
        densities = [
            (counts[name][index] / widths[index]) if widths[index] > 0.0 else 0.0
            for index in range(len(widths))
        ]
        errors = [
            (math.sqrt(errs2[name][index]) / widths[index]) if widths[index] > 0.0 else 0.0
            for index in range(len(widths))
        ]
        for index, binwrap in enumerate(hist.bins()):
            set_bin_val_err(binwrap, densities[index], errors[index])
        objects[path] = hist

    write_yoda_gz(objects, output)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", dest="logs", action="append", required=True,
                        help="Herwig log file containing POWHEG_RAW_MOM lines. May be given multiple times.")
    parser.add_argument("--output", required=True,
                        help="Output YODA or YODA.GZ file.")
    parser.add_argument("--analysis", default="MC_DIS_BREIT",
                        help="Analysis path to write inside YODA (default: MC_DIS_BREIT).")
    parser.add_argument("--legend", default="Raw POWHEG",
                        help="Legend annotation for the output objects.")
    parser.add_argument("--channel", choices=RAW_CHANNELS,
                        help="If set, keep only events from the selected POWHEG channel.")
    parser.add_argument("--cross-section-pb", type=float, default=None,
                        help="Total sample cross section in pb. If omitted, output is normalized per dumped event.")
    parser.add_argument("--sumw", type=float, default=None,
                        help="Total event sumW used for normalization. Defaults to the number of dumped events.")
    parser.add_argument("--summary-csv",
                        help="Optional CSV path for POWHEG_RAW_SUMMARY rows.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logs = [str(Path(path)) for path in args.logs]
    summary_rows = list(iter_raw_summaries(logs)) if args.summary_csv else []
    events = list(iter_raw_events(logs))
    if args.channel:
        events = [event for event in events if event.proc == args.channel]
    if not events:
        counts = {name: [0.0] * (len(edges) - 1) for name, edges in BINNING.items()}
        errs2 = {name: [0.0] * (len(edges) - 1) for name, edges in BINNING.items()}
        write_yoda(args.output, counts, errs2, args.analysis, args.legend)
        if args.summary_csv:
            write_summary_csv(args.summary_csv, summary_rows)
        print("[raw-powheg] parsed events: 0")
        print("[raw-powheg] accepted events after MC_DIS_BREIT cuts: 0")
        print("[raw-powheg] cutflow: parsed=0 pass_q2=0 pass_y=0 pass_pt=0 pass_rapidity=0 accepted=0")
        if args.channel:
            print(f"[raw-powheg] channel filter: {args.channel}")
        print("[raw-powheg] no POWHEG_RAW_MOM lines found; wrote zero-valued YODA")
        if args.summary_csv:
            print(f"[raw-powheg] wrote summary CSV {args.summary_csv}")
        print(f"[raw-powheg] wrote {args.output}")
        return 0

    sumw = args.sumw if args.sumw is not None else float(len(events))
    if sumw <= 0.0:
        raise SystemExit("sumW must be positive.")

    scale = 1.0 / sumw
    if args.cross_section_pb is not None:
        scale = args.cross_section_pb / sumw

    (
        counts,
        errs2,
        accepted,
        cutflow,
        channel_cutflows,
        xb_cutflows,
        channel_xb_cutflows,
        xb_pt_counts,
        channel_xb_pt_counts,
        xb_q2_counts,
        channel_xb_q2_counts,
        xb_y_counts,
        channel_xb_y_counts,
        xb_xp_counts,
        channel_xb_xp_counts,
        xp_q2_counts,
        channel_xp_q2_counts,
        xb_consistency,
        channel_xb_consistency,
        xb_means,
        channel_xb_means,
    ) = build_histograms(events, scale)
    write_yoda(args.output, counts, errs2, args.analysis, args.legend)
    if args.summary_csv:
        write_summary_csv(args.summary_csv, summary_rows)

    print(f"[raw-powheg] parsed events: {len(events)}")
    print(f"[raw-powheg] accepted events after MC_DIS_BREIT cuts: {accepted}")
    if args.channel:
        print(f"[raw-powheg] channel filter: {args.channel}")
    print(
        "[raw-powheg] cutflow:"
        f" parsed={cutflow.parsed}"
        f" pass_q2={cutflow.pass_q2}"
        f" pass_y={cutflow.pass_y}"
        f" pass_pt={cutflow.pass_pt}"
        f" pass_rapidity={cutflow.pass_rapidity}"
        f" accepted={cutflow.accepted}"
    )
    for proc in sorted(channel_cutflows):
        stats = channel_cutflows[proc]
        print(
            f"[raw-powheg] channel {proc}:"
            f" parsed={stats.parsed}"
            f" pass_q2={stats.pass_q2}"
            f" pass_y={stats.pass_y}"
            f" pass_pt={stats.pass_pt}"
            f" pass_rapidity={stats.pass_rapidity}"
            f" accepted={stats.accepted}"
        )
    xb_edges = BINNING["XBj"]
    print("[raw-powheg] xB-binned cutflow:")
    for index, stats in enumerate(xb_cutflows):
        if stats["parsed"] == 0:
            continue
        print(
            f"[raw-powheg]   xB {format_edges(xb_edges, index)}:"
            f" parsed={stats['parsed']}"
            f" pass_q2={stats['pass_q2']}"
            f" pass_y={stats['pass_y']}"
            f" pass_pt={stats['pass_pt']}"
            f" pass_rapidity={stats['pass_rapidity']}"
            f" accepted={stats['accepted']}"
        )
    for proc in RAW_CHANNELS:
        proc_cutflows = channel_xb_cutflows.get(proc, [])
        if not any(stats["parsed"] > 0 for stats in proc_cutflows):
            continue
        print(f"[raw-powheg] xB-binned cutflow ({proc}):")
        for index, stats in enumerate(proc_cutflows):
            if stats["parsed"] == 0:
                continue
            print(
                f"[raw-powheg]   xB {format_edges(xb_edges, index)}:"
                f" parsed={stats['parsed']}"
                f" pass_q2={stats['pass_q2']}"
                f" pass_y={stats['pass_y']}"
                f" pass_pt={stats['pass_pt']}"
                f" pass_rapidity={stats['pass_rapidity']}"
                f" accepted={stats['accepted']}"
            )
    print("[raw-powheg] xB vs raw-pT map (pass_pt counts, raw pair pT):")
    for xb_index, row in enumerate(xb_pt_counts):
        if not any(count > 0 for count in row):
            continue
        bins = ", ".join(
            f"{format_edges(RAW_PT_MAP_EDGES, pt_index)}={count}"
            for pt_index, count in enumerate(row)
            if count > 0
        )
        print(f"[raw-powheg]   xB {format_edges(xb_edges, xb_index)}: {bins}")
    for proc in RAW_CHANNELS:
        proc_counts = channel_xb_pt_counts.get(proc, [])
        if not any(any(count > 0 for count in row) for row in proc_counts):
            continue
        print(f"[raw-powheg] xB vs raw-pT map ({proc}, pass_pt counts):")
        for xb_index, row in enumerate(proc_counts):
            if not any(count > 0 for count in row):
                continue
            bins = ", ".join(
                f"{format_edges(RAW_PT_MAP_EDGES, pt_index)}={count}"
                for pt_index, count in enumerate(row)
                if count > 0
            )
            print(f"[raw-powheg]   xB {format_edges(xb_edges, xb_index)}: {bins}")
    print("[raw-powheg] xB vs Q2 map (parsed counts):")
    for xb_index, row in enumerate(xb_q2_counts):
        if not any(count > 0 for count in row):
            continue
        bins = ", ".join(
            f"{format_edges(RAW_Q2_MAP_EDGES, q2_index)}={count}"
            for q2_index, count in enumerate(row)
            if count > 0
        )
        print(f"[raw-powheg]   xB {format_edges(xb_edges, xb_index)}: {bins}")
    for proc in RAW_CHANNELS:
        proc_counts = channel_xb_q2_counts.get(proc, [])
        if not any(any(count > 0 for count in row) for row in proc_counts):
            continue
        print(f"[raw-powheg] xB vs Q2 map ({proc}, parsed counts):")
        for xb_index, row in enumerate(proc_counts):
            if not any(count > 0 for count in row):
                continue
            bins = ", ".join(
                f"{format_edges(RAW_Q2_MAP_EDGES, q2_index)}={count}"
                for q2_index, count in enumerate(row)
                if count > 0
            )
            print(f"[raw-powheg]   xB {format_edges(xb_edges, xb_index)}: {bins}")
    print("[raw-powheg] xB vs y map (parsed counts):")
    for xb_index, row in enumerate(xb_y_counts):
        if not any(count > 0 for count in row):
            continue
        bins = ", ".join(
            f"{format_edges(RAW_Y_MAP_EDGES, y_index)}={count}"
            for y_index, count in enumerate(row)
            if count > 0
        )
        print(f"[raw-powheg]   xB {format_edges(xb_edges, xb_index)}: {bins}")
    for proc in RAW_CHANNELS:
        proc_counts = channel_xb_y_counts.get(proc, [])
        if not any(any(count > 0 for count in row) for row in proc_counts):
            continue
        print(f"[raw-powheg] xB vs y map ({proc}, parsed counts):")
        for xb_index, row in enumerate(proc_counts):
            if not any(count > 0 for count in row):
                continue
            bins = ", ".join(
                f"{format_edges(RAW_Y_MAP_EDGES, y_index)}={count}"
                for y_index, count in enumerate(row)
                if count > 0
            )
            print(f"[raw-powheg]   xB {format_edges(xb_edges, xb_index)}: {bins}")
    print("[raw-powheg] xB vs xp map (parsed counts):")
    for xb_index, row in enumerate(xb_xp_counts):
        if not any(count > 0 for count in row):
            continue
        bins = ", ".join(
            f"{format_edges(RAW_XP_MAP_EDGES, xp_index)}={count}"
            for xp_index, count in enumerate(row)
            if count > 0
        )
        print(f"[raw-powheg]   xB {format_edges(xb_edges, xb_index)}: {bins}")
    for proc in RAW_CHANNELS:
        proc_counts = channel_xb_xp_counts.get(proc, [])
        if not any(any(count > 0 for count in row) for row in proc_counts):
            continue
        print(f"[raw-powheg] xB vs xp map ({proc}, parsed counts):")
        for xb_index, row in enumerate(proc_counts):
            if not any(count > 0 for count in row):
                continue
            bins = ", ".join(
                f"{format_edges(RAW_XP_MAP_EDGES, xp_index)}={count}"
                for xp_index, count in enumerate(row)
                if count > 0
            )
            print(f"[raw-powheg]   xB {format_edges(xb_edges, xb_index)}: {bins}")
    print("[raw-powheg] xp vs Q2 map (parsed counts):")
    for xp_index, row in enumerate(xp_q2_counts):
        if not any(count > 0 for count in row):
            continue
        bins = ", ".join(
            f"{format_edges(RAW_Q2_MAP_EDGES, q2_index)}={count}"
            for q2_index, count in enumerate(row)
            if count > 0
        )
        print(f"[raw-powheg]   xp {format_edges(RAW_XP_MAP_EDGES, xp_index)}: {bins}")
    for proc in RAW_CHANNELS:
        proc_counts = channel_xp_q2_counts.get(proc, [])
        if not any(any(count > 0 for count in row) for row in proc_counts):
            continue
        print(f"[raw-powheg] xp vs Q2 map ({proc}, parsed counts):")
        for xp_index, row in enumerate(proc_counts):
            if not any(count > 0 for count in row):
                continue
            bins = ", ".join(
                f"{format_edges(RAW_Q2_MAP_EDGES, q2_index)}={count}"
                for q2_index, count in enumerate(row)
                if count > 0
            )
            print(f"[raw-powheg]   xp {format_edges(RAW_XP_MAP_EDGES, xp_index)}: {bins}")
    print("[raw-powheg] xB-binned means (parsed events):")
    for index, stats in enumerate(xb_means):
        if stats["count"] <= 0.0:
            continue
        fields = [
            f"mean_raw_pT={stats['sum_pt']/stats['count']:.6f}",
        ]
        if stats["count_xp"] > 0.0:
            fields.append(f"mean_xp={stats['sum_xp']/stats['count_xp']:.6f}")
        if stats["count_zp"] > 0.0:
            fields.append(f"mean_zp={stats['sum_zp']/stats['count_zp']:.6f}")
        print(f"[raw-powheg]   xB {format_edges(xb_edges, index)}: " + " ".join(fields))
    for proc in RAW_CHANNELS:
        proc_means = channel_xb_means.get(proc, [])
        if not any(stats["count"] > 0.0 for stats in proc_means):
            continue
        print(f"[raw-powheg] xB-binned means ({proc}, parsed events):")
        for index, stats in enumerate(proc_means):
            if stats["count"] <= 0.0:
                continue
            fields = [
                f"mean_raw_pT={stats['sum_pt']/stats['count']:.6f}",
            ]
            if stats["count_xp"] > 0.0:
                fields.append(f"mean_xp={stats['sum_xp']/stats['count_xp']:.6f}")
            if stats["count_zp"] > 0.0:
                fields.append(f"mean_zp={stats['sum_zp']/stats['count_zp']:.6f}")
            print(f"[raw-powheg]   xB {format_edges(xb_edges, index)}: " + " ".join(fields))
    if xb_consistency["count"] > 0.0:
        mean_abs = xb_consistency["sum_abs"] / xb_consistency["count"]
        mean_rel = xb_consistency["sum_rel"] / xb_consistency["count"]
        print(
            "[raw-powheg] xB consistency:"
            f" count={int(xb_consistency['count'])}"
            f" mean_abs={mean_abs:.6e}"
            f" max_abs={xb_consistency['max_abs']:.6e}"
            f" mean_rel={mean_rel:.6e}"
            f" max_rel={xb_consistency['max_rel']:.6e}"
        )
        for proc in RAW_CHANNELS:
            stats = channel_xb_consistency.get(proc)
            if not stats or stats["count"] <= 0.0:
                continue
            mean_abs = stats["sum_abs"] / stats["count"]
            mean_rel = stats["sum_rel"] / stats["count"]
            print(
                f"[raw-powheg] xB consistency ({proc}):"
                f" count={int(stats['count'])}"
                f" mean_abs={mean_abs:.6e}"
                f" max_abs={stats['max_abs']:.6e}"
                f" mean_rel={mean_rel:.6e}"
                f" max_rel={stats['max_rel']:.6e}"
            )
    else:
        print("[raw-powheg] xB consistency: unavailable (missing s in POWHEG_RAW_MOM lines)")
    print(f"[raw-powheg] normalization scale per dumped event: {scale:.16e}")
    if args.summary_csv:
        print(f"[raw-powheg] wrote summary CSV {args.summary_csv}")
    print(f"[raw-powheg] wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
