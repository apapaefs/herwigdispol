#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
from contextlib import contextmanager
from dataclasses import dataclass
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time
from typing import Callable, Iterator, Sequence, TextIO

from herwig_fo.alphas import MatchboxAlphaS
from herwig_fo.born import DISPoint, born_event, born_weight, xb_from_q2_y
from herwig_fo.compare import compare_yoda
from herwig_fo.config import (
    ACTIVE_QUARK_FLAVORS,
    CUT_WINDOWS,
    DISCuts,
    RunConfig,
    cuts_for_window,
    HELICITIES,
    NC_SETUPS,
    PDF_PROFILES,
    default_run_config,
    require_helicity,
    require_pdf_profile,
    require_setup,
)
from herwig_fo.hepmc import HepMCEvent, write_hepmc3
from herwig_fo.nlo import herwig_xp_power_map, nlo_terms
from herwig_fo.pdfs import PDFPair, load_pdf_pair
from herwig_fo.real import local_real_weight_ratios, real_event
from herwig_fo.rivet import build_mc_dis_breit, run_rivet
from rivet_scale_plot_postprocess import (
    rewrite_no_scale_ratio_plot_scripts,
    rewrite_scale_envelope_plot_scripts,
)


RIVET_HELICITY_WEIGHT_NAMES = (
    "HERWIG_DIS_PP",
    "HERWIG_DIS_PM",
    "HERWIG_DIS_MP",
    "HERWIG_DIS_MM",
    "HERWIG_DIS_SIGMA",
    "HERWIG_DIS_DELTA_LL",
)
POLARIZED_REAL_WEIGHT_MODES = ("projected-nlo", "local-real")
DEFAULT_POLARIZED_REAL_WEIGHT_MODE = "local-real"
REAL_CONTRIBUTIONS = frozenset(("qcdc", "qcdc_counter", "bgf", "bgf_counter"))
DEFAULT_STANDALONE_POLDIS_REFERENCE = Path(
    "workflow/dispol/campaigns/rivetweights_noshowerMac03/"
    "analysis/_rivetplot_inputs/reference.scale.reference_all.sanitized.yoda.gz"
)
DEFAULT_STANDALONE_CAMPAIGN_MATCH = (
    "MC_DIS_BREIT/"
    "(Q2|XBj|Y|pT1|pT2|Mjj|Eta|Zeta|"
    "DQ2|DXBj|DY|DpT1|DpT2|DMjj|DEta|DZeta|"
    "ALLQ2|ALLXBj|ALLY|ALLpT1|ALLpT2|ALLMjj|ALLEta|ALLZeta|"
    "Q2PreCut|XBjPreCut|YPreCut|pT1PreCut|pT2PreCut|"
    "DQ2PreCut|DXBjPreCut|DYPreCut|DpT1PreCut|DpT2PreCut|"
    "ALLQ2PreCut|ALLXBjPreCut|ALLYPreCut|ALLpT1PreCut|ALLpT2PreCut)$"
)
SCALE_VARIATION_FACTORS = {
    "nominal": 1.0,
    "ScaleFactorDown": 0.5,
    "ScaleFactorUp": 2.0,
}


@dataclass(frozen=True)
class StandaloneShardSpec:
    index: int
    start: int
    events: int
    seed: int


@dataclass(frozen=True)
class StandaloneShardTask:
    variation_name: str
    scale_factor: float
    spec: StandaloneShardSpec
    work_dir: Path
    hepmc_path: Path
    raw_yoda_path: Path
    stdout_path: Path
    stderr_path: Path
    command: tuple[str, ...]


@dataclass(frozen=True)
class StandaloneShardResult:
    task: StandaloneShardTask
    returncode: int


def add_physics_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--setup", choices=NC_SETUPS, default="ALL", help="Neutral-current exchange setup.")
    parser.add_argument("--helicity", choices=tuple(HELICITIES), default="00", help="Beam helicity label.")
    parser.add_argument("--pdf-profile", choices=tuple(PDF_PROFILES), default="nnpdf_paired", help="PDF profile.")
    parser.add_argument(
        "--window",
        choices=tuple(CUT_WINDOWS),
        default="validation",
        help="Named DIS cut window; 'plain66' uses Q2min=49 GeV^2 from the reference campaign.",
    )
    parser.add_argument("--q2-min", type=float, default=None, help="Override the DIS lower Q2 cut in GeV^2.")
    parser.add_argument("--q2-max", type=float, default=None, help="Override the DIS upper Q2 cut in GeV^2.")
    parser.add_argument("--y-min", type=float, default=None, help="Override the DIS lower y cut.")
    parser.add_argument("--y-max", type=float, default=None, help="Override the DIS upper y cut.")
    parser.add_argument(
        "--q2-sampling",
        choices=("log", "linear"),
        default="log",
        help="Q2 integration variable for generated fixed-order seeds.",
    )
    parser.add_argument(
        "--sampler",
        choices=("halton", "random"),
        default="halton",
        help="Unit-cube sampler for generated fixed-order seeds.",
    )
    parser.add_argument("--toy-pdfs", action="store_true", help="Use deterministic toy PDFs for unit/smoke tests.")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed for the standalone integrator.")
    parser.add_argument(
        "--scale-factor",
        type=float,
        default=1.0,
        help="Common renormalization/factorization scale factor; mu2 = Q2 * factor^2.",
    )
    parser.add_argument("--jetinput", default="MEPARTONS", choices=("MEPARTONS",), help="Rivet jet input mode.")
    parser.add_argument("--no-progress", action="store_true", help="Suppress progress and elapsed-time messages.")
    parser.add_argument(
        "--powheg",
        action="store_true",
        help=argparse.SUPPRESS,
    )


def add_polarized_run_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--setup", choices=NC_SETUPS, default="ALL", help="Neutral-current exchange setup.")
    parser.add_argument("--pdf-profile", choices=tuple(PDF_PROFILES), default="nnpdf_paired", help="PDF profile.")
    parser.add_argument(
        "--window",
        choices=tuple(CUT_WINDOWS),
        default="validation",
        help="Named DIS cut window; 'plain66' uses Q2min=49 GeV^2 from the reference campaign.",
    )
    parser.add_argument("--q2-min", type=float, default=None, help="Override the DIS lower Q2 cut in GeV^2.")
    parser.add_argument("--q2-max", type=float, default=None, help="Override the DIS upper Q2 cut in GeV^2.")
    parser.add_argument("--y-min", type=float, default=None, help="Override the DIS lower y cut.")
    parser.add_argument("--y-max", type=float, default=None, help="Override the DIS upper y cut.")
    parser.add_argument(
        "--q2-sampling",
        choices=("log", "linear"),
        default="log",
        help="Q2 integration variable for generated fixed-order seeds.",
    )
    parser.add_argument(
        "--sampler",
        choices=("halton", "random"),
        default="halton",
        help="Unit-cube sampler for generated fixed-order seeds.",
    )
    parser.add_argument("--toy-pdfs", action="store_true", help="Use deterministic toy PDFs for unit/smoke tests.")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed for the standalone integrator.")
    parser.add_argument(
        "--scale-factor",
        type=float,
        default=1.0,
        help="Common renormalization/factorization scale factor; mu2 = Q2 * factor^2.",
    )
    parser.add_argument("--jetinput", default="MEPARTONS", choices=("MEPARTONS",), help="Rivet jet input mode.")
    parser.add_argument("--no-progress", action="store_true", help="Suppress progress and elapsed-time messages.")
    parser.add_argument(
        "--powheg",
        action="store_true",
        help=argparse.SUPPRESS,
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Standalone Python Herwig fixed-order NC DIS NLO validator.")
    sub = parser.add_subparsers(dest="command", required=True)

    generate = sub.add_parser("generate", help="Generate a signed HepMC3 fixed-order event stream.")
    add_physics_options(generate)
    generate.add_argument("--events", type=int, default=1000, help="Number of Monte Carlo phase-space seeds.")
    generate.add_argument("--output", type=Path, required=True, help="Output HepMC3 ASCII path.")

    integrate = sub.add_parser("integrate", help="Deterministically integrate standalone LO/NLO cross sections.")
    add_physics_options(integrate)
    integrate.add_argument("--order", choices=("LO", "NLO", "both"), default="NLO", help="Perturbative order to print.")
    integrate.add_argument(
        "--observable",
        choices=("helicity", "sigma0", "delta"),
        default="sigma0",
        help="Observable projection: selected helicity, unpolarized sigma0, or double-spin delta.",
    )
    integrate.add_argument("--q2-bins", type=int, default=80, help="Midpoint bins in the Q2 integration variable.")
    integrate.add_argument("--y-bins", type=int, default=80, help="Midpoint bins in y.")
    integrate.add_argument("--xp-bins", type=int, default=96, help="Midpoint bins in the Herwig xp random variable for NLO.")

    rivet = sub.add_parser("rivet", help="Run MC_DIS_BREIT:JETINPUT=MEPARTONS on an existing HepMC3 file.")
    rivet.add_argument("--hepmc", type=Path, required=True, help="Input HepMC3 ASCII path.")
    rivet.add_argument("--output", type=Path, required=True, help="Output YODA path.")
    rivet.add_argument("--analysis-plugin", type=Path, default=None, help="Existing MC_DIS_BREIT Rivet plugin.")
    rivet.add_argument("--build-analysis", action="store_true", help="Build MC_DIS_BREIT before running Rivet; this is the default when no plugin is supplied.")
    rivet.add_argument("--jetinput", default="MEPARTONS", choices=("MEPARTONS",), help="Rivet jet input mode.")
    rivet.add_argument("--rivetweights", action="store_true", help="Run MC_DIS_BREIT with RIVETWEIGHTS=YES.")
    rivet.add_argument("--no-progress", action="store_true", help="Suppress progress and elapsed-time messages.")

    compare = sub.add_parser("compare", help="Compare Herwig FO YODA output with POLDIS reference YODA.")
    compare.add_argument("--herwig", type=Path, required=True, help="Herwig FO YODA file.")
    compare.add_argument("--reference", type=Path, required=True, help="POLDIS reference YODA file.")
    compare.add_argument("--max-pull", type=float, default=3.0, help="Allowed combined-stat pull per bin.")

    polarize = sub.add_parser(
        "polarize",
        help="Build sigma, D*, and ALL* MC_DIS_BREIT YODA objects from helicity-resolved Rivet outputs.",
    )
    polarize.add_argument("--setup", choices=NC_SETUPS, default="ALL", help="Neutral-current exchange setup.")
    polarize.add_argument("--analysis", default="MC_DIS_BREIT", help="Clean output analysis path prefix.")
    polarize.add_argument("--00", dest="zero", type=Path, default=None, help="YODA file for the 00 helicity stream.")
    polarize.add_argument("--PP", type=Path, required=True, help="YODA file for the PP helicity stream.")
    polarize.add_argument("--PM", type=Path, required=True, help="YODA file for the PM helicity stream.")
    polarize.add_argument("--MP", type=Path, default=None, help="YODA file for the MP helicity stream.")
    polarize.add_argument("--MM", type=Path, default=None, help="YODA file for the MM helicity stream.")
    polarize.add_argument("--output", type=Path, required=True, help="Combined polarized YODA output path.")
    polarize.add_argument("--no-progress", action="store_true", help="Suppress progress and elapsed-time messages.")

    polarized_full = sub.add_parser(
        "polarized-full",
        help="Generate all needed helicity streams, run Rivet, and build sigma/D*/ALL* YODA output.",
    )
    add_polarized_run_options(polarized_full)
    polarized_full.add_argument("--events", type=int, default=1000, help="Number of Monte Carlo phase-space seeds per helicity.")
    polarized_full.add_argument("--work-dir", type=Path, required=True, help="Directory for per-helicity HepMC/YODA products.")
    polarized_full.add_argument("--output", type=Path, required=True, help="Combined polarized YODA output path.")
    polarized_full.add_argument(
        "--variation-output",
        type=Path,
        default=None,
        help="Variation-tagged YODA output path when --scale-variations is used.",
    )
    polarized_full.add_argument("--analysis", default="MC_DIS_BREIT", help="Clean output analysis path prefix.")
    polarized_full.add_argument("--analysis-plugin", type=Path, default=None, help="Existing MC_DIS_BREIT Rivet plugin.")
    polarized_full.add_argument("--build-analysis", action="store_true", help="Build MC_DIS_BREIT before running Rivet; this is the default when no plugin is supplied.")
    polarized_full.add_argument(
        "--scale-variations",
        action="store_true",
        help="Also run ScaleFactorDown/Up and write a variation-tagged YODA file.",
    )
    polarized_full.add_argument(
        "--independent-helicities",
        action="store_true",
        help="Diagnostic path: run separate PP/PM/MP/MM streams instead of correlated named weights.",
    )
    polarized_full.add_argument(
        "--polarized-real-weight-mode",
        choices=POLARIZED_REAL_WEIGHT_MODES,
        default=DEFAULT_POLARIZED_REAL_WEIGHT_MODE,
        help=(
            "How correlated polarized named weights are attached to real/counterevent pieces. "
            "'local-real' uses the fixed-order local real-emission spin kernels; "
            "'projected-nlo' is the Herwig RIVETWEIGHTS-style Born-projected NLO helicity-ratio diagnostic."
        ),
    )

    polarized_shard = sub.add_parser(
        "polarized-shard",
        help=argparse.SUPPRESS,
    )
    add_polarized_run_options(polarized_shard)
    polarized_shard.add_argument("--events", type=int, required=True, help=argparse.SUPPRESS)
    polarized_shard.add_argument("--work-dir", type=Path, required=True, help=argparse.SUPPRESS)
    polarized_shard.add_argument("--hepmc", type=Path, required=True, help=argparse.SUPPRESS)
    polarized_shard.add_argument("--raw-yoda", type=Path, required=True, help=argparse.SUPPRESS)
    polarized_shard.add_argument("--analysis", default="MC_DIS_BREIT", help=argparse.SUPPRESS)
    polarized_shard.add_argument("--analysis-plugin", type=Path, default=None, help=argparse.SUPPRESS)
    polarized_shard.add_argument("--build-analysis", action="store_true", help=argparse.SUPPRESS)
    polarized_shard.add_argument("--variation-name", default="nominal", help=argparse.SUPPRESS)
    polarized_shard.add_argument("--delete-hepmc-after-rivet", action="store_true", help=argparse.SUPPRESS)
    polarized_shard.add_argument(
        "--polarized-real-weight-mode",
        choices=POLARIZED_REAL_WEIGHT_MODES,
        default=DEFAULT_POLARIZED_REAL_WEIGHT_MODE,
        help=argparse.SUPPRESS,
    )

    plot_prep = sub.add_parser("plot-prep", help="Create plot-safe symmetric-error YODAs for standalone/POLDIS overlays.")
    plot_prep.add_argument("--herwig", type=Path, required=True, help="Standalone Herwig YODA input.")
    plot_prep.add_argument("--poldis", type=Path, required=True, help="POLDIS reference YODA input.")
    plot_prep.add_argument("--herwig-output", type=Path, required=True, help="Plot-safe Herwig YODA output.")
    plot_prep.add_argument("--poldis-output", type=Path, required=True, help="Plot-safe POLDIS YODA output.")
    plot_prep.add_argument("--analysis", default="MC_DIS_BREIT", help="Analysis path prefix.")
    plot_prep.add_argument(
        "--herwig-scale-envelope",
        action="store_true",
        help="Treat the Herwig input as a variation-tagged YODA and build the scale-envelope plot object first.",
    )

    full = sub.add_parser("full", help="Generate HepMC3, run Rivet, and compare to a POLDIS YODA reference.")
    add_physics_options(full)
    full.add_argument("--events", type=int, default=1000, help="Number of Monte Carlo phase-space seeds.")
    full.add_argument("--hepmc", type=Path, required=True, help="Generated HepMC3 path.")
    full.add_argument("--yoda", type=Path, required=True, help="Generated Herwig FO YODA path.")
    full.add_argument("--reference", type=Path, required=True, help="POLDIS reference YODA file.")
    full.add_argument("--analysis-plugin", type=Path, default=None, help="Existing MC_DIS_BREIT Rivet plugin.")
    full.add_argument("--build-analysis", action="store_true", help="Build MC_DIS_BREIT before running Rivet; this is the default when no plugin is supplied.")
    full.add_argument("--max-pull", type=float, default=3.0, help="Allowed combined-stat pull per bin.")

    standalone_campaign = sub.add_parser(
        "standalone-campaign",
        help="Run standalone correlated FO generation, plot-prep, and POLDIS overlay plots in one command.",
    )
    add_polarized_run_options(standalone_campaign)
    standalone_campaign.add_argument("-t", "--tag", required=True, help="Campaign tag under --base-dir/campaigns.")
    standalone_campaign.add_argument("--base-dir", type=Path, default=Path("workflow/dispol"), help="DISPOL workflow base directory.")
    standalone_campaign.add_argument("--campaign-dir", type=Path, default=None, help="Override the full campaign output directory.")
    standalone_campaign.add_argument("--events", type=int, default=20000, help="Number of Monte Carlo phase-space seeds.")
    standalone_campaign.add_argument("--jobs", type=int, default=1, help="Parallel shard jobs for standalone FO generation.")
    standalone_campaign.add_argument("--shards", type=int, default=1, help="Number of generation/Rivet shards to split --events over.")
    standalone_campaign.add_argument(
        "--keep-shard-hepmc",
        action="store_true",
        help="Keep per-shard HepMC3 files. By default parallel shards delete them after Rivet succeeds.",
    )
    standalone_campaign.add_argument(
        "--yoda-merge-tool",
        default="auto",
        help="YODA merge tool for parallel shards: auto, yodamerge, rivet-merge, or an explicit executable.",
    )
    standalone_campaign.add_argument(
        "--shard-monitor-interval",
        type=float,
        default=30.0,
        help="Seconds between live shard monitor updates during parallel standalone generation.",
    )
    standalone_campaign.add_argument(
        "--shard-status-file",
        type=Path,
        default=None,
        help="TSV status file for parallel standalone shards; defaults to campaign work/shard_status.tsv.",
    )
    standalone_campaign.add_argument(
        "--poldis-reference",
        type=Path,
        default=DEFAULT_STANDALONE_POLDIS_REFERENCE,
        help="POLDIS reference YODA/YODA.GZ for overlays.",
    )
    standalone_campaign.add_argument("--analysis", default="MC_DIS_BREIT", help="Clean output analysis path prefix.")
    standalone_campaign.add_argument("--analysis-plugin", type=Path, default=None, help="Existing MC_DIS_BREIT Rivet plugin.")
    standalone_campaign.add_argument("--build-analysis", action="store_true", help="Build MC_DIS_BREIT before running Rivet; this is the default when no plugin is supplied.")
    standalone_campaign.add_argument(
        "--polarized-real-weight-mode",
        choices=POLARIZED_REAL_WEIGHT_MODES,
        default=DEFAULT_POLARIZED_REAL_WEIGHT_MODE,
        help="Correlated polarized real/counterevent named-weight mode.",
    )
    standalone_campaign.add_argument(
        "--no-scale-variations",
        dest="scale_variations",
        action="store_false",
        default=True,
        help="Disable nominal/ScaleFactorDown/ScaleFactorUp generation and plot the nominal only.",
    )
    standalone_campaign.add_argument("--rivet-mkhtml", type=Path, default=None, help="Path to rivet-mkhtml; defaults to PATH lookup.")
    standalone_campaign.add_argument("--rivet-mkhtml-safe", type=Path, default=None, help="Path to rivet_mkhtml_safe.py.")
    standalone_campaign.add_argument("--analysis-rivet-dir", type=Path, default=None, help="Directory containing MC_DIS_BREIT.plot and Rivet analysis sources.")
    standalone_campaign.add_argument("--plot-config", type=Path, default=None, help="Rivet plot config file.")
    standalone_campaign.add_argument("--plots-dir", type=Path, default=None, help="Output plot directory.")
    standalone_campaign.add_argument("--output-stem", default=None, help="Stem for standalone YODA products.")
    standalone_campaign.add_argument("--poldis-output-stem", default=None, help="Stem for plot-safe POLDIS YODA.")
    standalone_campaign.add_argument("--match", default=DEFAULT_STANDALONE_CAMPAIGN_MATCH, help="rivet-mkhtml --match expression.")
    standalone_campaign.add_argument("--herwig-title", default="Standalone Python Herwig FO", help="Title for the standalone curve.")
    standalone_campaign.add_argument("--poldis-title", default="POLDIS NLO", help="Title for the POLDIS curve.")
    standalone_campaign.add_argument("--reflabel", default="POLDIS NLO", help="rivet-mkhtml reference label.")
    standalone_campaign.add_argument("--ratiolabel", default="Python FO/POLDIS", help="rivet-mkhtml ratio label.")
    standalone_campaign.add_argument("--mplconfigdir", type=Path, default=Path("/private/tmp/mplconfig"), help="MPLCONFIGDIR for rivet-mkhtml.")
    standalone_campaign.add_argument("--quiet-mkhtml", action="store_true", help="Do not pass --verbose to rivet-mkhtml.")
    return parser.parse_args(argv)


def required_polarized_helicities(setup: str) -> tuple[str, ...]:
    setup = require_setup(setup)
    if setup == "GAMMA":
        return ("00", "PP", "PM")
    return ("00", "PP", "PM", "MP", "MM")


def _format_duration(seconds: float) -> str:
    rounded = max(0, int(round(seconds)))
    if rounded < 60:
        return f"{rounded}s"
    minutes, sec = divmod(rounded, 60)
    if minutes < 60:
        return f"{minutes}m{sec:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h{minutes:02d}m{sec:02d}s"


class ProgressReporter:
    def __init__(
        self,
        label: str,
        total: int,
        unit: str = "items",
        stream: TextIO | None = None,
        clock: Callable[[], float] = time.monotonic,
        min_interval: float = 10.0,
        fraction_interval: float = 0.05,
    ) -> None:
        self.label = label
        self.total = max(0, int(total))
        self.unit = unit
        self.stream = stream if stream is not None else sys.stderr
        self.clock = clock
        self.min_interval = min_interval
        self.fraction_interval = fraction_interval
        self._start_time: float | None = None
        self._last_time: float | None = None
        self._next_fraction = fraction_interval
        self._last_completed = 0

    def _write(self, message: str) -> None:
        if self.stream is None:
            return
        print(message, file=self.stream, flush=True)

    def start(self) -> None:
        now = self.clock()
        self._start_time = now
        self._last_time = now
        self._last_completed = 0
        self._next_fraction = self.fraction_interval
        self._write(f"[{self.label}] start: {self.total} {self.unit}")

    def update(self, completed: int, force: bool = False) -> None:
        if self.total <= 0:
            return
        if self._start_time is None:
            self.start()
        now = self.clock()
        completed = min(max(0, int(completed)), self.total)
        fraction = completed / self.total
        time_due = self._last_time is None or (now - self._last_time) >= self.min_interval
        fraction_due = fraction >= self._next_fraction or completed == self.total
        if not force and not time_due and not fraction_due:
            return
        start_time = self._start_time if self._start_time is not None else now
        elapsed = max(now - start_time, 0.0)
        rate = completed / elapsed if elapsed > 0.0 else 0.0
        remaining = self.total - completed
        eta = remaining / rate if rate > 0.0 else math.inf
        eta_text = _format_duration(eta) if math.isfinite(eta) else "unknown"
        self._write(
            f"[{self.label}] {100.0 * fraction:5.1f}% "
            f"({completed}/{self.total} {self.unit}), "
            f"{rate:.2f} {self.unit}/s, ETA {eta_text}"
        )
        self._last_time = now
        self._last_completed = completed
        while self._next_fraction <= fraction:
            self._next_fraction += self.fraction_interval

    def done(self) -> None:
        if self._start_time is None:
            self.start()
        now = self.clock()
        start_time = self._start_time if self._start_time is not None else now
        elapsed = max(now - start_time, 0.0)
        rate = self.total / elapsed if elapsed > 0.0 else 0.0
        self._write(
            f"[{self.label}] done: {self.total} {self.unit} "
            f"in {_format_duration(elapsed)} ({rate:.2f} {self.unit}/s)"
        )


def _progress_stream(args: argparse.Namespace) -> TextIO | None:
    return None if bool(getattr(args, "no_progress", False)) else sys.stderr


def _progress_message(args: argparse.Namespace, message: str) -> None:
    stream = _progress_stream(args)
    if stream is not None:
        print(message, file=stream, flush=True)


@contextmanager
def _timed_step(args: argparse.Namespace, label: str) -> Iterator[None]:
    stream = _progress_stream(args)
    start = time.monotonic()
    if stream is not None:
        print(f"[{label}] start", file=stream, flush=True)
    try:
        yield
    finally:
        if stream is not None:
            print(f"[{label}] done in {_format_duration(time.monotonic() - start)}", file=stream, flush=True)


def _numbered_step(args: argparse.Namespace, index: int, total: int, label: str) -> None:
    command = getattr(args, "command", "step")
    _progress_message(args, f"[{command}] step {index}/{total}: {label}")


def _standalone_shard_specs(total_events: int, shard_count: int, seed: int) -> list[StandaloneShardSpec]:
    if total_events <= 0:
        raise SystemExit("--events must be positive.")
    if shard_count <= 0:
        raise SystemExit("--shards must be positive.")
    effective_shards = min(int(shard_count), int(total_events))
    base, remainder = divmod(int(total_events), effective_shards)
    specs: list[StandaloneShardSpec] = []
    start = 0
    for index in range(effective_shards):
        events = base + (1 if index < remainder else 0)
        specs.append(StandaloneShardSpec(index=index, start=start, events=events, seed=int(seed) + start))
        start += events
    return specs


def _parallel_generation_requested(args: argparse.Namespace) -> bool:
    return int(getattr(args, "jobs", 1)) > 1 or int(getattr(args, "shards", 1)) > 1


def _standalone_shard_key(task: StandaloneShardTask) -> tuple[str, int]:
    return (task.variation_name, task.spec.index)


def _write_shard_status_file(
    status_path: Path,
    tasks: Sequence[StandaloneShardTask],
    states: dict[tuple[str, int], str],
    returncodes: dict[tuple[str, int], int],
    start_times: dict[tuple[str, int], float],
    end_times: dict[tuple[str, int], float],
    *,
    now: float | None = None,
) -> None:
    timestamp = time.monotonic() if now is None else float(now)
    status_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["variation\tshard\tevents\tseed\tstatus\treturncode\telapsed_s\traw_yoda\tstdout\tstderr"]
    for task in tasks:
        key = _standalone_shard_key(task)
        status = states.get(key, "pending")
        returncode = returncodes.get(key)
        elapsed = ""
        if key in start_times:
            end = end_times.get(key, timestamp)
            elapsed = f"{max(0.0, end - start_times[key]):.3f}"
        lines.append(
            "\t".join(
                (
                    task.variation_name,
                    str(task.spec.index),
                    str(task.spec.events),
                    str(task.spec.seed),
                    status,
                    "" if returncode is None else str(returncode),
                    elapsed,
                    str(task.raw_yoda_path),
                    str(task.stdout_path),
                    str(task.stderr_path),
                )
            )
        )
    status_path.write_text("\n".join(lines) + "\n")


def _shard_monitor_message(
    *,
    total: int,
    finished: int,
    ok: int,
    failed: int,
    running: int,
    pending: int,
    started_at: float,
    now: float,
    status_path: Path,
) -> str:
    elapsed = max(0.0, now - started_at)
    rate = finished / elapsed if elapsed > 0.0 else 0.0
    eta = (total - finished) / rate if rate > 0.0 else math.inf
    eta_text = _format_duration(eta) if math.isfinite(eta) else "unknown"
    return (
        f"[standalone-campaign] shards: finished {finished}/{total} | ok {ok} | failed {failed} | "
        f"running {running} | pending {pending} | elapsed {_format_duration(elapsed)} | "
        f"ETA {eta_text} | status {status_path}"
    )


def run_config_from_args(args: argparse.Namespace) -> RunConfig:
    base = default_run_config()
    cuts = cuts_for_window(getattr(args, "window", "validation"))
    q2_min = cuts.q2_min if getattr(args, "q2_min", None) is None else float(args.q2_min)
    q2_max = cuts.q2_max if getattr(args, "q2_max", None) is None else float(args.q2_max)
    y_min = cuts.y_min if getattr(args, "y_min", None) is None else float(args.y_min)
    y_max = cuts.y_max if getattr(args, "y_max", None) is None else float(args.y_max)
    if not (0.0 < q2_min < q2_max):
        raise SystemExit("Require 0 < q2-min < q2-max.")
    if not (0.0 < y_min < y_max < 1.0):
        raise SystemExit("Require 0 < y-min < y-max < 1.")
    if float(getattr(args, "scale_factor", 1.0)) <= 0.0:
        raise SystemExit("Require --scale-factor > 0.")
    return RunConfig(
        beams=base.beams,
        cuts=DISCuts(q2_min=q2_min, q2_max=q2_max, y_min=y_min, y_max=y_max),
        pdf_profile=base.pdf_profile,
        setups=base.setups,
        helicities=base.helicities,
        alpha_s_mz=base.alpha_s_mz,
        mz=base.mz,
    )


def _q2_volume_and_jacobian(q2: float, cuts: DISCuts, sampling: str) -> tuple[float, float]:
    if sampling == "linear":
        return cuts.q2_max - cuts.q2_min, 1.0
    if sampling == "log":
        return math.log(cuts.q2_max / cuts.q2_min), q2
    raise ValueError(f"Unsupported Q2 sampling '{sampling}'.")


def _van_der_corput(index: int, base: int) -> float:
    if index <= 0:
        raise ValueError("Van der Corput indices are one-based in this sampler.")
    inverse = 1.0 / base
    factor = inverse
    value = 0.0
    while index:
        index, digit = divmod(index, base)
        value += digit * factor
        factor *= inverse
    return min(max(value, 1.0e-15), 1.0 - 1.0e-15)


def _seed_unit_coordinates(index: int, rng: random.Random, sampler: str, seed: int) -> tuple[float, float, float, float, float]:
    if sampler == "random":
        return (rng.random(), rng.random(), rng.random(), rng.random(), rng.random())
    if sampler == "halton":
        halton_index = max(1, seed) + index + 1
        return (
            _van_der_corput(halton_index, 2),
            _van_der_corput(halton_index, 3),
            _van_der_corput(halton_index, 5),
            _van_der_corput(halton_index, 7),
            _van_der_corput(halton_index, 11),
        )
    raise ValueError(f"Unsupported sampler '{sampler}'.")


def _coordinates_from_unit(q2_u: float, y_u: float, run: RunConfig, q2_sampling: str) -> tuple[float, float, float, float]:
    q2_volume, _ = _q2_volume_and_jacobian(run.cuts.q2_min, run.cuts, q2_sampling)
    if q2_sampling == "linear":
        q2 = run.cuts.q2_min + q2_u * q2_volume
    elif q2_sampling == "log":
        q2 = math.exp(math.log(run.cuts.q2_min) + q2_u * q2_volume)
    else:
        raise ValueError(f"Unsupported Q2 sampling '{q2_sampling}'.")
    y = run.cuts.y_min + y_u * (run.cuts.y_max - run.cuts.y_min)
    x_b = xb_from_q2_y(q2, y, run.beams)
    _, q2_jacobian = _q2_volume_and_jacobian(q2, run.cuts, q2_sampling)
    return q2, y, x_b, q2_jacobian


def _sample_coordinates(rng: random.Random, run: RunConfig, q2_sampling: str = "log") -> tuple[float, float, float, float]:
    return _coordinates_from_unit(rng.random(), rng.random(), run, q2_sampling)


def _point_from_coordinates(
    q2: float,
    y: float,
    x_b: float,
    flavor: int,
    setup: str,
    helicity: str,
    run: RunConfig | None = None,
    scale_factor: float = 1.0,
) -> DISPoint:
    if run is None:
        run = default_run_config()
    lepton_pol, hadron_pol = require_helicity(helicity)
    mu2 = q2 * scale_factor * scale_factor
    return DISPoint(
        q2=q2,
        y=y,
        x_b=x_b,
        flavor=flavor,
        setup=require_setup(setup),
        lepton_pol=lepton_pol,
        hadron_pol=hadron_pol,
        lepton_id=run.beams.lepton_id,
        mu2=mu2,
    )


def _sample_point(rng: random.Random, setup: str, helicity: str, run: RunConfig | None = None) -> DISPoint:
    if run is None:
        run = default_run_config()
    q2, y, x_b, _ = _sample_coordinates(rng, run)
    return _point_from_coordinates(q2, y, x_b, rng.choice(ACTIVE_QUARK_FLAVORS), setup, helicity, run)


@dataclass(frozen=True)
class HelicityWeightBundle:
    nominal: float
    pp: float
    pm: float
    mp: float
    mm: float
    sigma: float
    delta_ll: float

    def named_multipliers(self) -> dict[str, float]:
        if abs(self.nominal) <= 1.0e-300 or not math.isfinite(self.nominal):
            return {name: 0.0 for name in RIVET_HELICITY_WEIGHT_NAMES}
        return {
            "HERWIG_DIS_PP": self.pp / self.nominal,
            "HERWIG_DIS_PM": self.pm / self.nominal,
            "HERWIG_DIS_MP": self.mp / self.nominal,
            "HERWIG_DIS_MM": self.mm / self.nominal,
            "HERWIG_DIS_SIGMA": self.sigma / self.nominal,
            "HERWIG_DIS_DELTA_LL": self.delta_ll / self.nominal,
        }


def single_helicity_contribution_weight(
    helicity: str,
    q2: float,
    y: float,
    x_b: float,
    flavor: int,
    setup: str,
    contribution: str,
    base: float,
    x_p: float,
    z_p: float,
    phi: float,
    pdfs: PDFPair,
    alpha_s: Callable[[float], float],
    run: RunConfig | None = None,
    scale_factor: float = 1.0,
    x_p_jacobian: float = 1.0,
) -> float:
    if run is None:
        run = default_run_config()
    point = _point_from_coordinates(q2, y, x_b, flavor, setup, helicity, run, scale_factor=scale_factor)
    base_weight = base * born_weight(point, pdfs)
    contribution_key = contribution.lower()
    if contribution_key == "inclusive":
        return base_weight * nlo_terms(point, x_p=x_p, pdfs=pdfs, alpha_s=alpha_s, jacobian=x_p_jacobian).total
    local_real = local_real_weight_ratios(point, x_p, z_p, phi, pdfs, alpha_s)
    if contribution_key == "qcdc":
        return base_weight * local_real.qcdc
    if contribution_key == "qcdc_counter":
        return -base_weight * local_real.qcdc
    if contribution_key == "bgf":
        return base_weight * local_real.bgf
    if contribution_key == "bgf_counter":
        return -base_weight * local_real.bgf
    raise ValueError(f"Unsupported fixed-order contribution '{contribution}'.")


def correlated_helicity_weights(
    q2: float,
    y: float,
    x_b: float,
    flavor: int,
    setup: str,
    contribution: str,
    base: float,
    x_p: float,
    z_p: float,
    phi: float,
    pdfs: PDFPair,
    alpha_s: Callable[[float], float],
    run: RunConfig | None = None,
    scale_factor: float = 1.0,
    x_p_jacobian: float = 1.0,
) -> HelicityWeightBundle:
    pp = single_helicity_contribution_weight(
        "PP", q2, y, x_b, flavor, setup, contribution, base, x_p, z_p, phi, pdfs, alpha_s, run=run, scale_factor=scale_factor, x_p_jacobian=x_p_jacobian
    )
    pm = single_helicity_contribution_weight(
        "PM", q2, y, x_b, flavor, setup, contribution, base, x_p, z_p, phi, pdfs, alpha_s, run=run, scale_factor=scale_factor, x_p_jacobian=x_p_jacobian
    )
    mp = single_helicity_contribution_weight(
        "MP", q2, y, x_b, flavor, setup, contribution, base, x_p, z_p, phi, pdfs, alpha_s, run=run, scale_factor=scale_factor, x_p_jacobian=x_p_jacobian
    )
    mm = single_helicity_contribution_weight(
        "MM", q2, y, x_b, flavor, setup, contribution, base, x_p, z_p, phi, pdfs, alpha_s, run=run, scale_factor=scale_factor, x_p_jacobian=x_p_jacobian
    )
    sigma = 0.25 * math.fsum((pp, pm, mp, mm))
    delta_ll = 0.25 * math.fsum((pp, mm, -pm, -mp))
    return HelicityWeightBundle(
        nominal=sigma,
        pp=pp,
        pm=pm,
        mp=mp,
        mm=mm,
        sigma=sigma,
        delta_ll=delta_ll,
    )


def projected_nlo_correlated_helicity_weights(
    q2: float,
    y: float,
    x_b: float,
    flavor: int,
    setup: str,
    contribution: str,
    base: float,
    x_p: float,
    z_p: float,
    phi: float,
    pdfs: PDFPair,
    alpha_s: Callable[[float], float],
    run: RunConfig | None = None,
    scale_factor: float = 1.0,
    x_p_jacobian: float = 1.0,
) -> HelicityWeightBundle:
    """Attach Herwig RIVETWEIGHTS-style helicity ratios to local real carriers."""
    local = correlated_helicity_weights(
        q2,
        y,
        x_b,
        flavor,
        setup,
        contribution,
        base,
        x_p,
        z_p,
        phi,
        pdfs,
        alpha_s,
        run=run,
        scale_factor=scale_factor,
        x_p_jacobian=x_p_jacobian,
    )
    if contribution.lower() not in REAL_CONTRIBUTIONS:
        return local

    projected = correlated_helicity_weights(
        q2,
        y,
        x_b,
        flavor,
        setup,
        "inclusive",
        base,
        x_p,
        z_p,
        phi,
        pdfs,
        alpha_s,
        run=run,
        scale_factor=scale_factor,
        x_p_jacobian=x_p_jacobian,
    )
    if abs(projected.sigma) <= 1.0e-300 or not math.isfinite(projected.sigma):
        return local
    scale = local.nominal / projected.sigma
    sigma = projected.sigma * scale
    return HelicityWeightBundle(
        nominal=local.nominal,
        pp=projected.pp * scale,
        pm=projected.pm * scale,
        mp=projected.mp * scale,
        mm=projected.mm * scale,
        sigma=sigma,
        delta_ll=projected.delta_ll * scale,
    )


def polarized_real_weight_mode_from_args(args: argparse.Namespace) -> str:
    mode = getattr(args, "polarized_real_weight_mode", DEFAULT_POLARIZED_REAL_WEIGHT_MODE)
    if mode not in POLARIZED_REAL_WEIGHT_MODES:
        raise SystemExit(f"Unsupported --polarized-real-weight-mode '{mode}'.")
    return mode


def _correlated_bundle_for_mode(
    mode: str,
    q2: float,
    y: float,
    x_b: float,
    flavor: int,
    setup: str,
    contribution: str,
    base: float,
    x_p: float,
    z_p: float,
    phi: float,
    pdfs: PDFPair,
    alpha_s: Callable[[float], float],
    run: RunConfig | None = None,
    scale_factor: float = 1.0,
    x_p_jacobian: float = 1.0,
) -> HelicityWeightBundle:
    if mode == "projected-nlo":
        return projected_nlo_correlated_helicity_weights(
            q2,
            y,
            x_b,
            flavor,
            setup,
            contribution,
            base,
            x_p,
            z_p,
            phi,
            pdfs,
            alpha_s,
            run=run,
            scale_factor=scale_factor,
            x_p_jacobian=x_p_jacobian,
        )
    if mode == "local-real":
        return correlated_helicity_weights(
            q2,
            y,
            x_b,
            flavor,
            setup,
            contribution,
            base,
            x_p,
            z_p,
            phi,
            pdfs,
            alpha_s,
            run=run,
            scale_factor=scale_factor,
            x_p_jacobian=x_p_jacobian,
        )
    raise ValueError(f"Unsupported polarized real weight mode '{mode}'.")


def _event_with_correlated_weights(
    event: HepMCEvent,
    bundle: HelicityWeightBundle,
    weight_mode: str | None = None,
) -> HepMCEvent:
    weighted = event.with_weights(bundle.named_multipliers())
    if weight_mode is None:
        return weighted
    attributes = dict(weighted.attributes)
    attributes["HERWIG_FO_POL_WEIGHT_MODE"] = weight_mode
    return HepMCEvent(
        event_number=weighted.event_number,
        weight=weighted.weight,
        incoming=weighted.incoming,
        outgoing=weighted.outgoing,
        attributes=attributes,
        weights=weighted.weights,
    )


def generate_events(args: argparse.Namespace, progress: ProgressReporter | None = None) -> tuple[list[HepMCEvent], float]:
    if args.powheg:
        raise SystemExit("POWHEG-selected generation is intentionally unavailable in this Python fixed-order validator.")
    if args.events <= 0:
        raise SystemExit("--events must be positive.")

    run = run_config_from_args(args)
    profile = require_pdf_profile(args.pdf_profile)
    try:
        pdfs = load_pdf_pair(profile, use_toy=args.toy_pdfs)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    alpha_s = MatchboxAlphaS(run.alpha_s_mz, run.mz).alpha_s
    rng = random.Random(args.seed)
    q2_volume, _ = _q2_volume_and_jacobian(run.cuts.q2_min, run.cuts, args.q2_sampling)
    phase_space_volume = q2_volume * (run.cuts.y_max - run.cuts.y_min)
    events: list[HepMCEvent] = []
    event_number = 0

    if progress is not None:
        progress.start()
    for seed_index in range(args.events):
        q2_u, y_u, xp_u, z_u, phi_u = _seed_unit_coordinates(seed_index, rng, args.sampler, args.seed)
        q2, y, x_b, q2_jacobian = _coordinates_from_unit(q2_u, y_u, run, args.q2_sampling)
        x_p, x_p_jacobian = herwig_xp_power_map(xp_u, x_b)
        z_p = min(max(z_u, 1.0e-12), 1.0 - 1.0e-12)
        phi = 6.283185307179586 * phi_u
        for flavor in ACTIVE_QUARK_FLAVORS:
            point = _point_from_coordinates(
                q2,
                y,
                x_b,
                flavor,
                args.setup,
                args.helicity,
                run,
                scale_factor=float(getattr(args, "scale_factor", 1.0)),
            )
            born_pb = born_weight(point, pdfs)
            terms = nlo_terms(point, x_p=x_p, pdfs=pdfs, alpha_s=alpha_s, jacobian=x_p_jacobian)
            base = phase_space_volume * q2_jacobian * born_pb * x_p_jacobian / args.events

            inclusive_weight = base * terms.total
            events.append(born_event(point, inclusive_weight, event_number, run.beams))
            event_number += 1

            local_real = local_real_weight_ratios(point, x_p, z_p, phi, pdfs, alpha_s)
            qcdc_weight = base * local_real.qcdc
            events.append(real_event(point, "QCDC", x_p, z_p, phi, qcdc_weight, event_number, run.beams))
            event_number += 1
            events.append(born_event(point, -qcdc_weight, event_number, run.beams))
            event_number += 1

            bgf_weight = base * local_real.bgf
            events.append(real_event(point, "BGF", x_p, z_p, phi, bgf_weight, event_number, run.beams))
            event_number += 1
            events.append(born_event(point, -bgf_weight, event_number, run.beams))
            event_number += 1

        if progress is not None:
            progress.update(seed_index + 1)

    if progress is not None:
        progress.done()
    return events, sum(event.weight for event in events)


def generate_correlated_events(args: argparse.Namespace, progress: ProgressReporter | None = None) -> tuple[list[HepMCEvent], float]:
    if args.powheg:
        raise SystemExit("POWHEG-selected generation is intentionally unavailable in this Python fixed-order validator.")
    if args.events <= 0:
        raise SystemExit("--events must be positive.")

    run = run_config_from_args(args)
    profile = require_pdf_profile(args.pdf_profile)
    try:
        pdfs = load_pdf_pair(profile, use_toy=args.toy_pdfs)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    alpha_s = MatchboxAlphaS(run.alpha_s_mz, run.mz).alpha_s
    rng = random.Random(args.seed)
    q2_volume, _ = _q2_volume_and_jacobian(run.cuts.q2_min, run.cuts, args.q2_sampling)
    phase_space_volume = q2_volume * (run.cuts.y_max - run.cuts.y_min)
    scale_factor = float(getattr(args, "scale_factor", 1.0))
    weight_mode = polarized_real_weight_mode_from_args(args)
    events: list[HepMCEvent] = []
    event_number = 0

    if progress is not None:
        progress.start()
    for seed_index in range(args.events):
        q2_u, y_u, xp_u, z_u, phi_u = _seed_unit_coordinates(seed_index, rng, args.sampler, args.seed)
        q2, y, x_b, q2_jacobian = _coordinates_from_unit(q2_u, y_u, run, args.q2_sampling)
        x_p, x_p_jacobian = herwig_xp_power_map(xp_u, x_b)
        z_p = min(max(z_u, 1.0e-12), 1.0 - 1.0e-12)
        phi = 6.283185307179586 * phi_u
        base = phase_space_volume * q2_jacobian * x_p_jacobian / args.events
        for flavor in ACTIVE_QUARK_FLAVORS:
            carrier_point = _point_from_coordinates(
                q2,
                y,
                x_b,
                flavor,
                args.setup,
                "00",
                run,
                scale_factor=scale_factor,
            )

            inclusive = _correlated_bundle_for_mode(
                weight_mode,
                q2,
                y,
                x_b,
                flavor,
                args.setup,
                "inclusive",
                base,
                x_p,
                z_p,
                phi,
                pdfs,
                alpha_s,
                run=run,
                scale_factor=scale_factor,
                x_p_jacobian=x_p_jacobian,
            )
            events.append(
                _event_with_correlated_weights(
                    born_event(carrier_point, inclusive.nominal, event_number, run.beams),
                    inclusive,
                    weight_mode,
                )
            )
            event_number += 1

            qcdc = _correlated_bundle_for_mode(
                weight_mode,
                q2,
                y,
                x_b,
                flavor,
                args.setup,
                "qcdc",
                base,
                x_p,
                z_p,
                phi,
                pdfs,
                alpha_s,
                run=run,
                scale_factor=scale_factor,
                x_p_jacobian=x_p_jacobian,
            )
            events.append(
                _event_with_correlated_weights(
                    real_event(carrier_point, "QCDC", x_p, z_p, phi, qcdc.nominal, event_number, run.beams),
                    qcdc,
                    weight_mode,
                )
            )
            event_number += 1

            qcdc_counter = _correlated_bundle_for_mode(
                weight_mode,
                q2,
                y,
                x_b,
                flavor,
                args.setup,
                "qcdc_counter",
                base,
                x_p,
                z_p,
                phi,
                pdfs,
                alpha_s,
                run=run,
                scale_factor=scale_factor,
                x_p_jacobian=x_p_jacobian,
            )
            events.append(
                _event_with_correlated_weights(
                    born_event(carrier_point, qcdc_counter.nominal, event_number, run.beams),
                    qcdc_counter,
                    weight_mode,
                )
            )
            event_number += 1

            bgf = _correlated_bundle_for_mode(
                weight_mode,
                q2,
                y,
                x_b,
                flavor,
                args.setup,
                "bgf",
                base,
                x_p,
                z_p,
                phi,
                pdfs,
                alpha_s,
                run=run,
                scale_factor=scale_factor,
                x_p_jacobian=x_p_jacobian,
            )
            events.append(
                _event_with_correlated_weights(
                    real_event(carrier_point, "BGF", x_p, z_p, phi, bgf.nominal, event_number, run.beams),
                    bgf,
                    weight_mode,
                )
            )
            event_number += 1

            bgf_counter = _correlated_bundle_for_mode(
                weight_mode,
                q2,
                y,
                x_b,
                flavor,
                args.setup,
                "bgf_counter",
                base,
                x_p,
                z_p,
                phi,
                pdfs,
                alpha_s,
                run=run,
                scale_factor=scale_factor,
                x_p_jacobian=x_p_jacobian,
            )
            events.append(
                _event_with_correlated_weights(
                    born_event(carrier_point, bgf_counter.nominal, event_number, run.beams),
                    bgf_counter,
                    weight_mode,
                )
            )
            event_number += 1

        if progress is not None:
            progress.update(seed_index + 1)

    if progress is not None:
        progress.done()
    return events, sum(event.weight for event in events)


def _midpoints(count: int) -> tuple[float, ...]:
    if count <= 0:
        raise SystemExit("Integration bin counts must be positive.")
    return tuple((index + 0.5) / count for index in range(count))


def _observable_components(observable: str, helicity: str) -> tuple[tuple[str, float], ...]:
    if observable == "helicity":
        return ((helicity.upper(), 1.0),)
    if observable == "sigma0":
        return (("00", 1.0),)
    if observable == "delta":
        return (("PP", 0.25), ("MM", 0.25), ("PM", -0.25), ("MP", -0.25))
    raise ValueError(f"Unsupported observable '{observable}'.")


def integrate_helicity(
    run: RunConfig,
    setup: str,
    helicity: str,
    order: str,
    pdfs: PDFPair,
    alpha_s,
    q2_sampling: str,
    q2_bins: int,
    y_bins: int,
    xp_bins: int,
    scale_factor: float = 1.0,
) -> float:
    order = order.upper()
    if order not in ("LO", "NLO"):
        raise ValueError(f"Unsupported order '{order}'.")
    q2_points = _midpoints(q2_bins)
    y_points = _midpoints(y_bins)
    xp_points = _midpoints(xp_bins) if order == "NLO" else ()
    q2_volume, _ = _q2_volume_and_jacobian(run.cuts.q2_min, run.cuts, q2_sampling)
    q2_step = q2_volume / q2_bins
    y_step = (run.cuts.y_max - run.cuts.y_min) / y_bins
    total = 0.0

    for q2_u in q2_points:
        if q2_sampling == "linear":
            q2 = run.cuts.q2_min + q2_u * q2_volume
            q2_jacobian = 1.0
        elif q2_sampling == "log":
            q2 = math.exp(math.log(run.cuts.q2_min) + q2_u * q2_volume)
            q2_jacobian = q2
        else:
            raise ValueError(f"Unsupported Q2 sampling '{q2_sampling}'.")
        q2_weight = q2_step * q2_jacobian
        for y_u in y_points:
            y = run.cuts.y_min + y_u * (run.cuts.y_max - run.cuts.y_min)
            x_b = xb_from_q2_y(q2, y, run.beams)
            if not (0.0 < x_b < 1.0):
                continue
            cell_weight = q2_weight * y_step
            for flavor in ACTIVE_QUARK_FLAVORS:
                point = _point_from_coordinates(q2, y, x_b, flavor, setup, helicity, run, scale_factor=scale_factor)
                born_pb = born_weight(point, pdfs)
                if order == "LO":
                    total += cell_weight * born_pb
                    continue
                for xp_u in xp_points:
                    x_p, x_p_jacobian = herwig_xp_power_map(xp_u, x_b)
                    terms = nlo_terms(point, x_p=x_p, pdfs=pdfs, alpha_s=alpha_s, jacobian=x_p_jacobian)
                    total += cell_weight * born_pb * x_p_jacobian * terms.total / xp_bins

    return total


def integrate_observable(
    run: RunConfig,
    setup: str,
    observable: str,
    helicity: str,
    order: str,
    pdfs: PDFPair,
    alpha_s,
    q2_sampling: str,
    q2_bins: int,
    y_bins: int,
    xp_bins: int,
    scale_factor: float = 1.0,
) -> float:
    total = 0.0
    for label, coefficient in _observable_components(observable, helicity):
        total += coefficient * integrate_helicity(
            run,
            setup,
            label,
            order,
            pdfs,
            alpha_s,
            q2_sampling,
            q2_bins,
            y_bins,
            xp_bins,
            scale_factor=scale_factor,
        )
    return total


def command_generate(args: argparse.Namespace) -> int:
    progress = None
    if _progress_stream(args) is not None:
        progress = ProgressReporter(f"generate {args.setup}/{args.helicity}", args.events, unit="seeds")
    events, sumw = generate_events(args, progress=progress)
    with _timed_step(args, f"write HepMC {args.output}"):
        xsec = write_hepmc3(args.output, events, cross_section_pb=sumw)
    print(f"Wrote {len(events)} signed events to {args.output}")
    print(f"Set HepMC GenCrossSection to signed sumW = {xsec:.12e} pb")
    return 0


def command_integrate(args: argparse.Namespace) -> int:
    if args.powheg:
        raise SystemExit("POWHEG-selected generation is intentionally unavailable in this Python fixed-order validator.")
    run = run_config_from_args(args)
    profile = require_pdf_profile(args.pdf_profile)
    try:
        pdfs = load_pdf_pair(profile, use_toy=args.toy_pdfs)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    alpha_s = MatchboxAlphaS(run.alpha_s_mz, run.mz).alpha_s
    setup = require_setup(args.setup)
    orders = ("LO", "NLO") if args.order == "both" else (args.order,)
    print(
        "cuts: "
        f"Q2=[{run.cuts.q2_min:g}, {run.cuts.q2_max:g}] GeV^2, "
        f"y=[{run.cuts.y_min:g}, {run.cuts.y_max:g}], "
        f"Q2 sampling={args.q2_sampling}"
    )
    for order in orders:
        value = integrate_observable(
            run,
            setup,
            args.observable,
            args.helicity,
            order,
            pdfs,
            alpha_s,
            args.q2_sampling,
            args.q2_bins,
            args.y_bins,
            args.xp_bins,
            scale_factor=float(getattr(args, "scale_factor", 1.0)),
        )
        label = args.helicity.upper() if args.observable == "helicity" else args.observable
        print(f"{setup} {order} {label} = {value:.12f} pb")
    return 0


def command_rivet(args: argparse.Namespace) -> int:
    plugin = args.analysis_plugin
    if args.build_analysis or plugin is None:
        with _timed_step(args, "build MC_DIS_BREIT plugin"):
            plugin = build_mc_dis_breit(output_dir=(args.output.parent / "rivet-plugin"))
    with _timed_step(args, f"run Rivet {args.hepmc}"):
        run_rivet(args.hepmc, args.output, jetinput=args.jetinput, analysis_plugin=plugin, rivetweights=args.rivetweights)
    print(f"Wrote Rivet YODA output to {args.output}")
    return 0


def command_compare(args: argparse.Namespace) -> int:
    summary = compare_yoda(args.herwig, args.reference, max_pull=args.max_pull)
    print(f"max |pull| = {summary.max_abs_pull:.3g}")
    if summary.missing_in_herwig:
        print(f"missing in Herwig FO: {len(summary.missing_in_herwig)} objects")
    if not summary.ok:
        return 2
    print("YODA comparison passed within combined statistical uncertainty.")
    return 0


def _namespace_with(args: argparse.Namespace, **updates) -> argparse.Namespace:
    data = vars(args).copy()
    data.update(updates)
    return argparse.Namespace(**data)


def command_polarize(args: argparse.Namespace) -> int:
    try:
        from analyze_DIS_polarized import build_dis_polarized_objects, write_yoda_gz
    except ImportError as exc:
        raise SystemExit("Could not import YODA/Rivet helpers. Load a Rivet/YODA Python environment first.") from exc

    with _timed_step(args, "combine helicity YODAs"):
        objects = build_dis_polarized_objects(
            setup=args.setup,
            zero_path=str(args.zero) if args.zero is not None else None,
            pp_path=str(args.PP),
            pm_path=str(args.PM),
            mp_path=str(args.MP) if args.MP is not None else None,
            mm_path=str(args.MM) if args.MM is not None else None,
            analysis=args.analysis,
        )
    with _timed_step(args, f"write polarized YODA {args.output}"):
        write_yoda_gz(objects, str(args.output))
    delta_count = sum(1 for path in objects if f"/{args.analysis}/D" in path)
    asymmetry_count = sum(1 for path in objects if f"/{args.analysis}/ALL" in path)
    print(f"Wrote polarized YODA output to {args.output}")
    print(f"Included {delta_count} D* histograms and {asymmetry_count} ALL* histograms.")
    return 0


def _command_polarized_full_independent(args: argparse.Namespace) -> int:
    if args.powheg:
        raise SystemExit("POWHEG-selected generation is intentionally unavailable in this Python fixed-order validator.")
    if args.events <= 0:
        raise SystemExit("--events must be positive.")

    args.work_dir.mkdir(parents=True, exist_ok=True)
    plugin = args.analysis_plugin
    helicities = required_polarized_helicities(args.setup)
    total_steps = len(helicities) * 3 + 1 + (1 if (args.build_analysis or plugin is None) else 0)
    step = 0
    if args.build_analysis or plugin is None:
        step += 1
        _numbered_step(args, step, total_steps, "build MC_DIS_BREIT plugin")
        with _timed_step(args, "build MC_DIS_BREIT plugin"):
            plugin = build_mc_dis_breit(output_dir=(args.work_dir / "rivet-plugin"))

    yoda_paths: dict[str, Path] = {}
    for helicity in helicities:
        hepmc_path = args.work_dir / f"standalone_fo_{args.setup.lower()}_{helicity.lower()}.hepmc3"
        yoda_path = args.work_dir / f"standalone_fo_{args.setup.lower()}_{helicity.lower()}.yoda"
        gen_args = _namespace_with(args, helicity=helicity, output=hepmc_path)
        step += 1
        _numbered_step(args, step, total_steps, f"generate {args.setup}/{helicity}")
        progress = None
        if _progress_stream(args) is not None:
            progress = ProgressReporter(f"generate {args.setup}/{helicity}", args.events, unit="seeds")
        events, sumw = generate_events(gen_args, progress=progress)
        step += 1
        _numbered_step(args, step, total_steps, f"write HepMC {helicity}")
        with _timed_step(args, f"write HepMC {helicity}"):
            write_hepmc3(hepmc_path, events, cross_section_pb=sumw)
        step += 1
        _numbered_step(args, step, total_steps, f"run Rivet {helicity}")
        with _timed_step(args, f"run Rivet {helicity}"):
            run_rivet(hepmc_path, yoda_path, jetinput=args.jetinput, analysis_plugin=plugin)
        yoda_paths[helicity] = yoda_path
        print(f"{helicity}: wrote {len(events)} events, sumW={sumw:.12e} pb, YODA={yoda_path}")

    step += 1
    _numbered_step(args, step, total_steps, "combine polarized YODA")
    return command_polarize(
        argparse.Namespace(
            setup=args.setup,
            analysis=args.analysis,
            zero=yoda_paths.get("00"),
            PP=yoda_paths.get("PP"),
            PM=yoda_paths.get("PM"),
            MP=yoda_paths.get("MP"),
            MM=yoda_paths.get("MM"),
            output=args.output,
            no_progress=getattr(args, "no_progress", False),
        )
    )


def _default_variation_output_path(output: Path) -> Path:
    name = output.name
    if name.endswith(".yoda.gz"):
        return output.with_name(name[:-8] + ".scale_variations.yoda.gz")
    if name.endswith(".yoda"):
        return output.with_name(name[:-5] + ".scale_variations.yoda")
    return output.with_name(name + ".scale_variations.yoda")


def _derive_correlated_yoda(input_path: Path, output_path: Path, analysis: str) -> tuple[int, int]:
    try:
        import yoda
        from analyze_DIS_polarized import (
            build_correlated_helicity_derived_objects,
            normalize_analysis_path_objects,
            write_yoda_gz,
        )
    except ImportError as exc:
        raise SystemExit("Could not import YODA/Rivet helpers. Load a Rivet/YODA Python environment first.") from exc

    objects = normalize_analysis_path_objects(yoda.read(str(input_path)), analysis=analysis, keep_raw=False)
    derived = build_correlated_helicity_derived_objects(objects, analysis=analysis)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_yoda_gz(derived, str(output_path))
    delta_count = sum(1 for path in derived if f"/{analysis}/D" in path)
    asymmetry_count = sum(1 for path in derived if f"/{analysis}/ALL" in path)
    return delta_count, asymmetry_count


def _write_variation_tagged_yoda(
    inputs: dict[str, Path],
    output_path: Path,
    analysis: str,
) -> int:
    try:
        import yoda
        from analyze_DIS_polarized import clone_objects_with_variation, write_yoda_gz
    except ImportError as exc:
        raise SystemExit("Could not import YODA/Rivet helpers. Load a Rivet/YODA Python environment first.") from exc

    tagged: dict[str, object] = {}
    for variation_name, path in inputs.items():
        objects = yoda.read(str(path))
        if variation_name == "nominal":
            tagged.update(objects)
            continue
        tagged.update(clone_objects_with_variation(objects, variation_name))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_yoda_gz(tagged, str(output_path))
    return len(tagged)


def _run_correlated_polarized_variation(
    args: argparse.Namespace,
    plugin: Path,
    variation_name: str,
    scale_factor: float,
) -> tuple[Path, float, int]:
    suffix = "" if variation_name == "nominal" else f"_{variation_name.lower()}"
    hepmc_path = args.work_dir / f"standalone_fo_{args.setup.lower()}_correlated{suffix}.hepmc3"
    raw_yoda_path = args.work_dir / f"standalone_fo_{args.setup.lower()}_correlated{suffix}.raw.yoda"
    clean_yoda_path = args.output if variation_name == "nominal" else (
        args.work_dir / f"standalone_fo_{args.setup.lower()}_correlated{suffix}.yoda"
    )
    _, sumw, event_count = _run_correlated_polarized_raw_variation(
        args,
        plugin,
        variation_name,
        scale_factor,
        hepmc_path,
        raw_yoda_path,
    )
    with _timed_step(args, f"derive correlated YODA {variation_name}"):
        _derive_correlated_yoda(raw_yoda_path, clean_yoda_path, args.analysis)
    print(
        f"{variation_name}: wrote {event_count} correlated events, "
        f"sumW={sumw:.12e} pb, scale-factor={scale_factor:g}, "
        f"polarized-real-mode={polarized_real_weight_mode_from_args(args)}, YODA={clean_yoda_path}"
    )
    return clean_yoda_path, sumw, event_count


def _run_correlated_polarized_raw_variation(
    args: argparse.Namespace,
    plugin: Path,
    variation_name: str,
    scale_factor: float,
    hepmc_path: Path,
    raw_yoda_path: Path,
) -> tuple[Path, float, int]:
    gen_args = _namespace_with(args, output=hepmc_path, helicity="00", scale_factor=scale_factor)

    progress = None
    if _progress_stream(args) is not None:
        progress = ProgressReporter(
            f"generate {args.setup}/correlated/{variation_name}",
            args.events,
            unit="seeds",
        )
    events, sumw = generate_correlated_events(gen_args, progress=progress)
    with _timed_step(args, f"write correlated HepMC {variation_name}"):
        write_hepmc3(hepmc_path, events, cross_section_pb=sumw)
    with _timed_step(args, f"run correlated Rivet {variation_name}"):
        run_rivet(hepmc_path, raw_yoda_path, jetinput=args.jetinput, analysis_plugin=plugin, rivetweights=True)
    return raw_yoda_path, sumw, len(events)


def command_polarized_shard(args: argparse.Namespace) -> int:
    if args.powheg:
        raise SystemExit("POWHEG-selected generation is intentionally unavailable in this Python fixed-order validator.")
    if args.events <= 0:
        raise SystemExit("--events must be positive.")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    plugin = args.analysis_plugin
    if args.build_analysis or plugin is None:
        with _timed_step(args, "build MC_DIS_BREIT plugin"):
            plugin = build_mc_dis_breit(output_dir=(args.work_dir / "rivet-plugin"))
    _, sumw, event_count = _run_correlated_polarized_raw_variation(
        args,
        plugin,
        args.variation_name,
        float(args.scale_factor),
        args.hepmc,
        args.raw_yoda,
    )
    if args.delete_hepmc_after_rivet:
        try:
            args.hepmc.unlink()
        except FileNotFoundError:
            pass
    print(
        f"{args.variation_name}: wrote {event_count} correlated events, "
        f"sumW={sumw:.12e} pb, raw YODA={args.raw_yoda}"
    )
    return 0


def command_polarized_full(args: argparse.Namespace) -> int:
    if args.independent_helicities:
        if args.scale_variations:
            raise SystemExit("--scale-variations is available for the correlated default path, not --independent-helicities.")
        return _command_polarized_full_independent(args)
    if args.powheg:
        raise SystemExit("POWHEG-selected generation is intentionally unavailable in this Python fixed-order validator.")
    if args.events <= 0:
        raise SystemExit("--events must be positive.")

    args.work_dir.mkdir(parents=True, exist_ok=True)
    plugin = args.analysis_plugin
    variations = SCALE_VARIATION_FACTORS if args.scale_variations else {"nominal": float(args.scale_factor)}
    total_steps = len(variations) + (1 if (args.build_analysis or plugin is None) else 0)
    if args.scale_variations:
        total_steps += 1
    step = 0
    if args.build_analysis or plugin is None:
        step += 1
        _numbered_step(args, step, total_steps, "build MC_DIS_BREIT plugin")
        with _timed_step(args, "build MC_DIS_BREIT plugin"):
            plugin = build_mc_dis_breit(output_dir=(args.work_dir / "rivet-plugin"))

    yoda_paths: dict[str, Path] = {}
    for variation_name, scale_factor in variations.items():
        step += 1
        _numbered_step(args, step, total_steps, f"generate correlated {variation_name}")
        clean_yoda_path, _, _ = _run_correlated_polarized_variation(args, plugin, variation_name, scale_factor)
        yoda_paths[variation_name] = clean_yoda_path

    if args.scale_variations:
        variation_output = args.variation_output or _default_variation_output_path(args.output)
        step += 1
        _numbered_step(args, step, total_steps, f"write variation-tagged YODA {variation_output}")
        with _timed_step(args, f"write variation-tagged YODA {variation_output}"):
            object_count = _write_variation_tagged_yoda(yoda_paths, variation_output, args.analysis)
        print(f"Wrote variation-tagged YODA output to {variation_output} with {object_count} objects.")

    print(f"Wrote correlated polarized YODA output to {args.output}")
    return 0


def command_full(args: argparse.Namespace) -> int:
    progress = None
    if _progress_stream(args) is not None:
        progress = ProgressReporter(f"generate {args.setup}/{args.helicity}", args.events, unit="seeds")
    events, sumw = generate_events(args, progress=progress)
    with _timed_step(args, f"write HepMC {args.hepmc}"):
        write_hepmc3(args.hepmc, events, cross_section_pb=sumw)
    plugin = args.analysis_plugin
    if args.build_analysis or plugin is None:
        with _timed_step(args, "build MC_DIS_BREIT plugin"):
            plugin = build_mc_dis_breit(output_dir=(args.yoda.parent / "rivet-plugin"))
    with _timed_step(args, f"run Rivet {args.hepmc}"):
        run_rivet(args.hepmc, args.yoda, jetinput=args.jetinput, analysis_plugin=plugin)
    with _timed_step(args, f"compare {args.yoda}"):
        summary = compare_yoda(args.yoda, args.reference, max_pull=args.max_pull)
    print(f"Wrote {len(events)} events to {args.hepmc}")
    print(f"Wrote Rivet YODA output to {args.yoda}")
    print(f"max |pull| = {summary.max_abs_pull:.3g}")
    return 0 if summary.ok else 2


def command_plot_prep(args: argparse.Namespace) -> int:
    try:
        from analyze_DIS_polarized import (
            build_scale_envelope_plot_yoda,
            harmonize_plot_yoda_bin_grids,
            write_plot_safe_yoda,
        )
    except ImportError as exc:
        raise SystemExit("Could not import YODA/Rivet helpers. Load a Rivet/YODA Python environment first.") from exc

    if args.herwig_scale_envelope:
        herwig_count, herwig_dropped, _ = build_scale_envelope_plot_yoda(args.herwig, args.herwig_output)
    else:
        herwig_count, herwig_dropped = write_plot_safe_yoda(args.herwig, args.herwig_output, analysis=args.analysis)
    poldis_count, poldis_dropped = write_plot_safe_yoda(args.poldis, args.poldis_output, analysis=args.analysis)
    updated_files, updated_objects = harmonize_plot_yoda_bin_grids([args.herwig_output, args.poldis_output])
    print(
        f"Wrote plot-safe Herwig YODA to {args.herwig_output} "
        f"({herwig_count} objects, {herwig_dropped} skipped)."
    )
    print(
        f"Wrote plot-safe POLDIS YODA to {args.poldis_output} "
        f"({poldis_count} objects, {poldis_dropped} skipped)."
    )
    if updated_objects:
        print(f"Harmonized {updated_objects} object grids across {updated_files} files.")
    return 0


def _repo_root_from_script() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "workflow" / "dispol").exists() and (parent / "README.md").exists():
            return parent
        if (parent / "DISPOL").exists() and (parent / "README.md").exists():
            return parent
    return Path(__file__).resolve().parents[2]


def _resolve_path(path: Path, *, repo_root: Path, prefer_repo: bool = False) -> Path:
    if path.is_absolute():
        return path
    cwd_path = (Path.cwd() / path).resolve()
    repo_path = (repo_root / path).resolve()
    if prefer_repo and repo_path.exists():
        return repo_path
    if cwd_path.exists():
        return cwd_path
    if repo_path.exists():
        return repo_path
    return cwd_path


def _format_suffix_value(value: float) -> str:
    return f"{value:g}".replace(".", "p").replace("-", "m")


def _standalone_campaign_suffix(args: argparse.Namespace) -> str:
    run = run_config_from_args(args)
    q2_min = run.cuts.q2_min
    if abs(q2_min - 100.0) <= 1.0e-12:
        return "q2gt100"
    return f"q2gt{_format_suffix_value(q2_min)}"


def _prepend_env_path(env: dict[str, str], name: str, paths: Sequence[Path]) -> None:
    entries = [str(path) for path in paths if path is not None]
    existing = env.get(name)
    if existing:
        entries.append(existing)
    env[name] = os.pathsep.join(entries)


def _standalone_campaign_paths(args: argparse.Namespace) -> dict[str, Path]:
    repo_root = _repo_root_from_script()
    base_dir = _resolve_path(args.base_dir, repo_root=repo_root, prefer_repo=True)
    campaign_dir = args.campaign_dir.resolve() if args.campaign_dir is not None else (base_dir / "campaigns" / args.tag)
    analysis_dir = campaign_dir / "analysis"
    work_dir = campaign_dir / "work"
    suffix = _standalone_campaign_suffix(args)
    output_stem = args.output_stem or f"standalone_fo_{args.setup.lower()}_correlated_{suffix}"
    poldis_stem = args.poldis_output_stem or f"poldis_reference_{args.setup.lower()}_{suffix}"
    nominal_yoda = analysis_dir / f"{output_stem}.yoda"
    variation_yoda = analysis_dir / f"{output_stem}.scale_variations.yoda"
    if args.shard_status_file is None:
        shard_status_file = work_dir / "shard_status.tsv"
    elif args.shard_status_file.is_absolute():
        shard_status_file = args.shard_status_file
    else:
        shard_status_file = campaign_dir / args.shard_status_file
    herwig_plot_yoda = analysis_dir / (
        f"{output_stem}.scale_envelope_plot.yoda" if args.scale_variations else f"{output_stem}.plot_safe.yoda"
    )
    poldis_plot_yoda = analysis_dir / f"{poldis_stem}.plot_safe.yoda"
    plots_dir = args.plots_dir.resolve() if args.plots_dir is not None else (
        campaign_dir / "plots_mc_vs_poldis_all_nlo_with_polarized"
    )
    analysis_rivet_dir = (
        args.analysis_rivet_dir.resolve()
        if args.analysis_rivet_dir is not None
        else base_dir / "analyses" / "rivet" / "dis"
    )
    plot_config = args.plot_config.resolve() if args.plot_config is not None else analysis_rivet_dir / f"{args.analysis}.plot"
    poldis_reference = _resolve_path(args.poldis_reference, repo_root=repo_root, prefer_repo=True)
    rivet_mkhtml_safe = (
        args.rivet_mkhtml_safe.resolve()
        if args.rivet_mkhtml_safe is not None
        else Path(__file__).resolve().with_name("rivet_mkhtml_safe.py")
    )
    analysis_plugin = args.analysis_plugin
    if analysis_plugin is not None:
        analysis_plugin = _resolve_path(analysis_plugin, repo_root=repo_root)
    return {
        "repo_root": repo_root,
        "base_dir": base_dir,
        "campaign_dir": campaign_dir,
        "analysis_dir": analysis_dir,
        "work_dir": work_dir,
        "nominal_yoda": nominal_yoda,
        "variation_yoda": variation_yoda,
        "shard_status_file": shard_status_file,
        "herwig_plot_yoda": herwig_plot_yoda,
        "poldis_plot_yoda": poldis_plot_yoda,
        "plots_dir": plots_dir,
        "analysis_rivet_dir": analysis_rivet_dir,
        "plot_config": plot_config,
        "poldis_reference": poldis_reference,
        "rivet_mkhtml_safe": rivet_mkhtml_safe,
        "analysis_plugin": analysis_plugin,
    }


def _choose_yoda_merge_tool(preference: str) -> str | None:
    if preference == "auto":
        for candidate in ("yodamerge", "rivet-merge"):
            resolved = shutil.which(candidate)
            if resolved is not None:
                return resolved
        return None
    explicit = Path(preference)
    if explicit.exists():
        return str(explicit)
    return shutil.which(preference)


def _merge_yoda_inputs(
    inputs: Sequence[Path],
    output_path: Path,
    *,
    base_dir: Path,
    merge_tool_preference: str,
    label: str,
) -> None:
    unique_inputs = sorted({str(path.resolve()) for path in inputs})
    if not unique_inputs:
        raise SystemExit(f"No YODA inputs found for {label}.")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if len(unique_inputs) == 1:
        shutil.copy2(unique_inputs[0], output_path)
        return

    merge_tool = _choose_yoda_merge_tool(merge_tool_preference)
    if merge_tool is None:
        raise SystemExit("No YODA merge tool found. Load Rivet/YODA or pass --yoda-merge-tool.")

    tool_name = Path(merge_tool).name
    command = [merge_tool]
    if tool_name == "rivet-merge":
        command.append("-e")
    command.extend([*unique_inputs, "-o", str(output_path)])
    proc = subprocess.run(command, cwd=base_dir, text=True, capture_output=True)
    log_stem = output_path.with_suffix("")
    if proc.stdout:
        log_stem.with_name(log_stem.name + ".merge.stdout").write_text(proc.stdout)
    if proc.stderr:
        log_stem.with_name(log_stem.name + ".merge.stderr").write_text(proc.stderr)
    if proc.returncode != 0:
        raise SystemExit(f"YODA merge failed for {label} with rc={proc.returncode}: {' '.join(command)}")


def _append_optional_float_arg(command: list[str], option: str, value: float | None) -> None:
    if value is not None:
        command.extend([option, f"{float(value):.17g}"])


def _standalone_shard_task(
    args: argparse.Namespace,
    paths: dict[str, Path],
    *,
    plugin: Path,
    variation_name: str,
    scale_factor: float,
    spec: StandaloneShardSpec,
) -> StandaloneShardTask:
    variation_slug = variation_name.lower()
    shard_dir = paths["work_dir"] / "shards" / variation_slug / f"shard_{spec.index:04d}"
    hepmc_path = shard_dir / f"standalone_fo_{args.setup.lower()}_{variation_slug}_shard_{spec.index:04d}.hepmc3"
    raw_yoda_path = shard_dir / f"standalone_fo_{args.setup.lower()}_{variation_slug}_shard_{spec.index:04d}.raw.yoda"
    stdout_path = shard_dir / "stdout.log"
    stderr_path = shard_dir / "stderr.log"

    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "polarized-shard",
        "--setup",
        args.setup,
        "--pdf-profile",
        args.pdf_profile,
        "--window",
        args.window,
        "--q2-sampling",
        args.q2_sampling,
        "--sampler",
        args.sampler,
        "--seed",
        str(spec.seed),
        "--scale-factor",
        f"{scale_factor:.17g}",
        "--jetinput",
        args.jetinput,
        "--events",
        str(spec.events),
        "--work-dir",
        str(shard_dir),
        "--hepmc",
        str(hepmc_path),
        "--raw-yoda",
        str(raw_yoda_path),
        "--analysis",
        args.analysis,
        "--analysis-plugin",
        str(plugin),
        "--variation-name",
        variation_name,
        "--polarized-real-weight-mode",
        args.polarized_real_weight_mode,
        "--no-progress",
    ]
    _append_optional_float_arg(command, "--q2-min", args.q2_min)
    _append_optional_float_arg(command, "--q2-max", args.q2_max)
    _append_optional_float_arg(command, "--y-min", args.y_min)
    _append_optional_float_arg(command, "--y-max", args.y_max)
    if args.toy_pdfs:
        command.append("--toy-pdfs")
    if not args.keep_shard_hepmc:
        command.append("--delete-hepmc-after-rivet")

    return StandaloneShardTask(
        variation_name=variation_name,
        scale_factor=scale_factor,
        spec=spec,
        work_dir=shard_dir,
        hepmc_path=hepmc_path,
        raw_yoda_path=raw_yoda_path,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        command=tuple(command),
    )


def _run_standalone_shard_task(task: StandaloneShardTask) -> StandaloneShardResult:
    task.work_dir.mkdir(parents=True, exist_ok=True)
    with task.stdout_path.open("w") as stdout, task.stderr_path.open("w") as stderr:
        proc = subprocess.run(list(task.command), text=True, stdout=stdout, stderr=stderr)
    return StandaloneShardResult(task=task, returncode=proc.returncode)


def _run_parallel_standalone_generation(args: argparse.Namespace, paths: dict[str, Path]) -> None:
    if args.jobs <= 0:
        raise SystemExit("--jobs must be positive.")
    if args.shard_monitor_interval <= 0.0:
        raise SystemExit("--shard-monitor-interval must be positive.")
    specs = _standalone_shard_specs(args.events, args.shards, args.seed)
    paths["work_dir"].mkdir(parents=True, exist_ok=True)

    plugin = paths["analysis_plugin"]
    if args.build_analysis or plugin is None:
        with _timed_step(args, "build MC_DIS_BREIT plugin"):
            plugin = build_mc_dis_breit(output_dir=(paths["work_dir"] / "rivet-plugin"))

    variations = SCALE_VARIATION_FACTORS if args.scale_variations else {"nominal": float(args.scale_factor)}
    tasks = [
        _standalone_shard_task(
            args,
            paths,
            plugin=plugin,
            variation_name=variation_name,
            scale_factor=scale_factor,
            spec=spec,
        )
        for variation_name, scale_factor in variations.items()
        for spec in specs
    ]

    worker_count = max(1, min(int(args.jobs), len(tasks)))
    _progress_message(
        args,
        f"[standalone-campaign] launching {len(tasks)} raw Rivet shards "
        f"({len(specs)} shards x {len(variations)} scale point(s), max_workers={worker_count})",
    )
    progress = None
    if _progress_stream(args) is not None:
        progress = ProgressReporter("standalone raw shards", len(tasks), unit="shards")
        progress.start()

    failures: list[StandaloneShardResult] = []
    states = {_standalone_shard_key(task): "pending" for task in tasks}
    returncodes: dict[tuple[str, int], int] = {}
    start_times: dict[tuple[str, int], float] = {}
    end_times: dict[tuple[str, int], float] = {}
    status_path = paths["shard_status_file"]
    started_at = time.monotonic()

    def write_monitor(force_message: bool = False) -> None:
        now = time.monotonic()
        _write_shard_status_file(status_path, tasks, states, returncodes, start_times, end_times, now=now)
        if force_message or _progress_stream(args) is not None:
            finished = sum(1 for state in states.values() if state in {"completed", "failed"})
            ok = sum(1 for state in states.values() if state == "completed")
            failed = sum(1 for state in states.values() if state == "failed")
            running = sum(1 for state in states.values() if state == "running")
            pending = sum(1 for state in states.values() if state == "pending")
            _progress_message(
                args,
                _shard_monitor_message(
                    total=len(tasks),
                    finished=finished,
                    ok=ok,
                    failed=failed,
                    running=running,
                    pending=pending,
                    started_at=started_at,
                    now=now,
                    status_path=status_path,
                ),
            )

    completed = 0
    pending_tasks = list(tasks)
    active: dict[concurrent.futures.Future[StandaloneShardResult], StandaloneShardTask] = {}

    def launch_until_full(executor: concurrent.futures.Executor) -> None:
        while pending_tasks and len(active) < worker_count:
            task = pending_tasks.pop(0)
            key = _standalone_shard_key(task)
            states[key] = "running"
            start_times[key] = time.monotonic()
            active[executor.submit(_run_standalone_shard_task, task)] = task

    write_monitor(force_message=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=worker_count) as executor:
        launch_until_full(executor)
        write_monitor(force_message=True)
        while active:
            done, _ = concurrent.futures.wait(
                set(active),
                timeout=float(args.shard_monitor_interval),
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            if not done:
                write_monitor(force_message=True)
                continue
            for future in done:
                task = active.pop(future)
                key = _standalone_shard_key(task)
                try:
                    result = future.result()
                    returncode = result.returncode
                except Exception as exc:
                    task.stderr_path.parent.mkdir(parents=True, exist_ok=True)
                    with task.stderr_path.open("a") as err:
                        err.write(f"\ninternal shard runner failure: {exc}\n")
                    result = StandaloneShardResult(task=task, returncode=1)
                    returncode = 1
                end_times[key] = time.monotonic()
                returncodes[key] = returncode
                states[key] = "completed" if returncode == 0 else "failed"
                completed += 1
                if progress is not None:
                    progress.update(completed, force=True)
                if returncode != 0:
                    failures.append(result)
            launch_until_full(executor)
            write_monitor(force_message=True)
    if progress is not None:
        progress.done()
    write_monitor(force_message=True)

    if failures:
        first = failures[0].task
        raise SystemExit(
            f"{len(failures)} standalone shard(s) failed. First failure: "
            f"{first.variation_name} shard {first.spec.index}, logs: {first.stdout_path}, {first.stderr_path}"
        )

    yoda_paths: dict[str, Path] = {}
    for variation_name in variations:
        raw_inputs = [task.raw_yoda_path for task in tasks if task.variation_name == variation_name]
        missing = [path for path in raw_inputs if not path.exists()]
        if missing:
            raise SystemExit(f"Missing raw shard YODA for {variation_name}: {missing[0]}")
        suffix = "" if variation_name == "nominal" else f"_{variation_name.lower()}"
        merged_raw = paths["work_dir"] / f"standalone_fo_{args.setup.lower()}_correlated{suffix}.raw.merged.yoda"
        _progress_message(args, f"[standalone-campaign] merging {len(raw_inputs)} raw YODA shards for {variation_name}")
        _merge_yoda_inputs(
            raw_inputs,
            merged_raw,
            base_dir=paths["base_dir"],
            merge_tool_preference=args.yoda_merge_tool,
            label=f"{variation_name} raw shard YODA",
        )
        clean_yoda = paths["nominal_yoda"] if variation_name == "nominal" else (
            paths["work_dir"] / f"standalone_fo_{args.setup.lower()}_correlated{suffix}.yoda"
        )
        _progress_message(args, f"[standalone-campaign] deriving correlated polarized objects for {variation_name}")
        _derive_correlated_yoda(merged_raw, clean_yoda, args.analysis)
        yoda_paths[variation_name] = clean_yoda

    if args.scale_variations:
        object_count = _write_variation_tagged_yoda(yoda_paths, paths["variation_yoda"], args.analysis)
        _progress_message(
            args,
            f"[standalone-campaign] wrote variation-tagged YODA {paths['variation_yoda']} with {object_count} objects",
        )


def command_standalone_campaign(args: argparse.Namespace) -> int:
    if args.powheg:
        raise SystemExit("POWHEG-selected generation is intentionally unavailable in this Python fixed-order validator.")
    if args.events <= 0:
        raise SystemExit("--events must be positive.")
    if args.jobs <= 0:
        raise SystemExit("--jobs must be positive.")
    if args.shards <= 0:
        raise SystemExit("--shards must be positive.")

    paths = _standalone_campaign_paths(args)
    paths["analysis_dir"].mkdir(parents=True, exist_ok=True)
    paths["work_dir"].mkdir(parents=True, exist_ok=True)
    if not paths["poldis_reference"].exists():
        raise SystemExit(f"Missing POLDIS reference: {paths['poldis_reference']}")

    total_steps = 4
    if _parallel_generation_requested(args):
        _numbered_step(args, 1, total_steps, "run parallel correlated standalone raw shards")
        _run_parallel_standalone_generation(args, paths)
    else:
        _numbered_step(args, 1, total_steps, "run correlated standalone polarized-full")
        polarized_args = argparse.Namespace(
            command="polarized-full",
            setup=args.setup,
            pdf_profile=args.pdf_profile,
            window=args.window,
            q2_min=args.q2_min,
            q2_max=args.q2_max,
            y_min=args.y_min,
            y_max=args.y_max,
            q2_sampling=args.q2_sampling,
            sampler=args.sampler,
            toy_pdfs=args.toy_pdfs,
            seed=args.seed,
            scale_factor=args.scale_factor,
            jetinput=args.jetinput,
            no_progress=args.no_progress,
            powheg=args.powheg,
            events=args.events,
            work_dir=paths["work_dir"],
            output=paths["nominal_yoda"],
            variation_output=paths["variation_yoda"],
            analysis=args.analysis,
            analysis_plugin=paths["analysis_plugin"],
            build_analysis=args.build_analysis,
            scale_variations=args.scale_variations,
            independent_helicities=False,
            polarized_real_weight_mode=args.polarized_real_weight_mode,
        )
        rc = command_polarized_full(polarized_args)
        if rc != 0:
            return rc

    _numbered_step(args, 2, total_steps, "prepare plot-safe Herwig/POLDIS YODAs")
    plot_prep_args = argparse.Namespace(
        herwig=paths["variation_yoda"] if args.scale_variations else paths["nominal_yoda"],
        poldis=paths["poldis_reference"],
        herwig_output=paths["herwig_plot_yoda"],
        poldis_output=paths["poldis_plot_yoda"],
        analysis=args.analysis,
        herwig_scale_envelope=args.scale_variations,
    )
    rc = command_plot_prep(plot_prep_args)
    if rc != 0:
        return rc

    rivet_mkhtml = args.rivet_mkhtml
    if rivet_mkhtml is None:
        resolved = shutil.which("rivet-mkhtml")
        if resolved is None:
            raise SystemExit("Could not find rivet-mkhtml on PATH. Load the Herwig/Rivet module or pass --rivet-mkhtml.")
        rivet_mkhtml = Path(resolved)
    if not paths["rivet_mkhtml_safe"].exists():
        raise SystemExit(f"Missing rivet_mkhtml_safe.py wrapper: {paths['rivet_mkhtml_safe']}")
    if not paths["plot_config"].exists():
        raise SystemExit(f"Missing Rivet plot config: {paths['plot_config']}")

    args.mplconfigdir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(args.mplconfigdir)
    analysis_paths = [paths["work_dir"] / "rivet-plugin", paths["analysis_rivet_dir"]]
    if paths["analysis_plugin"] is not None:
        analysis_paths.insert(0, paths["analysis_plugin"].parent)
    _prepend_env_path(env, "RIVET_ANALYSIS_PATH", analysis_paths)
    _prepend_env_path(env, "RIVET_DATA_PATH", [paths["analysis_rivet_dir"], paths["base_dir"]])
    _prepend_env_path(env, "RIVET_PLOT_PATH", [paths["analysis_rivet_dir"]])

    command = [
        sys.executable,
        str(paths["rivet_mkhtml_safe"]),
        str(rivet_mkhtml),
        f"{paths['herwig_plot_yoda']}:Title={args.herwig_title}",
        f"{paths['poldis_plot_yoda']}:Title={args.poldis_title}",
        "-c",
        str(paths["plot_config"]),
        "--reflabel",
        args.reflabel,
        "--ratiolabel",
        args.ratiolabel,
        "--match",
        args.match,
        "-o",
        str(paths["plots_dir"]),
    ]
    if not args.quiet_mkhtml:
        command.insert(-2, "--verbose")

    _numbered_step(args, 3, total_steps, "make Rivet HTML plots")
    with _timed_step(args, f"run rivet-mkhtml {paths['plots_dir']}"):
        subprocess.run(command, check=True, env=env)

    _numbered_step(args, 4, total_steps, "postprocess Rivet ratio/scale plots")
    with _timed_step(args, f"postprocess plot scripts {paths['plots_dir']}"):
        old_mplconfigdir = os.environ.get("MPLCONFIGDIR")
        os.environ["MPLCONFIGDIR"] = str(args.mplconfigdir)
        try:
            if args.scale_variations:
                patched_count, rerendered_count = rewrite_scale_envelope_plot_scripts(
                    paths["plots_dir"],
                    python_executable=sys.executable,
                    reference_label=args.reflabel,
                )
            else:
                patched_count, rerendered_count = rewrite_no_scale_ratio_plot_scripts(
                    paths["plots_dir"],
                    python_executable=sys.executable,
                )
        finally:
            if old_mplconfigdir is None:
                os.environ.pop("MPLCONFIGDIR", None)
            else:
                os.environ["MPLCONFIGDIR"] = old_mplconfigdir
    _progress_message(
        args,
        f"[standalone-campaign] postprocessed Rivet scripts: patched {patched_count}, rerendered {rerendered_count}",
    )

    index = paths["plots_dir"] / args.analysis / "index.html"
    print(f"Wrote standalone campaign products to {paths['campaign_dir']}")
    print(f"Nominal YODA: {paths['nominal_yoda']}")
    if args.scale_variations:
        print(f"Scale variations YODA: {paths['variation_yoda']}")
    print(f"Plot index: {index}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.command == "generate":
        return command_generate(args)
    if args.command == "integrate":
        return command_integrate(args)
    if args.command == "rivet":
        return command_rivet(args)
    if args.command == "compare":
        return command_compare(args)
    if args.command == "polarize":
        return command_polarize(args)
    if args.command == "polarized-full":
        return command_polarized_full(args)
    if args.command == "polarized-shard":
        return command_polarized_shard(args)
    if args.command == "full":
        return command_full(args)
    if args.command == "plot-prep":
        return command_plot_prep(args)
    if args.command == "standalone-campaign":
        return command_standalone_campaign(args)
    raise SystemExit(f"Unknown command {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
