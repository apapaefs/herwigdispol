#!/usr/bin/env python3.10
"""
Compare alpha_s(Q^2) from multiple LHAPDF PDF sets over the active DISPOL range.

By default, this compares the polarized Herwig-side set and the unpolarized
reference set, and adds the Herwig Matchbox NLOAlphaS curve used by the DIS
cards, normalized through its input value alpha_s(M_Z) = 0.118:

  - BDSSV24-NNLO
  - PDF4LHC15_nnlo_100_pdfas
  - Herwig Matchbox NLOAlphaS (input alpha_s(MZ)=0.118)

The POLDIS-internal DSSV replica can still be compared later by passing an
external alpha_s(Q^2) export into a separate plotting step, but it is not
directly loadable from standalone Python LHAPDF in this environment.

Usage:
  python3.10 compare_pdf_alphas.py
  python3.10 compare_pdf_alphas.py --reference herwig
  python3.10 compare_pdf_alphas.py --input-alpha-s 0.118 --show
"""

from __future__ import annotations

import argparse
import math
import os
import tempfile
from pathlib import Path
from typing import Callable, Dict, List, Sequence


DEFAULT_LHAPDF_DATA_PATH = "/opt/homebrew/Cellar/lhapdf/6.5.4/share/LHAPDF"
DEFAULT_PDFS = [
    "BDSSV24-NNLO",
    "PDF4LHC15_nnlo_100_pdfas",
]
DEFAULT_REFERENCE = "BDSSV24-NNLO"
DEFAULT_INPUT_ALPHA_S = 0.118
DEFAULT_INPUT_SCALE = 91.188
DEFAULT_Q2_MIN = 49.0
DEFAULT_Q2_MAX = 2500.0
DEFAULT_NUM_POINTS = 300
DEFAULT_Q2_VALUES = [49.0, 100.0, 500.0, 1500.0, 2500.0]
HERWIG_REFERENCE_TOKEN = "herwig"


def parse_q2_values(text: str) -> List[float]:
    values: List[float] = []
    for chunk in text.split(","):
        item = chunk.strip()
        if not item:
            continue
        try:
            value = float(item)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"Could not parse Q^2 value '{item}'") from exc
        if value <= 0.0:
            raise argparse.ArgumentTypeError("Q^2 values must be positive")
        values.append(value)
    if not values:
        raise argparse.ArgumentTypeError("Provide at least one comma-separated Q^2 value")
    return values


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdf",
        action="append",
        dest="pdfs",
        help="Repeatable LHAPDF set name. If omitted, uses the DISPOL defaults.",
    )
    parser.add_argument(
        "--reference",
        help=(
            "Reference curve for the ratio panel. Defaults to "
            f"{DEFAULT_REFERENCE} when present, otherwise the first requested curve. "
            f"Use '{HERWIG_REFERENCE_TOKEN}' to select the Herwig-style curve."
        ),
    )
    parser.add_argument(
        "--input-alpha-s",
        type=float,
        default=DEFAULT_INPUT_ALPHA_S,
        help="Target Matchbox NLOAlphaS input_alpha_s at the input scale.",
    )
    parser.add_argument(
        "--input-scale",
        type=float,
        default=DEFAULT_INPUT_SCALE,
        help="Matchbox NLOAlphaS input_scale in GeV.",
    )
    parser.add_argument(
        "--alpha-s-mz",
        dest="input_alpha_s",
        type=float,
        help="Backward-compatible alias for --input-alpha-s.",
    )
    parser.add_argument(
        "--mz",
        dest="input_scale",
        type=float,
        help="Backward-compatible alias for --input-scale.",
    )
    parser.add_argument(
        "--no-herwig-curve",
        action="store_true",
        help="Disable the extra Herwig Matchbox NLOAlphaS comparison curve.",
    )
    parser.add_argument("--q2-min", type=float, default=DEFAULT_Q2_MIN, help="Minimum Q^2 in GeV^2.")
    parser.add_argument("--q2-max", type=float, default=DEFAULT_Q2_MAX, help="Maximum Q^2 in GeV^2.")
    parser.add_argument(
        "--num-points",
        type=int,
        default=DEFAULT_NUM_POINTS,
        help="Number of log-spaced Q^2 points used for the curves.",
    )
    parser.add_argument(
        "--q2-values",
        type=parse_q2_values,
        default=list(DEFAULT_Q2_VALUES),
        help="Comma-separated Q^2 values for the printed table.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output PNG path. Defaults to compare_pdf_alphas.png next to this script.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure interactively in addition to saving it.",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if args.q2_min <= 0.0 or args.q2_max <= 0.0:
        raise SystemExit("Q^2 bounds must be positive.")
    if args.q2_min >= args.q2_max:
        raise SystemExit("Require q2-min < q2-max.")
    if args.num_points < 2:
        raise SystemExit("Require at least 2 grid points.")
    if args.input_alpha_s <= 0.0:
        raise SystemExit("Require input-alpha-s > 0.")
    if args.input_scale <= 0.0:
        raise SystemExit("Require input-scale > 0.")


def resolve_pdf_names(args: argparse.Namespace) -> List[str]:
    pdf_names = list(args.pdfs) if args.pdfs else list(DEFAULT_PDFS)
    seen = set()
    deduped: List[str] = []
    for name in pdf_names:
        if name not in seen:
            deduped.append(name)
            seen.add(name)
    if len(deduped) < 2:
        raise SystemExit("Provide at least two distinct PDF sets for comparison.")
    return deduped


def resolve_reference_name(
    series_names: Sequence[str], requested: str | None, herwig_name: str | None = None
) -> str:
    if requested:
        if herwig_name and requested.strip().lower() == HERWIG_REFERENCE_TOKEN:
            return herwig_name
        if requested not in series_names:
            raise SystemExit(
                f"Reference curve '{requested}' is not in the requested comparison list: "
                f"{', '.join(series_names)}"
            )
        return requested
    if DEFAULT_REFERENCE in series_names:
        return DEFAULT_REFERENCE
    return series_names[0]


def import_runtime_dependencies(show: bool):
    if "LHAPDF_DATA_PATH" not in os.environ:
        os.environ["LHAPDF_DATA_PATH"] = DEFAULT_LHAPDF_DATA_PATH
    if "MPLCONFIGDIR" not in os.environ:
        mpl_cache_dir = Path(tempfile.gettempdir()) / "herwigpol-matplotlib"
        mpl_cache_dir.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(mpl_cache_dir)
    if "XDG_CACHE_HOME" not in os.environ:
        xdg_cache_dir = Path(tempfile.gettempdir()) / "herwigpol-xdg-cache"
        xdg_cache_dir.mkdir(parents=True, exist_ok=True)
        os.environ["XDG_CACHE_HOME"] = str(xdg_cache_dir)

    try:
        import numpy as np
    except ImportError as exc:
        raise SystemExit("Could not import numpy. Use a Python environment with numpy installed.") from exc

    try:
        import matplotlib
    except ImportError as exc:
        raise SystemExit("Could not import matplotlib. Use a Python environment with matplotlib installed.") from exc

    if not show:
        matplotlib.use("Agg")

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("Could not import matplotlib.pyplot.") from exc

    try:
        import lhapdf
    except ImportError as exc:
        raise SystemExit(
            "Could not import the LHAPDF Python bindings. Use a Python with `lhapdf` available "
            "and set LHAPDF_DATA_PATH if needed."
        ) from exc

    return np, plt, lhapdf


def load_pdfs(lhapdf, pdf_names: Sequence[str]) -> Dict[str, object]:
    pdfs: Dict[str, object] = {}
    for name in pdf_names:
        try:
            pdfs[name] = lhapdf.mkPDF(name, 0)
        except Exception as exc:  # pragma: no cover - exception type depends on LHAPDF build
            raise SystemExit(
                f"Could not load LHAPDF set '{name}'. Check that it is installed under "
                f"LHAPDF_DATA_PATH={os.environ.get('LHAPDF_DATA_PATH', '<unset>')}."
            ) from exc
    return pdfs


def evaluate_alphas(np, pdfs: Dict[str, object], q2_grid: Sequence[float]) -> Dict[str, object]:
    values: Dict[str, object] = {}
    for name, pdf in pdfs.items():
        values[name] = np.array([pdf.alphasQ2(float(q2)) for q2 in q2_grid], dtype=float)
    return values


def pdf_value_at_q2(pdf, q2: float) -> float:
    return float(pdf.alphasQ2(float(q2)))


def herwig_series_name(input_alpha_s: float) -> str:
    return f"Herwig Matchbox NLOAlphaS (input alpha_s(MZ)={input_alpha_s:.3f})"


def matchbox_beta_coefficients(nf: int) -> tuple[float, float]:
    beta0 = (33.0 - 2.0 * nf) / (12.0 * math.pi)
    beta1 = (153.0 - 19.0 * nf) / (24.0 * math.pi * math.pi)
    return beta0, beta1


def bisect_root(function: Callable[[float], float], low: float, high: float, label: str) -> float:
    f_low = function(low)
    f_high = function(high)
    if f_low == 0.0:
        return low
    if f_high == 0.0:
        return high
    if f_low * f_high > 0.0:
        raise ValueError(
            f"Could not bracket root for {label}: f({low:.6g})={f_low:.6g}, f({high:.6g})={f_high:.6g}"
        )

    for _ in range(200):
        mid = 0.5 * (low + high)
        f_mid = function(mid)
        if f_mid == 0.0:
            return mid
        if f_low * f_mid <= 0.0:
            high = mid
            f_high = f_mid
        else:
            low = mid
            f_low = f_mid

    return 0.5 * (low + high)


def matchbox_large_scale_alpha(scale_q2: float, lambda2: float, nf: int, two_largeq_terms: bool = True) -> float:
    if scale_q2 <= lambda2:
        raise ValueError(
            f"Matchbox NLOAlphaS emulation requires Q^2 > lambda^2. Got Q^2={scale_q2:.6g}, lambda^2={lambda2:.6g}."
        )
    beta0, beta1 = matchbox_beta_coefficients(nf)
    slog = math.log(scale_q2 / lambda2)
    leading = 1.0 / (beta0 * slog)
    result = leading * (1.0 - (beta1 / (beta0 * beta0)) * math.log(slog) / slog)
    if two_largeq_terms:
        result += leading * (beta1 / (beta0 * beta0 * slog)) ** 2 * ((math.log(slog) - 0.5) ** 2 - 1.25)
    return result


def default_matchbox_quark_masses2() -> List[float]:
    masses = [0.0, 0.005, 0.0023, 0.095, 1.25, 4.2, 174.2]
    masses2 = [mass * mass for mass in masses]
    if masses2[1] > masses2[2]:
        masses2[1], masses2[2] = masses2[2], masses2[1]
    return masses2


def active_flavours(scale_q2: float, quark_masses2: Sequence[float], max_active_flavours: int) -> int:
    active = 0
    if scale_q2 > 0.0:
        while quark_masses2[active] < scale_q2:
            active += 1
            if active == max_active_flavours + 1:
                break
        active -= 1
    return active


def solve_lambda2_for_input_alpha(input_alpha_s: float, input_scale: float, nf: int) -> float:
    input_q2 = float(input_scale) * float(input_scale)

    def equation(lambda2: float) -> float:
        return matchbox_large_scale_alpha(input_q2, lambda2, nf) - input_alpha_s

    return bisect_root(equation, 1.0e-6, 1.0, f"lambda^2 for nf={nf} at input scale")


def solve_matched_lambda2(target_alpha: float, threshold_q2: float, nf: int) -> float:
    def equation(lambda2: float) -> float:
        return matchbox_large_scale_alpha(threshold_q2, lambda2, nf) - target_alpha

    return bisect_root(equation, 1.0e-6, 1.0, f"matched lambda^2 for nf={nf}")


def solve_matchbox_lambda_squared(
    input_alpha_s: float,
    input_scale: float,
    min_active_flavours: int = 3,
    max_active_flavours: int = 6,
) -> List[float]:
    quark_masses2 = default_matchbox_quark_masses2()
    lambdas2 = [0.0] * 7
    input_q2 = float(input_scale) * float(input_scale)
    active_at_input = active_flavours(input_q2, quark_masses2, max_active_flavours)

    lambdas2[active_at_input] = solve_lambda2_for_input_alpha(input_alpha_s, input_scale, active_at_input)

    below = active_at_input
    while below > min_active_flavours:
        threshold_q2 = quark_masses2[below]
        target_alpha = matchbox_large_scale_alpha(threshold_q2, lambdas2[below], below)
        lambdas2[below - 1] = solve_matched_lambda2(target_alpha, threshold_q2, below - 1)
        below -= 1

    above = active_at_input
    while above < max_active_flavours:
        threshold_q2 = quark_masses2[above + 1]
        target_alpha = matchbox_large_scale_alpha(threshold_q2, lambdas2[above], above)
        lambdas2[above + 1] = solve_matched_lambda2(target_alpha, threshold_q2, above + 1)
        above += 1

    for flavor in range(min_active_flavours):
        lambdas2[flavor] = lambdas2[min_active_flavours]
    for flavor in range(max_active_flavours + 1, 7):
        lambdas2[flavor] = lambdas2[max_active_flavours]

    return lambdas2


def matchbox_running_alpha(
    q2: float,
    lambda_squared: Sequence[float],
    quark_masses2: Sequence[float],
    min_active_flavours: int = 3,
    max_active_flavours: int = 6,
) -> float:
    nf = active_flavours(float(q2), quark_masses2, max_active_flavours)
    nf = max(min_active_flavours, min(max_active_flavours, nf))
    return matchbox_large_scale_alpha(float(q2), lambda_squared[nf], nf)


def add_herwig_curve(
    np,
    series_names: List[str],
    alpha_values: Dict[str, object],
    point_evaluators: Dict[str, Callable[[float], float]],
    q2_grid: Sequence[float],
    input_alpha_s: float,
    input_scale: float,
) -> tuple[str, List[float]]:
    name = herwig_series_name(input_alpha_s)
    quark_masses2 = default_matchbox_quark_masses2()
    lambda_squared = solve_matchbox_lambda_squared(input_alpha_s, input_scale)
    alpha_values[name] = np.array(
        [matchbox_running_alpha(float(q2), lambda_squared, quark_masses2) for q2 in q2_grid],
        dtype=float,
    )
    point_evaluators[name] = (
        lambda q2, solved_lambdas=lambda_squared, masses2=quark_masses2: matchbox_running_alpha(
            float(q2), solved_lambdas, masses2
        )
    )
    series_names.append(name)
    return name, lambda_squared


def format_q2(value: float) -> str:
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.3f}"


def format_value(value: float) -> str:
    return f"{value:.8f}"


def render_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    widths = [len(header) for header in headers]
    for row in rows:
        for idx, cell in enumerate(row):
            widths[idx] = max(widths[idx], len(cell))

    def render_row(row: Sequence[str]) -> str:
        cells = [cell.ljust(widths[idx]) for idx, cell in enumerate(row)]
        return "| " + " | ".join(cells) + " |"

    separator = "+-" + "-+-".join("-" * width for width in widths) + "-+"
    lines = [separator, render_row(headers), separator]
    lines.extend(render_row(row) for row in rows)
    lines.append(separator)
    return "\n".join(lines)


def build_table_rows(
    point_evaluators: Dict[str, Callable[[float], float]],
    series_names: Sequence[str],
    reference_name: str,
    q2_values: Sequence[float],
) -> List[List[str]]:
    rows: List[List[str]] = []
    for q2 in q2_values:
        sampled = {name: point_evaluators[name](q2) for name in series_names}
        reference_value = sampled[reference_name]
        row = [format_q2(q2)]
        for name in series_names:
            row.append(format_value(sampled[name]))
        for name in series_names:
            if name == reference_name:
                continue
            value = sampled[name]
            row.append(format_value(value / reference_value))
            row.append(f"{value - reference_value:+.8f}")
        rows.append(row)
    return rows


def make_plot(np, plt, q2_grid, alpha_values, series_names, reference_name: str, output_path: Path) -> None:
    fig, (ax_alpha, ax_ratio) = plt.subplots(
        2,
        1,
        figsize=(8.5, 8.0),
        sharex=True,
        gridspec_kw={"height_ratios": [2.0, 1.1]},
    )

    for name in series_names:
        linestyle = "--" if name.startswith("Herwig Matchbox NLOAlphaS") else "-"
        ax_alpha.plot(q2_grid, alpha_values[name], linewidth=2.0, linestyle=linestyle, label=name)
    ax_alpha.set_xscale("log")
    ax_alpha.set_ylabel(r"$\alpha_s(Q^2)$")
    ax_alpha.set_title(r"PDF $\alpha_s(Q^2)$ comparison over $49 \leq Q^2 \leq 2500\ \mathrm{GeV}^2$")
    ax_alpha.grid(True, which="both", alpha=0.3)
    ax_alpha.legend()

    reference_values = alpha_values[reference_name]
    for name in series_names:
        ratio = alpha_values[name] / reference_values
        linestyle = "--" if name.startswith("Herwig Matchbox NLOAlphaS") else "-"
        ax_ratio.plot(q2_grid, ratio, linewidth=2.0, linestyle=linestyle, label=f"{name} / {reference_name}")
    ax_ratio.axhline(1.0, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    ax_ratio.set_xscale("log")
    ax_ratio.set_xlabel(r"$Q^2\ [\mathrm{GeV}^2]$")
    ax_ratio.set_ylabel("ratio")
    ax_ratio.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    validate_args(args)
    pdf_names = resolve_pdf_names(args)
    output_path = args.output or (Path(__file__).resolve().parent / "compare_pdf_alphas.png")

    np, plt, lhapdf = import_runtime_dependencies(show=args.show)
    pdfs = load_pdfs(lhapdf, pdf_names)

    q2_grid = np.logspace(np.log10(args.q2_min), np.log10(args.q2_max), args.num_points)
    alpha_values = evaluate_alphas(np, pdfs, q2_grid)
    series_names = list(pdf_names)
    point_evaluators: Dict[str, Callable[[float], float]] = {
        name: (lambda q2, pdf=pdfs[name]: pdf_value_at_q2(pdf, q2)) for name in pdf_names
    }
    herwig_name = None
    herwig_lambdas2 = None
    if not args.no_herwig_curve:
        herwig_name, herwig_lambdas2 = add_herwig_curve(
            np,
            series_names,
            alpha_values,
            point_evaluators,
            q2_grid,
            args.input_alpha_s,
            args.input_scale,
        )
    reference_name = resolve_reference_name(series_names, args.reference, herwig_name)

    headers = [r"Q^2 [GeV^2]"] + [f"alpha_s({name})" for name in series_names]
    for name in series_names:
        if name == reference_name:
            continue
        headers.append(f"{name} / {reference_name}")
        headers.append(f"{name} - {reference_name}")

    rows = build_table_rows(point_evaluators, series_names, reference_name, args.q2_values)

    print("Using LHAPDF_DATA_PATH =", os.environ.get("LHAPDF_DATA_PATH", "<unset>"))
    print("PDF sets:", ", ".join(pdf_names))
    if herwig_name:
        print(f"Herwig-style curve: {herwig_name}")
        print(f"  input_scale = {args.input_scale:.4f} GeV, input_alpha_s = {args.input_alpha_s:.3f}")
        print("  uses matchbox::nlo_alpha_s with large_scale evaluation and threshold matching")
        print(
            "  solved lambdas [GeV] for nf=3,4,5,6 = "
            + ", ".join(f"{math.sqrt(herwig_lambdas2[nf]):.6f}" for nf in range(3, 7))
        )
    print("Reference:", reference_name)
    print(f"Q^2 range: {args.q2_min:g} to {args.q2_max:g} GeV^2")
    print()
    print(render_table(headers, rows))

    make_plot(np, plt, q2_grid, alpha_values, series_names, reference_name, output_path)
    print()
    print(f"Saved plot: {output_path}")

    if args.show:
        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
