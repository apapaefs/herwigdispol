#!/usr/bin/env python3
"""Build the polarized-shower paper plots from physical-helicity shards.

The spinphysical campaigns generate the four physical beam-helicity samples
independently, with separate positive- and negative-weight POWHEG runs.  This
script combines those samples at shard-block level, so correlations between
the numerator and denominator of an asymmetry or moment are retained in its
Monte Carlo uncertainty.

Only the observables used in the paper are read.  In particular, the script
does not use the correlated analysis weights intended for differential hard-
process validation: a physical-helicity hard process is needed to retain the
polarization information that initializes the polarized parton shower.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


ANALYSIS = "MC_DIS_PS"
SETUPS = ("SPINCOMP", "SPINHAD")
HELICITIES = ("PP", "PM", "MP", "MM")
ORDERS = ("POSNLO", "NEGNLO")
VARIANTS = (
    "RIVETPS-SPIN",
    "RIVETPS-NOSPIN",
    "RIVETPS-NOSPIN-UNPOL",
)
OBSERVABLES = (
    "pT3OverpT1",
    "DeltaPhiHardJ3",
    "Cos2DeltaPhiHardJ3NumPt3OverpT1",
    "Cos2DeltaPhiHardJ3DenPt3OverpT1",
)
LABELS = {
    "RIVETPS-SPIN": "Full",
    "RIVETPS-NOSPIN": "Born-Only",
    "RIVETPS-NOSPIN-UNPOL": "None",
}
STYLES = {
    "RIVETPS-SPIN": {"color": "#D62728", "marker": "o"},
    "RIVETPS-NOSPIN": {"color": "#1F77B4", "marker": "s"},
    "RIVETPS-NOSPIN-UNPOL": {"color": "#2CA02C", "marker": "^"},
}
SETUP_LABELS = {
    "SPINCOMP": "Before hadronization",
    "SPINHAD": "After hadronization",
}
REFERENCE_VARIANT = "RIVETPS-NOSPIN-UNPOL"
PT_RATIO_REBIN_EDGES = np.asarray(
    (0.0, 2 / 15, 4 / 15, 6 / 15, 8 / 15, 10 / 15, 12 / 15, 1.0),
    dtype=float,
)


def recover_shard_index(spec: Mapping[str, object]) -> int:
    """Recover the logical shard from the tag, including replacement runs.

    A few recovered spinphysical15 records retain the replacement process's
    temporary ``shard_index`` value.  The ``-sNNN`` tag suffix is the stable
    logical identity and is therefore authoritative.
    """

    match = re.search(r"-s(\d+)$", str(spec.get("tag", "")))
    if match is None:
        raise ValueError(f"Cannot recover a shard index from tag {spec.get('tag')!r}")
    return int(match.group(1))


def _object_path(objects: Mapping[str, object], analysis: str, observable: str) -> str:
    for candidate in (f"/{analysis}/{observable}", f"/RAW/{analysis}/{observable}"):
        if candidate in objects:
            return candidate

    matches = []
    for path in objects:
        parts = [part for part in path.split("/") if part]
        if parts and parts[0] == "RAW":
            parts = parts[1:]
        if len(parts) >= 2 and parts[0].split(":", 1)[0] == analysis and parts[1] == observable:
            matches.append(path)
    if not matches:
        raise KeyError(f"No {analysis}/{observable} object in YODA file")
    return sorted(matches)[0]


def _read_yoda_task(task):
    """Worker entry point: read all paper observables from one shard."""

    key, path, analysis, observables = task
    import yoda  # Imported here so numerical unit tests do not require YODA.

    objects = yoda.read(path)
    payload = {}
    for observable in observables:
        histogram = objects[_object_path(objects, analysis, observable)]
        bins = list(histogram.bins())
        if not bins:
            raise ValueError(f"Empty histogram {analysis}/{observable} in {path}")
        edges = [float(bins[0].xMin()), *[float(item.xMax()) for item in bins]]
        values = [
            float(item.val()) if hasattr(item, "val") else float(item.height())
            for item in bins
        ]
        payload[observable] = (edges, values)
    return key, payload


def select_shard_files(
    manifest: Mapping[str, object],
    setups: Sequence[str] = SETUPS,
    variants: Sequence[str] = VARIANTS,
) -> tuple[dict[tuple[str, str, str, str, int], str], int]:
    """Select and validate the physical-helicity production shards."""

    shard_count = int(manifest.get("shards_per_logical_run", 0))
    if shard_count <= 0:
        raise ValueError("Manifest has no positive shards_per_logical_run")

    selected: dict[tuple[str, str, str, str, int], str] = {}
    for record in manifest.get("finished", []):
        spec = record.get("spec", {})
        job = spec.get("job", {})
        setup = str(job.get("setup", ""))
        variant = str(job.get("analysis_variant", ""))
        order = str(job.get("order", ""))
        helicity = str(job.get("helicity", ""))
        if (
            setup not in setups
            or variant not in variants
            or order not in ORDERS
            or helicity not in HELICITIES
        ):
            continue

        shard_index = recover_shard_index(spec)
        key = (setup, variant, order, helicity, shard_index)
        files = record.get("yoda_files", [])
        if len(files) != 1:
            raise ValueError(f"Expected one YODA file for {key}, found {files!r}")
        if key in selected:
            raise ValueError(f"Duplicate logical shard {key}")
        selected[key] = str(files[0])

    expected = len(setups) * len(variants) * len(ORDERS) * len(HELICITIES) * shard_count
    if len(selected) != expected:
        raise ValueError(f"Selected {len(selected)} physical-helicity shards; expected {expected}")

    for key, path in selected.items():
        if not Path(path).is_file():
            raise FileNotFoundError(f"Missing YODA file for {key}: {path}")
    return selected, shard_count


def read_shards(
    selected: Mapping[tuple[str, str, str, str, int], str],
    workers: int,
    analysis: str = ANALYSIS,
    observables: Sequence[str] = OBSERVABLES,
) -> tuple[dict[tuple[str, str, str, str, int, str], np.ndarray], dict[str, np.ndarray]]:
    """Read each YODA file once and return values plus common bin edges."""

    tasks = [(key, path, analysis, tuple(observables)) for key, path in selected.items()]
    values: dict[tuple[str, str, str, str, int, str], np.ndarray] = {}
    edges: dict[str, np.ndarray] = {}
    print(f"Reading {len(tasks)} shard YODAs with {workers} workers ...", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = pool.map(_read_yoda_task, tasks, chunksize=8)
        for index, (key, payload) in enumerate(results, 1):
            for observable, (loaded_edges, loaded_values) in payload.items():
                loaded_edges_array = np.asarray(loaded_edges, dtype=float)
                if observable not in edges:
                    edges[observable] = loaded_edges_array
                elif not np.array_equal(edges[observable], loaded_edges_array):
                    raise ValueError(f"Mismatched {observable} binning in shard {key}")
                values[(*key, observable)] = np.asarray(loaded_values, dtype=float)
            if index % 800 == 0 or index == len(tasks):
                print(f"  read {index}/{len(tasks)}", flush=True)
    return values, edges


def build_physical_blocks(
    values: Mapping[tuple[str, str, str, str, int, str], np.ndarray],
    shard_count: int,
) -> dict[tuple[str, str, str, str], np.ndarray]:
    """Form POSNLO-NEGNLO and then sigma/Delta for every shard block."""

    result: dict[tuple[str, str, str, str], np.ndarray] = {}
    signs = {"PP": 1.0, "PM": -1.0, "MP": -1.0, "MM": 1.0}
    for setup in SETUPS:
        for variant in VARIANTS:
            for observable in OBSERVABLES:
                helicity_blocks = {}
                for helicity in HELICITIES:
                    positive = np.stack(
                        [
                            values[(setup, variant, "POSNLO", helicity, shard, observable)]
                            for shard in range(1, shard_count + 1)
                        ]
                    )
                    negative = np.stack(
                        [
                            values[(setup, variant, "NEGNLO", helicity, shard, observable)]
                            for shard in range(1, shard_count + 1)
                        ]
                    )
                    helicity_blocks[helicity] = positive - negative

                sigma = sum(helicity_blocks.values()) / 4.0
                delta = sum(signs[h] * helicity_blocks[h] for h in HELICITIES) / 4.0
                result[(setup, variant, observable, "sigma")] = sigma
                result[(setup, variant, observable, "delta")] = delta
    return result


def mean_and_error(blocks: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean and standard error from equal-statistics shard blocks."""

    array = np.asarray(blocks, dtype=float)
    if array.shape[0] < 2:
        raise ValueError("At least two shard blocks are required")
    return np.mean(array, axis=0), np.std(array, axis=0, ddof=1) / math.sqrt(array.shape[0])


def ratio_of_means(
    numerator_blocks: np.ndarray,
    denominator_blocks: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Ratio of block means with numerator-denominator covariance retained."""

    numerator = np.asarray(numerator_blocks, dtype=float)
    denominator = np.asarray(denominator_blocks, dtype=float)
    if numerator.shape != denominator.shape:
        raise ValueError("Ratio block arrays must have the same shape")
    denominator_mean = np.mean(denominator, axis=0)
    if np.any(denominator_mean == 0.0):
        raise ZeroDivisionError("A ratio denominator has zero block mean")
    ratio = np.mean(numerator, axis=0) / denominator_mean
    influence = (numerator - ratio * denominator) / denominator_mean
    error = np.std(influence, axis=0, ddof=1) / math.sqrt(numerator.shape[0])
    return ratio, error, influence


def independent_ratio(
    numerator: np.ndarray,
    numerator_error: np.ndarray,
    denominator: np.ndarray,
    denominator_error: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Ratio of separately generated configuration estimates."""

    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    ratio = numerator / denominator
    error = np.hypot(
        np.asarray(numerator_error, dtype=float) / denominator,
        numerator * np.asarray(denominator_error, dtype=float) / denominator**2,
    )
    return ratio, error


def independent_difference(
    value: np.ndarray,
    error: np.ndarray,
    reference: np.ndarray,
    reference_error: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Difference of separately generated configuration estimates."""

    return np.asarray(value) - np.asarray(reference), np.hypot(error, reference_error)


def integrate_rebin(blocks: np.ndarray, old_edges: np.ndarray, new_edges: np.ndarray) -> np.ndarray:
    """Integrate histogram heights into a new, possibly coarser binning."""

    old_edges = np.asarray(old_edges, dtype=float)
    new_edges = np.asarray(new_edges, dtype=float)
    array = np.asarray(blocks, dtype=float)
    if array.shape[-1] != len(old_edges) - 1:
        raise ValueError("Histogram array and old edges are inconsistent")

    rebinned = []
    for low, high in zip(new_edges[:-1], new_edges[1:]):
        overlaps = np.maximum(
            0.0,
            np.minimum(old_edges[1:], high) - np.maximum(old_edges[:-1], low),
        )
        if not np.isclose(np.sum(overlaps), high - low, rtol=0.0, atol=1e-12):
            raise ValueError(f"New bin [{low}, {high}] is not covered by the old binning")
        rebinned.append(array @ overlaps)
    return np.stack(rebinned, axis=-1)


def cumulative_tail(blocks: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Return width-weighted integrals above every lower bin edge."""

    weighted = np.asarray(blocks, dtype=float) * np.diff(np.asarray(edges, dtype=float))
    return np.flip(np.cumsum(np.flip(weighted, axis=-1), axis=-1), axis=-1)


def _curve_payload(values, errors, x, xerr=None):
    payload = {
        "x": np.asarray(x, dtype=float),
        "values": np.asarray(values, dtype=float),
        "errors": np.asarray(errors, dtype=float),
    }
    if xerr is not None:
        payload["xerr"] = np.asarray(xerr, dtype=float)
    return payload


def build_plot_data(blocks, edges):
    """Construct the five paper observables for both shower stages."""

    output = {}
    pt_edges = edges["pT3OverpT1"]
    phi_edges = edges["DeltaPhiHardJ3"]
    phi_x = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    phi_xerr = 0.5 * np.diff(phi_edges)
    moment_edges = edges["Cos2DeltaPhiHardJ3DenPt3OverpT1"]

    for setup in SETUPS:
        setup_data = {
            "ALLpT3OverpT1": {},
            "DeltaPhiHardJ3": {},
            "DDeltaPhiHardJ3": {},
            "Cos2DeltaPhiHardJ3Cumulative": {},
            "ALLCos2DeltaPhiHardJ3Cumulative": {},
        }
        for variant in VARIANTS:
            delta_pt = integrate_rebin(
                blocks[(setup, variant, "pT3OverpT1", "delta")],
                pt_edges,
                PT_RATIO_REBIN_EDGES,
            )
            sigma_pt = integrate_rebin(
                blocks[(setup, variant, "pT3OverpT1", "sigma")],
                pt_edges,
                PT_RATIO_REBIN_EDGES,
            )
            all_value, all_error, _ = ratio_of_means(delta_pt, sigma_pt)
            setup_data["ALLpT3OverpT1"][variant] = _curve_payload(
                all_value,
                all_error,
                0.5 * (PT_RATIO_REBIN_EDGES[:-1] + PT_RATIO_REBIN_EDGES[1:]),
                0.5 * np.diff(PT_RATIO_REBIN_EDGES),
            )

            sigma_phi = blocks[(setup, variant, "DeltaPhiHardJ3", "sigma")]
            delta_phi = blocks[(setup, variant, "DeltaPhiHardJ3", "delta")]
            sigma_value, sigma_error = mean_and_error(sigma_phi)
            delta_value, delta_error = mean_and_error(delta_phi)
            setup_data["DeltaPhiHardJ3"][variant] = _curve_payload(
                sigma_value, sigma_error, phi_x, phi_xerr
            )
            setup_data["DDeltaPhiHardJ3"][variant] = _curve_payload(
                delta_value, delta_error, phi_x, phi_xerr
            )

            sigma_num = cumulative_tail(
                blocks[(setup, variant, "Cos2DeltaPhiHardJ3NumPt3OverpT1", "sigma")],
                moment_edges,
            )
            delta_num = cumulative_tail(
                blocks[(setup, variant, "Cos2DeltaPhiHardJ3NumPt3OverpT1", "delta")],
                moment_edges,
            )
            sigma_den = cumulative_tail(
                blocks[(setup, variant, "Cos2DeltaPhiHardJ3DenPt3OverpT1", "sigma")],
                moment_edges,
            )
            c2_value, c2_error, _ = ratio_of_means(sigma_num, sigma_den)
            dc2_value, dc2_error, _ = ratio_of_means(delta_num, sigma_den)
            # Successive cumulative thresholds are strongly correlated.  Show
            # every second threshold to keep the visual comparison readable.
            shown = np.arange(0, len(moment_edges) - 1, 2)
            setup_data["Cos2DeltaPhiHardJ3Cumulative"][variant] = _curve_payload(
                c2_value[shown], c2_error[shown], moment_edges[:-1][shown]
            )
            setup_data["ALLCos2DeltaPhiHardJ3Cumulative"][variant] = _curve_payload(
                dc2_value[shown], dc2_error[shown], moment_edges[:-1][shown]
            )
        output[setup] = setup_data
    return output


def _matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10.0,
            "axes.labelsize": 10.0,
            "axes.titlesize": 10.5,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.25,
            "savefig.bbox": "tight",
        }
    )
    return plt


def _draw_curve(axis, curve, variant, label=True):
    style = STYLES[variant]
    axis.errorbar(
        curve["x"],
        curve["values"],
        yerr=curve["errors"],
        xerr=curve.get("xerr"),
        label=LABELS[variant] if label else None,
        color=style["color"],
        marker=style["marker"],
        linestyle="none",
        markersize=3.6,
        markerfacecolor="white",
        markeredgewidth=0.8,
        capsize=1.7,
        elinewidth=0.75,
        zorder=4,
    )


def _comparison_curve(curve, reference, mode):
    if mode == "ratio":
        values, errors = independent_ratio(
            curve["values"], curve["errors"], reference["values"], reference["errors"]
        )
    elif mode == "difference":
        values, errors = independent_difference(
            curve["values"], curve["errors"], reference["values"], reference["errors"]
        )
    else:
        raise ValueError(f"Unknown comparison mode {mode!r}")
    return _curve_payload(values, errors, curve["x"], curve.get("xerr"))


def _finite_extent(curves: Iterable[Mapping[str, np.ndarray]], include_zero=False):
    lows = []
    highs = []
    for curve in curves:
        values = np.asarray(curve["values"])
        errors = np.asarray(curve["errors"])
        finite = np.isfinite(values) & np.isfinite(errors)
        lows.extend((values[finite] - errors[finite]).tolist())
        highs.extend((values[finite] + errors[finite]).tolist())
    if include_zero:
        lows.append(0.0)
        highs.append(0.0)
    low, high = min(lows), max(highs)
    span = high - low
    margin = 0.12 * (span if span > 0 else max(abs(low), 1.0))
    return low - margin, high + margin


def plot_observable(output_dir, setup, observable, curves, png=False):
    """Draw one paper panel with a reference comparison underneath."""

    plt = _matplotlib()
    metadata = {
        "ALLpT3OverpT1": (
            r"Double-spin asymmetry $A_{LL}(r_3)$",
            r"$r_3=p^{\mathrm{jet3}}_{T,B}/p^{\mathrm{jet1}}_{T,B}$",
            r"$A_{LL}$",
            "ratio",
        ),
        "DeltaPhiHardJ3": (
            "Hard-plane 3rd-jet azimuth",
            r"$\Delta\phi_{\mathrm{hard},3}$",
            r"$\mathrm{d}\sigma/\mathrm{d}\Delta\phi_{\mathrm{hard},3}$ [pb]",
            "ratio",
        ),
        "DDeltaPhiHardJ3": (
            "Polarized hard-plane 3rd-jet azimuth",
            r"$\Delta\phi_{\mathrm{hard},3}$",
            r"$\mathrm{d}\Delta\sigma/\mathrm{d}\Delta\phi_{\mathrm{hard},3}$ [pb]",
            "ratio",
        ),
        "Cos2DeltaPhiHardJ3Cumulative": (
            r"Hard-plane $C_2$ tail moment",
            r"$r_{\mathrm{cut}}$",
            r"$C_{2}^{\mathrm{hard}}(r_{\mathrm{cut}})$",
            "difference",
        ),
        "ALLCos2DeltaPhiHardJ3Cumulative": (
            r"Spin-dependent hard-plane $C_2$ tail moment",
            r"$r_{\mathrm{cut}}$",
            r"$\Delta C_{2}^{\mathrm{hard}}(r_{\mathrm{cut}})$",
            "difference",
        ),
    }
    title, xlabel, ylabel, comparison_mode = metadata[observable]
    figure, (upper, lower) = plt.subplots(
        2,
        1,
        figsize=(4.65, 4.35),
        gridspec_kw={"height_ratios": (3.2, 1.0), "hspace": 0.06},
        sharex=True,
    )

    for variant in VARIANTS:
        _draw_curve(upper, curves[variant], variant)
    upper.set_title(f"{title}\n{SETUP_LABELS[setup]}", pad=5)
    upper.set_ylabel(ylabel)
    upper.legend(loc="best", frameon=False, ncol=1)
    upper.tick_params(labelbottom=False)

    reference = curves[REFERENCE_VARIANT]
    baseline = 1.0 if comparison_mode == "ratio" else 0.0
    lower.axhline(baseline, color=STYLES[REFERENCE_VARIANT]["color"], linestyle=":", linewidth=1.0)
    for variant in VARIANTS[:-1]:
        _draw_curve(lower, _comparison_curve(curves[variant], reference, comparison_mode), variant, label=False)
    lower.set_xlabel(xlabel)
    lower.set_ylabel(
        "ratio to None" if comparison_mode == "ratio" else "difference\nfrom None",
        fontsize=8.0,
        labelpad=7,
    )

    for axis in (upper, lower):
        axis.set_axisbelow(True)
        axis.minorticks_on()
        axis.tick_params(which="major", length=4.0, width=0.8)
        axis.tick_params(which="minor", length=2.2, width=0.6)
        axis.grid(which="major", color="0.86", linewidth=0.55)
        axis.grid(which="minor", color="0.92", linewidth=0.35, linestyle=":")

    if observable == "ALLpT3OverpT1":
        upper.set_xlim(0.24, 1.02)
        upper.set_ylim(0.020, 0.045)
        lower.set_ylim(0.80, 1.15)
    elif observable in {"DeltaPhiHardJ3", "DDeltaPhiHardJ3"}:
        upper.set_xlim(-0.08, math.pi + 0.08)
        lower.set_ylim(0.80, 1.20)
        upper.set_ylim(*_finite_extent(curves.values(), include_zero=False))
    else:
        upper.set_xlim(-0.04, 1.02)
        upper.axhline(0.0, color="0.55", linewidth=0.7, zorder=1)
        upper.set_ylim(*_finite_extent(curves.values(), include_zero=True))
        differences = [
            _comparison_curve(curves[variant], reference, "difference")
            for variant in VARIANTS[:-1]
        ]
        low, high = _finite_extent(differences, include_zero=True)
        extent = max(abs(low), abs(high))
        lower.set_ylim(-extent, extent)

    upper.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)
    lower.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useMathText=True)

    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{setup.lower()}_{observable}"
    pdf_path = output_dir / f"{stem}.pdf"
    figure.savefig(
        pdf_path,
        metadata={
            "Title": f"{title} ({SETUP_LABELS[setup].lower()})",
            "Creator": Path(__file__).name,
        },
    )
    if png:
        figure.savefig(output_dir / f"{stem}.png", dpi=180)
    plt.close(figure)
    return pdf_path


def _json_curve(curve):
    return {
        key: np.asarray(value, dtype=float).tolist()
        for key, value in curve.items()
        if key in {"x", "xerr", "values", "errors"}
    }


def write_summary(path, manifest_path, manifest, plot_data, plot_paths, shard_count):
    """Record inputs, uncertainty convention, and numerical plot values."""

    payload = {
        "campaign": manifest.get("tag"),
        "manifest": str(manifest_path.resolve()),
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "physical_helicity_samples": list(HELICITIES),
        "nlo_combination": "POSNLO minus NEGNLO for each helicity",
        "shards_per_logical_run": shard_count,
        "uncertainty": (
            "equal-statistics shard-block standard errors; ratio-of-means influence "
            "functions retain numerator-denominator and cumulative-bin covariance"
        ),
        "configuration_comparisons": "independent-error propagation",
        "pt_ratio_rebin_edges": PT_RATIO_REBIN_EDGES.tolist(),
        "cumulative_threshold_note": "successive thresholds are statistically correlated",
        "plot_files": [str(path.resolve()) for path in plot_paths],
        "data": {},
    }
    for setup, observables in plot_data.items():
        payload["data"][setup] = {
            observable: {variant: _json_curve(curve) for variant, curve in curves.items()}
            for observable, curves in observables.items()
        }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="Completed spinphysical manifest.json")
    parser.add_argument("--output-dir", type=Path, required=True, help="Destination for paper PDFs")
    parser.add_argument("--workers", type=int, default=16, help="Parallel YODA readers")
    parser.add_argument("--png", action="store_true", help="Also write PNG previews")
    parser.add_argument(
        "--summary",
        type=Path,
        help="JSON summary path (default: OUTPUT_DIR/spinphysical_paper_plots.json)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    manifest = json.loads(args.manifest.read_text())
    selected, shard_count = select_shard_files(manifest)
    values, edges = read_shards(selected, args.workers)
    blocks = build_physical_blocks(values, shard_count)
    plot_data = build_plot_data(blocks, edges)

    plot_paths = []
    for setup in SETUPS:
        for observable, curves in plot_data[setup].items():
            path = plot_observable(args.output_dir, setup, observable, curves, png=args.png)
            plot_paths.append(path)
            print(f"Wrote {path}")

    summary_path = args.summary or args.output_dir / "spinphysical_paper_plots.json"
    write_summary(summary_path, args.manifest, manifest, plot_data, plot_paths, shard_count)
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
