#!/usr/bin/env python3
"""Dump selected MC_DIS_BREIT/MC_DIS_PS histogram bins from analyzed or campaign-resolved YODAs."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import yoda

from analyze_DIS_polarized import _resolve_input_path, build_dis_polarized_objects
from run_validation_campaign import find_order_yoda_inputs


def iter_bins(obj) -> Iterable[Tuple[int, float, float, float, float]]:
    for index, bin_obj in enumerate(obj.bins()):
        if hasattr(bin_obj, "val"):
            value = float(bin_obj.val())
            try:
                error = float(bin_obj.errAvg())
            except Exception:
                try:
                    error = 0.5 * (
                        abs(float(bin_obj.errMinus())) + abs(float(bin_obj.errPlus()))
                    )
                except Exception:
                    error = 0.0
        else:
            value = float(bin_obj.height())
            error = float(bin_obj.heightErr()) if hasattr(bin_obj, "heightErr") else 0.0
        yield (
            index,
            float(bin_obj.xMin()),
            float(bin_obj.xMax()),
            value,
            error,
        )


def dump_histogram_from_aos(aos: dict, source: str, analysis: str, label: str, xmin: float | None, xmax: float | None) -> None:
    yoda_path = _resolve_input_path(aos, analysis, label)
    hist = aos[yoda_path]

    print(f"# source={source}")
    print(f"# hist={yoda_path}")
    print("bin\txmin\txmax\theight\terror")
    for index, xlo, xhi, height, error in iter_bins(hist):
        if xmin is not None and xhi <= xmin:
            continue
        if xmax is not None and xlo >= xmax:
            continue
        print(f"{index}\t{xlo:.12g}\t{xhi:.12g}\t{height:.12g}\t{error:.12g}")
    print()


def dump_histogram(path: Path, analysis: str, label: str, xmin: float | None, xmax: float | None) -> None:
    aos = yoda.read(str(path))
    dump_histogram_from_aos(aos, str(path), analysis, label, xmin, xmax)


def build_campaign_objects(
    base_dir: Path,
    tag: str,
    setup: str,
    order: str,
    analysis: str,
    analysis_variant: str | None,
    scale_variation: str,
) -> Tuple[Dict[str, object], Dict[str, Path]]:
    hel_files = find_order_yoda_inputs(
        base_dir=base_dir,
        tag=tag,
        setup=setup,
        order=order,
        analysis_variant=analysis_variant,
        scale_variation=scale_variation,
    )
    objects = build_dis_polarized_objects(
        setup=setup,
        zero_path=str(hel_files["00"]) if "00" in hel_files else None,
        pp_path=str(hel_files["PP"]),
        pm_path=str(hel_files["PM"]),
        mp_path=str(hel_files["MP"]) if "MP" in hel_files else None,
        mm_path=str(hel_files["MM"]) if "MM" in hel_files else None,
        analysis=analysis,
    )
    return objects, hel_files


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("yodas", nargs="*", help="Input analyzed YODA or YODA.GZ files")
    parser.add_argument("--base-dir", type=Path, help="Base DISPOL directory containing campaigns/")
    parser.add_argument("--tag", help="Campaign tag, e.g. rivetfo47")
    parser.add_argument("--setup", default="ALL", help="DIS setup for campaign mode (default: ALL)")
    parser.add_argument(
        "--order",
        choices=("LO", "POSNLO", "NEGNLO", "NLO"),
        help="Order to resolve from the campaign manifest in campaign mode",
    )
    parser.add_argument(
        "--analysis-variant",
        default=None,
        help="Optional analysis variant override for campaign mode (default: manifest value)",
    )
    parser.add_argument(
        "--scale-variation",
        default="nominal",
        help="Scale-variation label for campaign mode (default: nominal)",
    )
    parser.add_argument("--analysis", default="MC_DIS_BREIT", help="Analysis prefix (default: MC_DIS_BREIT)")
    parser.add_argument("--label", required=True, help="Histogram label, e.g. Q2, DQ2, pT2")
    parser.add_argument("--xmin", type=float, help="Only print bins with xmax > xmin")
    parser.add_argument("--xmax", type=float, help="Only print bins with xmin < xmax")
    args = parser.parse_args()

    if args.yodas:
        for raw_path in args.yodas:
            dump_histogram(Path(raw_path), args.analysis, args.label, args.xmin, args.xmax)
        return 0

    if not args.base_dir or not args.tag or not args.order:
        parser.error("Provide either input yodas or campaign mode arguments: --base-dir, --tag, and --order.")

    objects, hel_files = build_campaign_objects(
        base_dir=args.base_dir,
        tag=args.tag,
        setup=args.setup,
        order=args.order,
        analysis=args.analysis,
        analysis_variant=args.analysis_variant,
        scale_variation=args.scale_variation,
    )
    source = f"campaign={args.tag} setup={args.setup} order={args.order}"
    print("# resolved-helicity-inputs")
    for hel, path in sorted(hel_files.items()):
        print(f"# {hel}={path}")
    print()
    dump_histogram_from_aos(objects, source, args.analysis, args.label, args.xmin, args.xmax)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
