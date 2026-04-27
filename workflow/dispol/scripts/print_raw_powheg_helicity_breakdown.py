#!/usr/bin/env python3.10
"""Print helicity-resolved raw-POWHEG selected contributions for one order/channel.

This is the follow-up diagnostic to ``print_raw_powheg_channel_pieces.py``:
it lets us inspect whether a channel-level polarized piece such as the low-Q2
BGF contribution to ``DQ2`` comes from a stable helicity pattern or from one
pathological helicity sector.

For a polarized histogram request like ``DQ2``, the raw sidecars only contain
the underlying selected base histogram ``Q2`` for each helicity run. This
helper therefore prints the base-helicity bin heights together with the signed
polarized combination terms and the resulting channel-level delta.

To make setup-to-setup subtraction robust, the helper always emits the same
full helicity schema:

  00, PP, PM, MP, MM, t_PP, t_PM, t_MP, t_MM, D, A

For ``GAMMA``, the photon-only symmetry relations are used to populate the
missing entries consistently in the printed table:

  MP -> PM
  MM -> PP

  ALL/Z/CC:
    D = (PP + MM - PM - MP) / 4

  GAMMA:
    D = (PP - PM) / 2

The denominator shown is the corresponding ``00`` selected histogram, so the
last column is the channel-level asymmetry ``D / 00`` for that order/channel.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Optional

from print_raw_powheg_channel_pieces import (
    combine_setup_measurements,
    format_float,
    load_helicity_measurements,
    resolve_hist_request,
)


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
        "--order",
        choices=("POSNLO", "NEGNLO"),
        required=True,
        help="Which real-emission order branch to inspect.",
    )
    parser.add_argument(
        "--channel",
        choices=("combined", "Compton", "BGF"),
        default="combined",
        help="Winner channel to inspect (default: combined).",
    )
    parser.add_argument(
        "--analysis",
        default="MC_DIS_BREIT",
        help="Analysis prefix in the raw YODA paths (default: MC_DIS_BREIT).",
    )
    parser.add_argument(
        "--hist",
        default="DQ2",
        help="Polarized histogram suffix to inspect (default: DQ2).",
    )
    parser.add_argument("--xmin", type=float, default=None, help="Only print bins with xmax > xmin.")
    parser.add_argument("--xmax", type=float, default=None, help="Only print bins with xmin < xmax.")
    parser.add_argument(
        "--precision",
        type=int,
        default=6,
        help="Digits after the decimal point for floating-point output (default: 6).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base_dir = args.base_dir.resolve()
    setup = args.setup.upper()
    hist_mode, base_label = resolve_hist_request(args.hist)
    if hist_mode == "all":
        raise ValueError(
            "Helicity breakdown is only meaningful for sigma or polarized D-histograms; "
            f"got {args.hist!r}."
        )

    channel: Optional[str]
    if args.channel == "combined":
        channel = None
    else:
        channel = args.channel

    measurements = load_helicity_measurements(
        base_dir=base_dir,
        tag=args.tag,
        setup=setup,
        order=args.order,
        analysis=args.analysis,
        base_label=base_label,
        channel=channel,
    )
    sigma, delta, _ratio = combine_setup_measurements(setup, measurements)

    print(f"Campaign: {base_dir / 'campaigns' / args.tag}")
    print(f"Setup: {setup}")
    print(f"Order: {args.order}")
    print(f"Channel: {args.channel}")
    print(f"Histogram request: /{args.analysis}/{args.hist}")
    print(f"Raw base histogram: /{args.analysis}/{base_label}")
    print("Columns:")
    print("  00/PP/PM/MP/MM are the selected raw-helicity base-histogram bin heights.")
    print("  t_* are the signed contributions entering the polarized delta combination.")
    print("  D is the resulting polarized channel piece for this order/channel.")
    print("  A = D / 00 uses the 00 selected channel piece as denominator.")
    print()

    print(
        "bin  range"
        "                             00          PP          PM          MP          MM"
        "        t_PP        t_PM        t_MP        t_MM           D           A"
    )

    pp = measurements["PP"]
    pm = measurements["PM"]
    zero = measurements["00"]
    mp = measurements.get("MP")
    mm = measurements.get("MM")

    for index in range(pp.nbins):
        xlo = pp.edges[index]
        xhi = pp.edges[index + 1]
        if args.xmin is not None and xhi <= args.xmin:
            continue
        if args.xmax is not None and xlo >= args.xmax:
            continue

        zero_val = zero.values[index]
        pp_val = pp.values[index]
        pm_val = pm.values[index]
        d_val = delta.values[index]
        a_val = float("nan") if zero_val == 0.0 else d_val / zero_val

        if setup == "GAMMA":
            # The raw GAMMA runs only exist for 00, PP, and PM. Populate the
            # missing entries with the photon-symmetry identifications so the
            # printed table has the same schema as ALL/Z/CC and can be safely
            # subtracted column-by-column in follow-up diagnostics.
            mp_val = pm_val
            mm_val = pp_val
        else:
            if mp is None or mm is None:
                raise ValueError(f"Setup {setup} requires MP and MM helicity inputs.")
            mp_val = mp.values[index]
            mm_val = mm.values[index]

        t_pp = 0.25 * pp_val
        t_pm = -0.25 * pm_val
        t_mp = -0.25 * mp_val
        t_mm = 0.25 * mm_val
        print(
            f"{index + 1:3d}  [{xlo:7.2f}, {xhi:7.2f})"
            f"  {format_float(zero_val, args.precision):>10}"
            f"  {format_float(pp_val, args.precision):>10}"
            f"  {format_float(pm_val, args.precision):>10}"
            f"  {format_float(mp_val, args.precision):>10}"
            f"  {format_float(mm_val, args.precision):>10}"
            f"  {format_float(t_pp, args.precision):>10}"
            f"  {format_float(t_pm, args.precision):>10}"
            f"  {format_float(t_mp, args.precision):>10}"
            f"  {format_float(t_mm, args.precision):>10}"
            f"  {format_float(d_val, args.precision):>10}"
            f"  {format_float(a_val, args.precision):>10}"
        )

    print()
    print(f"Integrated 00 sum = {format_float(math.fsum(zero.values), args.precision)}")
    print(f"Integrated D  sum = {format_float(math.fsum(delta.values), args.precision)}")
    sigma_sum = math.fsum(sigma.values)
    asym_sum = float("nan") if sigma_sum == 0.0 else math.fsum(delta.values) / sigma_sum
    print(f"Integrated A = {format_float(asym_sum, args.precision)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
