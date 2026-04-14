#!/usr/bin/env python3
"""Pointwise Born-level NC gamma/Z interference audit.

This script compares the neutral-current Born-level electroweak pieces used by
Herwig and POLDIS in two stages:

1. A raw gamma/Z interference audit of the flavor-resolved coupling-propagator
   coefficients before any PDF weighting.
2. A pointwise inclusive-DIS audit of the channel kernels on a fixed (Q2, y)
   grid, including the derived interference quantity
   sigma_LL,int = ALL - GAMMA - Z.

The raw stage is an exact direct comparison between the Herwig D/N coefficient
split and the POLDIS COMBPDFVEC/COMBPDFAX interference terms.

For the pointwise stage, the exact comparison is performed on sigma0 and
sigma_LL. POLDIS does not expose separate Born-level PP/PM/MP/MM inclusive
weights in the same direct way, so the printed helicity table reconstructs them
from (sigma0, sigma_LL) with vanishing single-spin terms. That reconstruction
is exact for the displayed sigma0 and sigma_LL and is meant only as a compact
view of the channel bookkeeping.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple


DEFAULT_SIN2_THETA_W = 0.221639970483740179
DEFAULT_ALPHA_EM = 0.00729927
DEFAULT_MZ = 91.1876
DEFAULT_GAMMA_Z = 2.4952
DEFAULT_MW = 80.450
DEFAULT_GAMMA_W = 2.085
DEFAULT_LEPTON_CHARGE = -1.0
DEFAULT_Q2_VALUES = (100.0, 500.0, 1500.0)
DEFAULT_Y_VALUES = (0.25, 0.40, 0.55)
DEFAULT_FLAVORS = ("u", "d", "ubar", "dbar")


@dataclass(frozen=True)
class EWInputs:
    sin2_theta_w: float
    alpha_em: float
    mz: float
    gamma_z: float
    mw: float
    gamma_w: float
    lepton_charge: float


@dataclass(frozen=True)
class Flavor:
    name: str
    herwig_pdg_id: int
    poldis_index: int
    species_charge: float
    is_up_type: bool


@dataclass
class NCCoefficients:
    D0: float = 0.0
    Dl: float = 0.0
    Dq: float = 0.0
    Dlq: float = 0.0
    N0: float = 0.0
    Nl: float = 0.0
    Nq: float = 0.0
    Nlq: float = 0.0

    def __add__(self, other: "NCCoefficients") -> "NCCoefficients":
        return NCCoefficients(
            D0=self.D0 + other.D0,
            Dl=self.Dl + other.Dl,
            Dq=self.Dq + other.Dq,
            Dlq=self.Dlq + other.Dlq,
            N0=self.N0 + other.N0,
            Nl=self.Nl + other.Nl,
            Nq=self.Nq + other.Nq,
            Nlq=self.Nlq + other.Nlq,
        )


@dataclass(frozen=True)
class POLDISChannelPieces:
    vec_gg: float
    vec_gz: float
    vec_zz: float
    ax_gz: float
    ax_zz: float


FLAVORS: Dict[str, Flavor] = {
    "u": Flavor("u", 2, 2, 2.0 / 3.0, True),
    "d": Flavor("d", 1, 1, -1.0 / 3.0, False),
    "ubar": Flavor("ubar", -2, -2, 2.0 / 3.0, True),
    "dbar": Flavor("dbar", -1, -1, -1.0 / 3.0, False),
}

HELICITY_COMBINATIONS: Tuple[Tuple[str, int, int], ...] = (
    ("PP", +1, +1),
    ("PM", +1, -1),
    ("MP", -1, +1),
    ("MM", -1, -1),
)

CHANNEL_ORDER = ("GAMMA", "Z", "ALL", "INT")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare Born-level NC gamma/Z interference coefficients "
        "between Herwig and POLDIS."
    )
    parser.add_argument(
        "--sin2-theta-w",
        type=float,
        default=DEFAULT_SIN2_THETA_W,
        help="Matched sin^2(theta_W) used by the current Herwig draft tables.",
    )
    parser.add_argument(
        "--alpha-em",
        type=float,
        default=DEFAULT_ALPHA_EM,
        help="Fixed alpha_EM input used in both setups.",
    )
    parser.add_argument("--mz", type=float, default=DEFAULT_MZ)
    parser.add_argument("--gamma-z", type=float, default=DEFAULT_GAMMA_Z)
    parser.add_argument("--mw", type=float, default=DEFAULT_MW)
    parser.add_argument("--gamma-w", type=float, default=DEFAULT_GAMMA_W)
    parser.add_argument(
        "--lepton-charge",
        type=float,
        default=DEFAULT_LEPTON_CHARGE,
        help="Incoming lepton charge: -1 for e-, +1 for e+.",
    )
    parser.add_argument(
        "--q2",
        nargs="+",
        type=float,
        default=list(DEFAULT_Q2_VALUES),
        help="Q^2 values in GeV^2 for the pointwise audit.",
    )
    parser.add_argument(
        "--y",
        nargs="+",
        type=float,
        default=list(DEFAULT_Y_VALUES),
        help="DIS y values for the pointwise audit.",
    )
    parser.add_argument(
        "--flavors",
        nargs="+",
        choices=sorted(FLAVORS),
        default=list(DEFAULT_FLAVORS),
        help="Flavor list for the audit grid.",
    )
    parser.add_argument(
        "--show-herwig-exact-helicities",
        action="store_true",
        help="Also print the exact Herwig helicity weights, including single-spin pieces.",
    )
    return parser.parse_args()


def build_ew_inputs(args: argparse.Namespace) -> EWInputs:
    return EWInputs(
        sin2_theta_w=args.sin2_theta_w,
        alpha_em=args.alpha_em,
        mz=args.mz,
        gamma_z=args.gamma_z,
        mw=args.mw,
        gamma_w=args.gamma_w,
        lepton_charge=args.lepton_charge,
    )


def theta_w(inputs: EWInputs) -> float:
    return math.asin(math.sqrt(inputs.sin2_theta_w))


def poldis_lepton_couplings(inputs: EWInputs) -> Tuple[float, float]:
    theta = theta_w(inputs)
    cvze = (-0.5 / math.sin(2.0 * theta) + math.tan(theta)) * (-inputs.lepton_charge)
    caze = -0.5 / math.sin(2.0 * theta)
    return cvze, caze


def poldis_quark_couplings(flavor: Flavor, inputs: EWInputs) -> Tuple[float, float, float]:
    theta = theta_w(inputs)
    if flavor.is_up_type:
        cvz = 0.5 / math.sin(2.0 * theta) - (2.0 / 3.0) * math.tan(theta)
        caz = 0.5 / math.sin(2.0 * theta)
    else:
        cvz = -0.5 / math.sin(2.0 * theta) + (1.0 / 3.0) * math.tan(theta)
        caz = -0.5 / math.sin(2.0 * theta)

    if flavor.poldis_index < 0:
        cvz = -cvz

    eq = flavor.species_charge if flavor.poldis_index > 0 else -flavor.species_charge
    return eq, cvz, caz


def herwig_quark_couplings(flavor: Flavor, inputs: EWInputs) -> Tuple[float, float, float]:
    s2 = inputs.sin2_theta_w
    if flavor.is_up_type:
        aq = 1.0
        vq = aq - 4.0 * flavor.species_charge * s2
    else:
        aq = -1.0
        vq = aq - 4.0 * flavor.species_charge * s2
    return flavor.species_charge, vq / 4.0, aq / 4.0


def herwig_lepton_couplings(inputs: EWInputs) -> Tuple[float, float, float]:
    ql = inputs.lepton_charge
    al = -1.0 if inputs.lepton_charge < 0.0 else 1.0
    vl = al - 4.0 * ql * inputs.sin2_theta_w
    return ql, vl / 4.0, al / 4.0


def herwig_eta_l(inputs: EWInputs) -> float:
    return 1.0 if inputs.lepton_charge < 0.0 else -1.0


def herwig_eta_q(flavor: Flavor) -> float:
    return -1.0 if flavor.herwig_pdg_id < 0 else 1.0


def propagator_factors(q2: float, inputs: EWInputs) -> Tuple[float, float]:
    den = (q2 + inputs.mz * inputs.mz) ** 2 + (inputs.mz * inputs.gamma_z) ** 2
    propz = q2 * q2 / den
    intgz = inputs.lepton_charge * (q2 + inputs.mz * inputs.mz) / q2
    return propz, intgz


def herwig_nc_components(flavor: Flavor, q2: float, inputs: EWInputs) -> Dict[str, NCCoefficients]:
    ql, cvl, cal = herwig_lepton_couplings(inputs)
    qq, cvq, caq = herwig_quark_couplings(flavor, inputs)
    eta_l = herwig_eta_l(inputs)
    eta_q = herwig_eta_q(flavor)

    sw = math.sqrt(inputs.sin2_theta_w)
    cw = math.sqrt(1.0 - inputs.sin2_theta_w)
    den = (q2 + inputs.mz * inputs.mz) ** 2 + (inputs.mz * inputs.gamma_z) ** 2
    k = 1.0 / (sw * cw) ** 2
    z_int = k * q2 * (q2 + inputs.mz * inputs.mz) / den
    z_sq = (k * k) * q2 * q2 / den

    gg = NCCoefficients(
        D0=ql * ql * qq * qq,
        Nlq=2.0 * ql * ql * qq * qq,
    )
    gz = NCCoefficients(
        D0=2.0 * z_int * ql * qq * cvl * cvq,
        Dl=eta_l * -2.0 * cal * (z_int * cvq * ql * qq),
        Dq=eta_q * -2.0 * caq * (z_int * cvl * ql * qq),
        Dlq=eta_l * eta_q * 2.0 * cal * caq * (z_int * ql * qq),
        N0=eta_l * eta_q * 4.0 * cal * caq * (z_int * ql * qq),
        Nl=eta_q * -4.0 * caq * (z_int * cvl * ql * qq),
        Nq=eta_l * -4.0 * cal * (z_int * cvq * ql * qq),
        Nlq=4.0 * z_int * cvl * cvq * ql * qq,
    )
    zz = NCCoefficients(
        D0=z_sq * (cvl * cvl + cal * cal) * (cvq * cvq + caq * caq),
        Dl=eta_l * -2.0 * cal * (z_sq * cvl * (caq * caq + cvq * cvq)),
        Dq=eta_q * -2.0 * caq * (z_sq * cvq * (cal * cal + cvl * cvl)),
        Dlq=eta_l * eta_q * 4.0 * cal * caq * z_sq * cvl * cvq,
        N0=eta_l * eta_q * 8.0 * cal * caq * z_sq * cvl * cvq,
        Nl=eta_q * -4.0 * caq * (z_sq * cvq * (cal * cal + cvl * cvl)),
        Nq=eta_l * -4.0 * cal * (z_sq * cvl * (caq * caq + cvq * cvq)),
        Nlq=2.0 * z_sq * (caq * caq + cvq * cvq) * (cal * cal + cvl * cvl),
    )
    return {"gg": gg, "gz": gz, "zz": zz}


def poldis_channel_pieces(flavor: Flavor, q2: float, inputs: EWInputs) -> POLDISChannelPieces:
    eq, cvz, caz = poldis_quark_couplings(flavor, inputs)
    cvze, caze = poldis_lepton_couplings(inputs)
    propz, intgz = propagator_factors(q2, inputs)

    return POLDISChannelPieces(
        vec_gg=eq * eq,
        vec_gz=2.0 * eq * propz * intgz * (cvz * cvze),
        vec_zz=propz * (cvz * cvz + caz * caz) * (cvze * cvze + caze * caze),
        ax_gz=2.0 * eq * propz * intgz * (caz * caze),
        ax_zz=propz * (4.0 * cvz * caz * cvze * caze),
    )


def y_plus(y: float) -> float:
    return 1.0 + (1.0 - y) ** 2


def y_minus(y: float) -> float:
    return 1.0 - (1.0 - y) ** 2


def ell_from_y(y: float) -> float:
    return 2.0 / y - 1.0


def herwig_poldis_kernel_factor(ell: float) -> float:
    return 2.0 / (1.0 + ell) ** 2


def sigma_weight(coeff: NCCoefficients, ell: float, pl: int, pq: int) -> float:
    return (1.0 + ell * ell) * (
        coeff.D0 + pl * coeff.Dl + pq * coeff.Dq + pl * pq * coeff.Dlq
    ) + ell * (
        coeff.N0 + pl * coeff.Nl + pq * coeff.Nq + pl * pq * coeff.Nlq
    )


def herwig_exact_channel_weights(
    coeff: NCCoefficients, ell: float
) -> Dict[str, float]:
    kernel = herwig_poldis_kernel_factor(ell)
    return {
        hel: kernel * sigma_weight(coeff, ell, pl, pq)
        for hel, pl, pq in HELICITY_COMBINATIONS
    }


def helicity_combo_sigma0(weights: Dict[str, float]) -> float:
    return 0.25 * (weights["PP"] + weights["PM"] + weights["MP"] + weights["MM"])


def helicity_combo_sigma_ll(weights: Dict[str, float]) -> float:
    return 0.25 * (weights["PP"] + weights["MM"] - weights["PM"] - weights["MP"])


def poldis_channel_observables(
    pieces: POLDISChannelPieces, channel: str, y: float
) -> Tuple[float, float]:
    yp = y_plus(y)
    ym = y_minus(y)

    if channel == "GAMMA":
        vec = pieces.vec_gg
        ax = 0.0
    elif channel == "Z":
        vec = pieces.vec_zz
        ax = pieces.ax_zz
    elif channel == "ALL":
        vec = pieces.vec_gg + pieces.vec_gz + pieces.vec_zz
        ax = pieces.ax_gz + pieces.ax_zz
    else:
        raise ValueError(f"Unknown channel {channel}")

    sigma0 = yp * vec + ym * ax
    sigma_ll = ym * vec + yp * ax
    return sigma0, sigma_ll


def reconstruct_helicities(sigma0: float, sigma_ll: float) -> Dict[str, float]:
    return {
        "PP": sigma0 + sigma_ll,
        "PM": sigma0 - sigma_ll,
        "MP": sigma0 - sigma_ll,
        "MM": sigma0 + sigma_ll,
    }


def safe_ratio(num: float, den: float) -> str:
    if abs(den) <= 1e-15:
        return "n/a"
    return f"{num / den:+.12f}"


def format_float(value: float) -> str:
    return f"{value:+.12e}"


def build_table(headers: Sequence[str], rows: Iterable[Sequence[str]]) -> str:
    rows = [list(row) for row in rows]
    widths = [len(h) for h in headers]
    for row in rows:
        for idx, cell in enumerate(row):
            widths[idx] = max(widths[idx], len(cell))

    def fmt_row(row: Sequence[str]) -> str:
        return "  ".join(cell.ljust(widths[idx]) for idx, cell in enumerate(row))

    sep = "  ".join("-" * width for width in widths)
    lines = [fmt_row(headers), sep]
    lines.extend(fmt_row(row) for row in rows)
    return "\n".join(lines)


def stage1_rows(flavor: Flavor, q2: float, inputs: EWInputs) -> Tuple[List[List[str]], List[float]]:
    herwig = herwig_nc_components(flavor, q2, inputs)["gz"]
    poldis = poldis_channel_pieces(flavor, q2, inputs)

    comparisons = (
        ("vector-even", herwig.D0, poldis.vec_gz),
        ("axial-even", 0.5 * herwig.N0, poldis.ax_gz),
        ("vector-spin", 0.5 * herwig.Nlq, poldis.vec_gz),
        ("axial-spin", herwig.Dlq, poldis.ax_gz),
    )

    rows: List[List[str]] = []
    deltas: List[float] = []
    for piece, hw, pd in comparisons:
        rows.append(
            [
                flavor.name,
                piece,
                format_float(hw),
                format_float(pd),
                safe_ratio(hw, pd),
                format_float(hw - pd),
            ]
        )
        if abs(pd) > 1e-15:
            deltas.append(abs(hw / pd - 1.0))
    return rows, deltas


def channel_coefficients(coeffs: Dict[str, NCCoefficients]) -> Dict[str, NCCoefficients]:
    gamma = coeffs["gg"]
    z_only = coeffs["zz"]
    all_nc = coeffs["gg"] + coeffs["gz"] + coeffs["zz"]
    interference = coeffs["gz"]
    return {
        "GAMMA": gamma,
        "Z": z_only,
        "ALL": all_nc,
        "INT": interference,
    }


def poldis_channel_summaries(
    pieces: POLDISChannelPieces, y: float
) -> Dict[str, Tuple[float, float]]:
    gamma = poldis_channel_observables(pieces, "GAMMA", y)
    z_only = poldis_channel_observables(pieces, "Z", y)
    all_nc = poldis_channel_observables(pieces, "ALL", y)
    interference = (
        all_nc[0] - gamma[0] - z_only[0],
        all_nc[1] - gamma[1] - z_only[1],
    )
    return {
        "GAMMA": gamma,
        "Z": z_only,
        "ALL": all_nc,
        "INT": interference,
    }


def stage2_point(
    flavor: Flavor,
    q2: float,
    y: float,
    inputs: EWInputs,
    show_exact_herwig_helicities: bool,
) -> Tuple[str, List[float]]:
    coeffs = channel_coefficients(herwig_nc_components(flavor, q2, inputs))
    pieces = poldis_channel_pieces(flavor, q2, inputs)
    poldis = poldis_channel_summaries(pieces, y)
    ell = ell_from_y(y)

    summary_rows: List[List[str]] = []
    helicity_rows: List[List[str]] = []
    exact_rows: List[List[str]] = []
    deltas: List[float] = []

    for channel in CHANNEL_ORDER:
        herwig_weights = herwig_exact_channel_weights(coeffs[channel], ell)
        herwig_sigma0 = helicity_combo_sigma0(herwig_weights)
        herwig_sigma_ll = helicity_combo_sigma_ll(herwig_weights)
        poldis_sigma0, poldis_sigma_ll = poldis[channel]

        summary_rows.append(
            [
                channel,
                format_float(herwig_sigma0),
                format_float(poldis_sigma0),
                safe_ratio(herwig_sigma0, poldis_sigma0),
                format_float(herwig_sigma_ll),
                format_float(poldis_sigma_ll),
                safe_ratio(herwig_sigma_ll, poldis_sigma_ll),
            ]
        )

        if abs(poldis_sigma0) > 1e-15:
            deltas.append(abs(herwig_sigma0 / poldis_sigma0 - 1.0))
        if abs(poldis_sigma_ll) > 1e-15:
            deltas.append(abs(herwig_sigma_ll / poldis_sigma_ll - 1.0))

        herwig_reco = reconstruct_helicities(herwig_sigma0, herwig_sigma_ll)
        poldis_reco = reconstruct_helicities(poldis_sigma0, poldis_sigma_ll)
        for hel in ("PP", "PM", "MP", "MM"):
            helicity_rows.append(
                [
                    channel,
                    hel,
                    format_float(herwig_reco[hel]),
                    format_float(poldis_reco[hel]),
                    safe_ratio(herwig_reco[hel], poldis_reco[hel]),
                ]
            )

        if show_exact_herwig_helicities:
            for hel in ("PP", "PM", "MP", "MM"):
                exact_rows.append([channel, hel, format_float(herwig_weights[hel])])

    lines = [
        f"Point: flavor={flavor.name}  Q2={q2:.3f} GeV^2  y={y:.3f}  ell={ell:.6f}",
        build_table(
            (
                "channel",
                "Herwig sigma0",
                "POLDIS sigma0",
                "ratio",
                "Herwig sigma_LL",
                "POLDIS sigma_LL",
                "ratio",
            ),
            summary_rows,
        ),
        "",
        "Reconstructed helicities from (sigma0, sigma_LL):",
        build_table(
            ("channel", "helicity", "Herwig(reco)", "POLDIS(reco)", "ratio"),
            helicity_rows,
        ),
    ]

    if show_exact_herwig_helicities:
        lines.extend(
            [
                "",
                "Exact Herwig helicities (includes single-spin pieces):",
                build_table(("channel", "helicity", "Herwig exact"), exact_rows),
            ]
        )

    return "\n".join(lines), deltas


def print_inputs(inputs: EWInputs, q2_values: Sequence[float], y_values: Sequence[float], flavors: Sequence[str]) -> None:
    print("NC Born gamma/Z interference audit")
    print("---------------------------------")
    print(
        "Matched inputs: "
        f"sin^2(theta_W)={inputs.sin2_theta_w:.18f}, "
        f"alpha_EM={inputs.alpha_em:.8f}, "
        f"mZ={inputs.mz:.4f}, GammaZ={inputs.gamma_z:.4f}, "
        f"mW={inputs.mw:.4f}, GammaW={inputs.gamma_w:.4f}, "
        f"LEPCHARGE={inputs.lepton_charge:+.1f}"
    )
    print(
        "Grid: "
        f"flavors={','.join(flavors)}  "
        f"Q2={','.join(f'{q:.3f}' for q in q2_values)}  "
        f"y={','.join(f'{val:.3f}' for val in y_values)}"
    )
    print()


def main() -> None:
    args = parse_args()
    inputs = build_ew_inputs(args)
    q2_values = tuple(args.q2)
    y_values = tuple(args.y)
    flavors = tuple(args.flavors)

    print_inputs(inputs, q2_values, y_values, flavors)

    stage1_deviations: List[float] = []
    print("Stage 1: raw gamma/Z interference coefficients")
    print("----------------------------------------------")
    for q2 in q2_values:
        rows: List[List[str]] = []
        for flavor_name in flavors:
            flavor_rows, deltas = stage1_rows(FLAVORS[flavor_name], q2, inputs)
            rows.extend(flavor_rows)
            stage1_deviations.extend(deltas)
        print(f"Q2 = {q2:.3f} GeV^2")
        print(
            build_table(
                ("flavor", "piece", "Herwig", "POLDIS", "ratio", "delta"),
                rows,
            )
        )
        print()

    stage2_deviations: List[float] = []
    print("Stage 2: pointwise channel kernels")
    print("----------------------------------")
    for flavor_name in flavors:
        flavor = FLAVORS[flavor_name]
        for q2 in q2_values:
            for y in y_values:
                report, deltas = stage2_point(
                    flavor,
                    q2,
                    y,
                    inputs,
                    show_exact_herwig_helicities=args.show_herwig_exact_helicities,
                )
                stage2_deviations.extend(deltas)
                print(report)
                print()

    max_stage1 = max(stage1_deviations) if stage1_deviations else 0.0
    max_stage2 = max(stage2_deviations) if stage2_deviations else 0.0
    print("Summary")
    print("-------")
    print(f"Max Stage-1 relative deviation: {max_stage1:.3e}")
    print(f"Max Stage-2 relative deviation: {max_stage2:.3e}")


if __name__ == "__main__":
    main()
