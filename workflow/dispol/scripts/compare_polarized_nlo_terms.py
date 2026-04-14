#!/usr/bin/env python3.10
"""
Pointwise polarized NLO audit for the common Herwig/POLDIS DIS kernels.

This script compares the live Herwig polarized NLO building blocks against a
flavor-resolved single-flavor analogue of the polarized POLDIS G2/G4/GL
inclusive coefficient-function decomposition.

The mapping is intentionally pragmatic:
  - Herwig `virtual correction + quark collinear` is compared to the POLDIS
    quark singular-like (+/delta) contribution.
  - Herwig `gluon collinear` is compared to zero on the POLDIS side, because
    the inclusive POLDIS gluon channel has no separate singular piece.
  - Herwig `realq` is compared to the POLDIS quark regular piece.
  - Herwig `realg` is compared to the POLDIS gluon regular piece.

So this is not a one-to-one local-subtraction audit. It is a pointwise audit of
the closest analytically meaningful pieces that can be compared before Monte
Carlo integration.
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

from compare_pdf_alphas import (
    DEFAULT_INPUT_ALPHA_S,
    DEFAULT_INPUT_SCALE,
    default_matchbox_quark_masses2,
    matchbox_running_alpha,
    solve_matchbox_lambda_squared,
)


DEFAULT_LHAPDF_DATA_PATH = "/opt/homebrew/Cellar/lhapdf/6.5.4/share/LHAPDF"
DEFAULT_SUM_PDF = "PDF4LHC15_nnlo_100_pdfas"
DEFAULT_DIFF_PDF = "BDSSV24-NNLO"
DEFAULT_CHANNEL = "ALL"
DEFAULT_Q2_VALUES = [49.0, 100.0, 500.0, 1500.0, 2500.0]
DEFAULT_Y_VALUES = [0.2, 0.4, 0.6]
DEFAULT_XP_FRACTIONS = [0.25, 0.5, 0.75]
BEAM_E_A_GEV = 18.0
BEAM_E_B_GEV = 275.0
DEFAULT_S = 4.0 * BEAM_E_A_GEV * BEAM_E_B_GEV


@dataclass(frozen=True)
class Flavor:
    label: str
    lhapdf_id: int
    species_charge: float
    is_up_type: bool


@dataclass(frozen=True)
class EWInputs:
    sin2_theta_w: float
    mz: float
    gamma_z: float
    lepton_charge: float


@dataclass(frozen=True)
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
class ChannelPieces:
    vec: float
    ax: float


@dataclass(frozen=True)
class HerwigHelicityTerms:
    born_abs: float
    virt_corr: float
    quark_collinear: float
    collq_k1: float
    collq_k2: float
    collq_odd: float
    gluon_collinear: float
    qcdc_mapped: float
    qcdc_even: float
    qcdc_aq: float
    bgf_mapped: float
    Pq: float
    Pq_m: float
    Pg_m: float
    dq_ratio: float
    dg_ratio: float
    deltaq_over_lo: float
    deltag_over_lo: float
    q_odd_response: float
    g_odd_response: float
    a_born: float
    a_q_mapped: float
    a_g_r2: float
    a_g_r3: float
    bgf_split_unpol: float
    bgf_split_odd: float
    bgf_legacy_unpol: float
    bgf_legacy_spin: float
    bgf_legacy_total: float


@dataclass(frozen=True)
class NCAuditPayload:
    d_even: float
    d_spin: float
    n_even: float
    n_spin: float
    q_unpolarized: float
    q_polarized: float
    g_unpolarized: float
    g_polarized: float
    q_odd_response: float
    g_odd_response: float
    born_factor: float
    real_denominator_factor: float
    mapped_denominator_ratio: float


@dataclass(frozen=True)
class AggregatedTerms:
    born_abs: float
    virtual_abs: float
    collq_abs: float
    collq_k1_abs: float
    collq_k2_abs: float
    collq_odd_abs: float
    quark_collinear: float
    gluon_collinear: float
    qcdc_mapped: float
    bgf_mapped: float
    quark_collinear_abs: float
    gluon_collinear_abs: float
    qcdc_mapped_abs: float
    qcdc_even_abs: float
    qcdc_aq_abs: float
    bgf_mapped_abs: float
    nlo_sum_abs: float
    inputs: Dict[str, float]


FLAVORS: List[Flavor] = [
    Flavor("u", 2, 2.0 / 3.0, True),
    Flavor("d", 1, -1.0 / 3.0, False),
    Flavor("ubar", -2, 2.0 / 3.0, True),
    Flavor("dbar", -1, -1.0 / 3.0, False),
]

HELICITY_COMBINATIONS: List[Tuple[str, int, int]] = [
    ("PP", +1, +1),
    ("PM", +1, -1),
    ("MP", -1, +1),
    ("MM", -1, -1),
]


def parse_float_list(text: str) -> List[float]:
    out: List[float] = []
    for chunk in text.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        value = float(chunk)
        if value <= 0.0:
            raise argparse.ArgumentTypeError("Values must be positive.")
        out.append(value)
    if not out:
        raise argparse.ArgumentTypeError("Provide at least one comma-separated value.")
    return out


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--channel",
        choices=("GAMMA", "Z", "ALL"),
        default=DEFAULT_CHANNEL,
        help="NC channel to audit.",
    )
    parser.add_argument(
        "--sum-pdf",
        default=DEFAULT_SUM_PDF,
        help="Unpolarized PDF used for the Herwig sum PDF.",
    )
    parser.add_argument(
        "--diff-pdf",
        default=DEFAULT_DIFF_PDF,
        help="Polarized PDF used for the Herwig diff PDF and the POLDIS polarized kernels.",
    )
    parser.add_argument(
        "--q2-values",
        type=parse_float_list,
        default=list(DEFAULT_Q2_VALUES),
        help="Comma-separated Q^2 values in GeV^2.",
    )
    parser.add_argument(
        "--y-values",
        type=parse_float_list,
        default=list(DEFAULT_Y_VALUES),
        help="Comma-separated inelasticity values.",
    )
    parser.add_argument(
        "--xp-fractions",
        type=parse_float_list,
        default=list(DEFAULT_XP_FRACTIONS),
        help="Comma-separated fractions used to place xp between x_B and 1.",
    )
    parser.add_argument(
        "--sin2-theta-w",
        type=float,
        default=0.221639970483740179,
        help="sin^2(theta_W) used in the matched old-table setup.",
    )
    parser.add_argument("--mz", type=float, default=91.1876, help="Z mass in GeV.")
    parser.add_argument("--gamma-z", type=float, default=2.4952, help="Z width in GeV.")
    parser.add_argument(
        "--input-alpha-s",
        type=float,
        default=DEFAULT_INPUT_ALPHA_S,
        help="Herwig Matchbox NLOAlphaS input value at the input scale.",
    )
    parser.add_argument(
        "--input-scale",
        type=float,
        default=DEFAULT_INPUT_SCALE,
        help="Herwig Matchbox NLOAlphaS input scale in GeV.",
    )
    parser.add_argument(
        "--s-hadronic",
        type=float,
        default=DEFAULT_S,
        help="ep hadronic s in GeV^2. Defaults to 4 * 18 GeV * 275 GeV.",
    )
    parser.add_argument(
        "--flavor",
        action="append",
        choices=[flavor.label for flavor in FLAVORS],
        help="Optional repeated flavor filter.",
    )
    parser.add_argument(
        "--bgf-debug",
        action="store_true",
        help="Print a compact Herwig-only BGF diagnostic block for the PP helicity.",
    )
    parser.add_argument(
        "--absolute-debug",
        action="store_true",
        help="Print absolute sigma_LL contribution tables alongside the normalized ratios.",
    )
    parser.add_argument(
        "--c4-diagnostic",
        action="store_true",
        help="Add the exact POLDIS-basis quark regular (c4-c2) piece as a diagnostic-only shift on the Herwig qcdc term.",
    )
    return parser.parse_args(argv)


def import_lhapdf():
    if "LHAPDF_DATA_PATH" not in os.environ:
        os.environ["LHAPDF_DATA_PATH"] = DEFAULT_LHAPDF_DATA_PATH
    try:
        import lhapdf
    except ImportError as exc:
        raise SystemExit(
            "Could not import lhapdf in python3.10. Set LHAPDF_DATA_PATH if needed "
            "and use a Python environment with LHAPDF bindings installed."
        ) from exc
    return lhapdf


def resolve_flavors(requested: Sequence[str] | None) -> List[Flavor]:
    if not requested:
        return list(FLAVORS)
    wanted = set(requested)
    return [flavor for flavor in FLAVORS if flavor.label in wanted]


def conjugate_flavor(flavor: Flavor) -> Flavor:
    for candidate in FLAVORS:
        if candidate.lhapdf_id == -flavor.lhapdf_id:
            return candidate
    raise ValueError(f"Could not find charge-conjugate flavor for {flavor.label}")


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
    if flavor.lhapdf_id < 0:
        cvz = -cvz
    eq = flavor.species_charge if flavor.lhapdf_id > 0 else -flavor.species_charge
    return eq, cvz, caz


def herwig_quark_couplings(flavor: Flavor, inputs: EWInputs) -> Tuple[float, float, float]:
    s2 = inputs.sin2_theta_w
    aq = 1.0 if flavor.is_up_type else -1.0
    vq = aq - 4.0 * flavor.species_charge * s2
    return flavor.species_charge, vq / 4.0, aq / 4.0


def herwig_lepton_couplings(inputs: EWInputs) -> Tuple[float, float, float]:
    ql = inputs.lepton_charge
    al = -1.0 if ql < 0.0 else 1.0
    vl = al - 4.0 * ql * inputs.sin2_theta_w
    return ql, vl / 4.0, al / 4.0


def herwig_eta_l(inputs: EWInputs) -> float:
    return 1.0 if inputs.lepton_charge < 0.0 else -1.0


def herwig_eta_q(flavor: Flavor) -> float:
    return -1.0 if flavor.lhapdf_id < 0 else 1.0


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
    return {"GAMMA": gg, "GZ": gz, "Z": zz}


def herwig_channel_coefficients(flavor: Flavor, q2: float, inputs: EWInputs, channel: str) -> NCCoefficients:
    pieces = herwig_nc_components(flavor, q2, inputs)
    if channel == "GAMMA":
        return pieces["GAMMA"]
    if channel == "Z":
        return pieces["Z"]
    if channel == "ALL":
        return pieces["GAMMA"] + pieces["GZ"] + pieces["Z"]
    raise ValueError(f"Unsupported channel {channel}")


def poldis_channel_pieces(flavor: Flavor, q2: float, inputs: EWInputs, channel: str) -> ChannelPieces:
    eq, cvz, caz = poldis_quark_couplings(flavor, inputs)
    cvze, caze = poldis_lepton_couplings(inputs)
    propz, intgz = propagator_factors(q2, inputs)
    gamma_vec = eq * eq
    gz_vec = 2.0 * eq * propz * intgz * (cvz * cvze)
    z_vec = propz * (cvz * cvz + caz * caz) * (cvze * cvze + caze * caze)
    gz_ax = 2.0 * eq * propz * intgz * (caz * caze)
    z_ax = propz * (4.0 * cvz * caz * cvze * caze)
    if channel == "GAMMA":
        return ChannelPieces(vec=gamma_vec, ax=0.0)
    if channel == "Z":
        return ChannelPieces(vec=z_vec, ax=z_ax)
    if channel == "ALL":
        return ChannelPieces(vec=gamma_vec + gz_vec + z_vec, ax=gz_ax + z_ax)
    raise ValueError(f"Unsupported channel {channel}")


def ell_from_y(y: float) -> float:
    return 2.0 / y - 1.0


def sigma_weight(coeff: NCCoefficients, ell: float, pl: int, pq: float) -> float:
    return (1.0 + ell * ell) * (
        coeff.D0 + pl * coeff.Dl + pq * coeff.Dq + pl * pq * coeff.Dlq
    ) + ell * (
        coeff.N0 + pl * coeff.Nl + pq * coeff.Nq + pl * pq * coeff.Nlq
    )


def herwig_kernel_prefactor(ell: float) -> float:
    return 2.0 / (1.0 + ell) ** 2


def clamp_pol(value: float) -> float:
    return max(-1.0, min(1.0, value))


def xfx_over_x(member, pid: int, x: float, q2: float) -> float:
    return member.xfxQ2(pid, x, q2) / x


def real_emission_denominator(coeff: NCCoefficients, pl: int, pq: float) -> float:
    return coeff.D0 + pl * coeff.Dl + pq * coeff.Dq + pl * pq * coeff.Dlq


def real_emission_numerator(coeff: NCCoefficients, pl: int, pq: float) -> float:
    return coeff.N0 + pl * coeff.Nl + pq * coeff.Nq + pl * pq * coeff.Nlq


def nc_audit_payload(
    coeff: NCCoefficients, ell: float, pl: int, pq_born: float, pq_mapped: float
) -> NCAuditPayload:
    d_even = coeff.D0 + pl * coeff.Dl
    d_spin = pq_born * (coeff.Dq + pl * coeff.Dlq)
    n_even = coeff.N0 + pl * coeff.Nl
    n_spin = pq_born * (coeff.Nq + pl * coeff.Nlq)
    born_factor = (1.0 + ell * ell) * (d_even + d_spin) + ell * (n_even + n_spin)
    real_denominator_factor = real_emission_denominator(coeff, pl, pq_born)
    mapped_denominator = real_emission_denominator(coeff, pl, pq_mapped)
    mapped_denominator_ratio = (
        mapped_denominator / real_denominator_factor
        if abs(real_denominator_factor) > 1.0e-30
        else 1.0
    )
    if abs(born_factor) <= 1.0e-30:
        return NCAuditPayload(
            d_even=d_even,
            d_spin=d_spin,
            n_even=n_even,
            n_spin=n_spin,
            q_unpolarized=1.0,
            q_polarized=0.0,
            g_unpolarized=1.0,
            g_polarized=0.0,
            q_odd_response=0.0,
            g_odd_response=0.0,
            born_factor=born_factor,
            real_denominator_factor=real_denominator_factor,
            mapped_denominator_ratio=mapped_denominator_ratio,
        )

    q_odd_response = (
        (1.0 + ell * ell) * (coeff.Dq + pl * coeff.Dlq)
        + ell * (coeff.Nq + pl * coeff.Nlq)
    ) / born_factor
    g_odd_response = (ell * (coeff.Nq + pl * coeff.Nlq)) / born_factor
    q_unpolarized = ((1.0 + ell * ell) * d_even + ell * n_even) / born_factor
    g_unpolarized = ((1.0 + ell * ell) * d_even) / born_factor
    return NCAuditPayload(
        d_even=d_even,
        d_spin=d_spin,
        n_even=n_even,
        n_spin=n_spin,
        q_unpolarized=q_unpolarized,
        q_polarized=pq_born * q_odd_response,
        g_unpolarized=g_unpolarized,
        g_polarized=pq_born * g_odd_response,
        q_odd_response=q_odd_response,
        g_odd_response=g_odd_response,
        born_factor=born_factor,
        real_denominator_factor=real_denominator_factor,
        mapped_denominator_ratio=mapped_denominator_ratio,
    )


def evaluate_herwig_helicity_terms(
    flavor: Flavor,
    q2: float,
    y: float,
    xp: float,
    s_hadronic: float,
    channel: str,
    inputs: EWInputs,
    sum_pdf,
    diff_pdf,
    herwig_alpha_s: float,
) -> Dict[str, HerwigHelicityTerms]:
    ell = ell_from_y(y)
    x_b = q2 / (y * s_hadronic)
    u = x_b / xp
    coeff = herwig_channel_coefficients(flavor, q2, inputs, channel)
    coeff_cc = herwig_channel_coefficients(conjugate_flavor(flavor), q2, inputs, channel)
    reduced_cf = 4.0 / 3.0 * herwig_alpha_s / (2.0 * math.pi)
    reduced_tr = 0.5 * herwig_alpha_s / (2.0 * math.pi)
    log_ratio = math.log((1.0 - xp) / xp)

    lo_pdf = xfx_over_x(sum_pdf, flavor.lhapdf_id, x_b, q2)
    if abs(lo_pdf) <= 1.0e-30:
        raise ZeroDivisionError("Encountered vanishing LO sum PDF at x_B.")

    q_pdf = xfx_over_x(sum_pdf, flavor.lhapdf_id, u, q2)
    g_pdf = xfx_over_x(sum_pdf, 21, u, q2)
    dlo_pdf = xfx_over_x(diff_pdf, flavor.lhapdf_id, x_b, q2)
    dq_pdf = xfx_over_x(diff_pdf, flavor.lhapdf_id, u, q2)
    dg_pdf = xfx_over_x(diff_pdf, 21, u, q2)

    q_ratio = q_pdf / lo_pdf
    g_ratio = g_pdf / lo_pdf
    ratio_floor = 1.0e-12
    min_dlo = max(ratio_floor, 1.0e-4 * abs(lo_pdf))

    terms: Dict[str, HerwigHelicityTerms] = {}
    for hel, pl, pz in HELICITY_COMBINATIONS:
        pq = clamp_pol(pz * dlo_pdf / lo_pdf) if abs(lo_pdf) > 1.0e-30 else 0.0
        pq_m = clamp_pol(pz * dq_pdf / q_pdf) if abs(q_pdf) > ratio_floor else 0.0
        pg_m = clamp_pol(pz * dg_pdf / g_pdf) if abs(g_pdf) > ratio_floor else 0.0
        deltaq_over_lo = pz * dq_pdf / lo_pdf if abs(lo_pdf) > ratio_floor else 0.0
        deltag_over_lo = pz * dg_pdf / lo_pdf if abs(lo_pdf) > ratio_floor else 0.0
        dq_ratio = dq_pdf / dlo_pdf if abs(dlo_pdf) > min_dlo else 0.0
        dg_ratio = dg_pdf / dlo_pdf if abs(dlo_pdf) > min_dlo else 0.0

        d_born = real_emission_denominator(coeff, pl, pq)
        n_born = real_emission_numerator(coeff, pl, pq)
        d_mapped = real_emission_denominator(coeff, pl, pq_m)
        n_mapped = real_emission_numerator(coeff, pl, pq_m)
        d_g_r2 = real_emission_denominator(coeff, pl, pg_m)
        n_g_r2 = real_emission_numerator(coeff, pl, pg_m)
        d_g_r3 = real_emission_denominator(coeff_cc, pl, -pg_m)
        n_g_r3 = real_emission_numerator(coeff_cc, pl, -pg_m)
        a_born = n_born / d_born if abs(d_born) > 1.0e-30 else 0.0
        a_q_mapped = n_mapped / d_mapped if abs(d_mapped) > 1.0e-30 else 0.0
        a_g_r2 = n_g_r2 / d_g_r2 if abs(d_g_r2) > 1.0e-30 else 0.0
        a_g_r3 = n_g_r3 / d_g_r3 if abs(d_g_r3) > 1.0e-30 else 0.0
        audit = nc_audit_payload(coeff, ell, pl, pq, pq_m)
        qcdc_den_ratio = audit.mapped_denominator_ratio
        born_kernel = herwig_kernel_prefactor(ell) * sigma_weight(coeff, ell, pl, pq)
        born_abs = lo_pdf * born_kernel

        q_u = audit.q_unpolarized
        q_p = audit.q_polarized
        g_u = audit.g_unpolarized
        g_p = audit.g_polarized
        q_odd_response = audit.q_odd_response
        g_odd_response = audit.g_odd_response

        virt_correction = reduced_cf * (
            -4.5
            - math.pi * math.pi / 3.0
            + 1.5 * math.log(1.0 / (1.0 - x_b))
            + 2.0 * math.log(1.0 - x_b) * 0.0
            + math.log(1.0 - x_b) ** 2
        )
        collg_unpol = reduced_tr / xp * g_ratio * (
            2.0 * xp * (1.0 - xp) + (xp * xp + (1.0 - xp) ** 2) * log_ratio
        )
        collg_odd = g_odd_response * reduced_tr / xp * deltag_over_lo * (
            2.0 * (1.0 - xp) + (2.0 * xp - 1.0) * log_ratio
        )

        k1 = (
            1.0
            - xp
            - 2.0 / (1.0 - xp) * math.log(xp)
            - (1.0 + xp) * math.log((1.0 - xp) / xp)
        )
        k2 = 2.0 / (1.0 - xp) * math.log(1.0 - xp) - 1.5 / (1.0 - xp)
        collq_k1_unpol = reduced_cf / xp * q_ratio * k1
        collq_k2_unpol = reduced_cf / xp * (q_ratio - xp) * k2
        collq_unpol = collq_k1_unpol + collq_k2_unpol
        collq_odd = q_odd_response * (
            reduced_cf / xp * deltaq_over_lo * k1
            + reduced_cf / xp * (deltaq_over_lo - pq * xp) * k2
        )

        born_factor = 1.0 + a_born * ell + ell * ell
        qcdc_even = (
            qcdc_den_ratio
            * reduced_cf
            / xp
            / born_factor
            * q_ratio
            * (2.0 + 2.0 * ell * ell - xp + 3.0 * xp * ell * ell)
        )
        qcdc_aq = (
            qcdc_den_ratio
            * reduced_cf
            / xp
            / born_factor
            * q_ratio
            * (a_q_mapped * ell * (2.0 * xp + 1.0))
        )
        qcdc_real = qcdc_even + qcdc_aq
        bgf_split_unpol = g_u * (
            -reduced_tr / xp * g_ratio * ((1.0 + ell * ell + 2.0 * (1.0 - 3.0 * ell * ell) * xp * (1.0 - xp)) / (1.0 + ell * ell))
        )
        bgf_split_odd = g_odd_response * (-2.0 * reduced_tr / xp * deltag_over_lo * (2.0 * xp - 1.0))
        bgf_real = bgf_split_unpol + bgf_split_odd
        bgf_legacy_unpol = -reduced_tr / xp * g_ratio * (
            1.0 + ell * ell + 2.0 * (1.0 - 3.0 * ell * ell) * xp * (1.0 - xp)
        )
        bgf_legacy_spin = -reduced_tr / xp * g_ratio * (
            2.0 * ell * (xp * xp * a_g_r2 + (1.0 - xp) * (1.0 - xp) * a_g_r3)
        )
        bgf_legacy_total = bgf_legacy_unpol + bgf_legacy_spin

        terms[hel] = HerwigHelicityTerms(
            born_abs=born_abs,
            virt_corr=virt_correction,
            quark_collinear=q_u * collq_unpol + collq_odd,
            collq_k1=q_u * collq_k1_unpol,
            collq_k2=q_u * collq_k2_unpol,
            collq_odd=collq_odd,
            gluon_collinear=g_u * collg_unpol + collg_odd,
            qcdc_mapped=qcdc_real,
            qcdc_even=qcdc_even,
            qcdc_aq=qcdc_aq,
            bgf_mapped=bgf_real,
            Pq=pq,
            Pq_m=pq_m,
            Pg_m=pg_m,
            dq_ratio=dq_ratio,
            dg_ratio=dg_ratio,
            deltaq_over_lo=deltaq_over_lo,
            deltag_over_lo=deltag_over_lo,
            q_odd_response=q_odd_response,
            g_odd_response=g_odd_response,
            a_born=a_born,
            a_q_mapped=a_q_mapped,
            a_g_r2=a_g_r2,
            a_g_r3=a_g_r3,
            bgf_split_unpol=bgf_split_unpol,
            bgf_split_odd=bgf_split_odd,
            bgf_legacy_unpol=bgf_legacy_unpol,
            bgf_legacy_spin=bgf_legacy_spin,
            bgf_legacy_total=bgf_legacy_total,
        )
    return terms


def ll_combine(values: Dict[str, float]) -> float:
    return 0.25 * (values["PP"] + values["MM"] - values["PM"] - values["MP"])


def aggregate_herwig_terms(
    flavor: Flavor,
    q2: float,
    y: float,
    xp: float,
    s_hadronic: float,
    channel: str,
    inputs: EWInputs,
    sum_pdf,
    diff_pdf,
    herwig_alpha_s: float,
) -> AggregatedTerms:
    hel_terms = evaluate_herwig_helicity_terms(
        flavor, q2, y, xp, s_hadronic, channel, inputs, sum_pdf, diff_pdf, herwig_alpha_s
    )
    born_values = {hel: term.born_abs for hel, term in hel_terms.items()}
    born_abs = ll_combine(born_values)
    if abs(born_abs) <= 1.0e-30:
        raise ZeroDivisionError("Encountered vanishing sigma_LL Born kernel.")

    quark_values = {
        hel: term.born_abs * (term.virt_corr + term.quark_collinear)
        for hel, term in hel_terms.items()
    }
    virtual_values = {hel: term.born_abs * term.virt_corr for hel, term in hel_terms.items()}
    collq_values = {hel: term.born_abs * term.quark_collinear for hel, term in hel_terms.items()}
    collq_k1_values = {hel: term.born_abs * term.collq_k1 for hel, term in hel_terms.items()}
    collq_k2_values = {hel: term.born_abs * term.collq_k2 for hel, term in hel_terms.items()}
    collq_odd_values = {hel: term.born_abs * term.collq_odd for hel, term in hel_terms.items()}
    gluon_coll_values = {
        hel: term.born_abs * term.gluon_collinear for hel, term in hel_terms.items()
    }
    qcdc_values = {
        hel: term.born_abs * term.qcdc_mapped for hel, term in hel_terms.items()
    }
    qcdc_even_values = {hel: term.born_abs * term.qcdc_even for hel, term in hel_terms.items()}
    qcdc_aq_values = {hel: term.born_abs * term.qcdc_aq for hel, term in hel_terms.items()}
    bgf_values = {
        hel: term.born_abs * term.bgf_mapped for hel, term in hel_terms.items()
    }
    # Keep the singular quark bucket explicitly split so the absolute audit can
    # check `virtual_abs + collq_abs == quark_collinear_abs`.
    virtual_abs = ll_combine(virtual_values)
    collq_abs = ll_combine(collq_values)
    # And keep the non-delta singular remainder split one level deeper so the
    # control-point audit can check `collq_k1_abs + collq_k2_abs + collq_odd_abs == collq_abs`.
    collq_k1_abs = ll_combine(collq_k1_values)
    collq_k2_abs = ll_combine(collq_k2_values)
    collq_odd_abs = ll_combine(collq_odd_values)
    quark_collinear_abs = ll_combine(quark_values)
    gluon_collinear_abs = ll_combine(gluon_coll_values)
    # Likewise for the mapped QCDC term:
    # `qcdc_even_abs + qcdc_aq_abs == qcdc_mapped_abs`.
    qcdc_mapped_abs = ll_combine(qcdc_values)
    qcdc_even_abs = ll_combine(qcdc_even_values)
    qcdc_aq_abs = ll_combine(qcdc_aq_values)
    bgf_mapped_abs = ll_combine(bgf_values)

    first = hel_terms["PP"]
    return AggregatedTerms(
        born_abs=born_abs,
        virtual_abs=virtual_abs,
        collq_abs=collq_abs,
        collq_k1_abs=collq_k1_abs,
        collq_k2_abs=collq_k2_abs,
        collq_odd_abs=collq_odd_abs,
        quark_collinear=quark_collinear_abs / born_abs,
        gluon_collinear=gluon_collinear_abs / born_abs,
        qcdc_mapped=qcdc_mapped_abs / born_abs,
        bgf_mapped=bgf_mapped_abs / born_abs,
        quark_collinear_abs=quark_collinear_abs,
        gluon_collinear_abs=gluon_collinear_abs,
        qcdc_mapped_abs=qcdc_mapped_abs,
        qcdc_even_abs=qcdc_even_abs,
        qcdc_aq_abs=qcdc_aq_abs,
        bgf_mapped_abs=bgf_mapped_abs,
        nlo_sum_abs=quark_collinear_abs + gluon_collinear_abs + qcdc_mapped_abs + bgf_mapped_abs,
        inputs={
            "Pq": first.Pq,
            "Pq_m": first.Pq_m,
            "Pg_m": first.Pg_m,
            "dqRatio": first.dq_ratio,
            "dgRatio": first.dg_ratio,
            "deltaqOverLo": first.deltaq_over_lo,
            "deltagOverLo": first.deltag_over_lo,
            "qOddResponse": first.q_odd_response,
            "gOddResponse": first.g_odd_response,
        },
    )


def poldis_coefficients_single_flavor(
    q2: float,
    x_b: float,
    y: float,
    xp: float,
    channel: str,
    flavor: Flavor,
    inputs: EWInputs,
    diff_pdf,
    poldis_alpha_s: float,
) -> AggregatedTerms:
    z = xp
    u = x_b / z
    pieces = poldis_channel_pieces(flavor, q2, inputs, channel)
    vec = pieces.vec
    ax = pieces.ax
    ym = 1.0 - (1.0 - y) ** 2
    yp = 1.0 + (1.0 - y) ** 2

    f1 = diff_pdf.xfxQ2(flavor.lhapdf_id, x_b, q2)
    fu = diff_pdf.xfxQ2(flavor.lhapdf_id, u, q2)
    gu = diff_pdf.xfxQ2(21, u, q2)

    sigma1 = f1 * vec
    sigma = fu * vec
    sigmaw1 = f1 * ax
    sigmaw = fu * ax
    gl = gu * vec

    g2_lo_raw = x_b * sigma1
    g4_lo_raw = x_b * sigmaw1
    born_abs_raw = ym * g2_lo_raw + yp * g4_lo_raw
    if abs(born_abs_raw) <= 1.0e-30:
        raise ZeroDivisionError("Encountered vanishing POLDIS sigma_LL Born kernel.")

    logmu = 0.0
    c2qs_reg1 = (4.0 / 3.0) * (
        -2.0 * (1.0 + z) * math.log(1.0 - z)
        - 2.0 * (1.0 + z * z) / (1.0 - z) * math.log(z)
        + 4.0
        + 2.0 * z
        - 2.0 * (1.0 + z) * logmu
    )
    c2qs_sing1 = (4.0 / 3.0) * (4.0 * math.log(1.0 - z) / (1.0 - z) - 3.0 / (1.0 - z) + 4.0 / (1.0 - z) * logmu)
    c2qs_delta1 = (4.0 / 3.0) * (
        -4.0 * math.pi * math.pi / 6.0
        - 9.0
        + 3.0 * logmu
        - 3.0 * math.log(1.0 - x_b)
        + 2.0 * math.log(1.0 - x_b) ** 2
        + 4.0 * math.log(1.0 - x_b) * logmu
    )
    c4qs_reg1 = (4.0 / 3.0) * (
        -2.0 * (1.0 + z) * math.log(1.0 - z)
        - 2.0 * (1.0 + z * z) / (1.0 - z) * math.log(z)
        + 6.0
        + 4.0 * z
        - 2.0 * (1.0 + z) * logmu
    )
    c4qs_sing1 = c2qs_sing1
    c4qs_delta1 = c2qs_delta1
    clqs1 = (4.0 / 3.0) * 4.0 * z
    c2g1 = 1.0 * 0.5 * (
        4.0 * (2.0 * z - 1.0) * (math.log(1.0 - z) - math.log(z))
        + 4.0 * (3.0 - 4.0 * z)
        + 4.0 * (2.0 * z - 1.0) * logmu
    )
    a2pi = poldis_alpha_s / (2.0 * math.pi)

    quark_regular_abs_raw = x_b * (a2pi / 2.0) * (
        ym * (c2qs_reg1 * sigma / z)
        + yp * (c4qs_reg1 * sigmaw / z)
        - y * y * (clqs1 * sigmaw / z)
    )
    qcdc_c4_minus_c2_abs_raw = x_b * (a2pi / 2.0) * (yp * ((c4qs_reg1 - c2qs_reg1) * sigmaw / z))
    virtual_abs_raw = x_b * (a2pi / 2.0) * (
        ym * (sigma1 * c2qs_delta1) + yp * (sigmaw1 * c4qs_delta1)
    )
    collq_abs_raw = x_b * (a2pi / 2.0) * (
        ym * (c2qs_sing1 * (sigma / z - sigma1))
        + yp * (c4qs_sing1 * (sigmaw / z - sigmaw1))
    )
    quark_singular_abs_raw = virtual_abs_raw + collq_abs_raw
    gluon_regular_abs_raw = x_b * (a2pi / 2.0) * (ym * gl * c2g1 / z)

    # The inclusive POLDIS analogue is built from x f(x)-style inputs, so the
    # raw absolute pieces carry a common extra x_B^2 relative to the Herwig
    # f(x)-normalized kernel. Remove only that overall factor here; the
    # normalized ratios below remain unchanged.
    absolute_rescale = 1.0 / (x_b * x_b)
    born_abs = born_abs_raw * absolute_rescale
    virtual_abs = virtual_abs_raw * absolute_rescale
    collq_abs = collq_abs_raw * absolute_rescale
    # The inclusive reference exposes only the total non-delta singular remainder,
    # not a one-to-one `k1/k2/odd` local split.
    collq_k1_abs = math.nan
    collq_k2_abs = math.nan
    collq_odd_abs = math.nan
    quark_singular_abs = quark_singular_abs_raw * absolute_rescale
    quark_regular_abs = quark_regular_abs_raw * absolute_rescale
    qcdc_c4_minus_c2_abs = qcdc_c4_minus_c2_abs_raw * absolute_rescale
    gluon_regular_abs = gluon_regular_abs_raw * absolute_rescale
    nlo_sum_abs = quark_singular_abs + quark_regular_abs + gluon_regular_abs

    return AggregatedTerms(
        born_abs=born_abs,
        virtual_abs=virtual_abs,
        collq_abs=collq_abs,
        collq_k1_abs=collq_k1_abs,
        collq_k2_abs=collq_k2_abs,
        collq_odd_abs=collq_odd_abs,
        quark_collinear=quark_singular_abs / born_abs,
        gluon_collinear=0.0,
        qcdc_mapped=quark_regular_abs / born_abs,
        bgf_mapped=gluon_regular_abs / born_abs,
        quark_collinear_abs=quark_singular_abs,
        gluon_collinear_abs=0.0,
        qcdc_mapped_abs=quark_regular_abs,
        qcdc_even_abs=math.nan,
        qcdc_aq_abs=math.nan,
        bgf_mapped_abs=gluon_regular_abs,
        nlo_sum_abs=nlo_sum_abs,
        inputs={
            "qcdc_c4_minus_c2_abs": qcdc_c4_minus_c2_abs,
            "qcdc_c4_minus_c2": qcdc_c4_minus_c2_abs / born_abs,
        },
    )


def format_float(value: float) -> str:
    return f"{value:.8e}"


def format_ratio(herwig: float, poldis: float) -> str:
    if abs(poldis) <= 1.0e-30:
        return "nan"
    return f"{herwig / poldis:.8e}"


def build_comparison_row(
    flavor_label: str,
    x_b: float,
    xp: float,
    term: str,
    herwig_value: float,
    poldis_value: float,
) -> List[str]:
    return [
        flavor_label,
        f"{x_b:.6f}",
        f"{xp:.6f}",
        term,
        format_float(herwig_value),
        format_float(poldis_value),
        format_ratio(herwig_value, poldis_value),
        format_float(herwig_value - poldis_value),
    ]


def render_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    widths = [len(header) for header in headers]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))
    sep = "  "
    head = sep.join(header.ljust(widths[idx]) for idx, header in enumerate(headers))
    rule = sep.join("-" * widths[idx] for idx in range(len(headers)))
    body = [sep.join(value.ljust(widths[idx]) for idx, value in enumerate(row)) for row in rows]
    return "\n".join([head, rule, *body])


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    lhapdf = import_lhapdf()
    sum_pdf = lhapdf.mkPDF(args.sum_pdf, 0)
    diff_pdf = lhapdf.mkPDF(args.diff_pdf, 0)

    ew_inputs = EWInputs(
        sin2_theta_w=args.sin2_theta_w,
        mz=args.mz,
        gamma_z=args.gamma_z,
        lepton_charge=-1.0,
    )
    quark_masses2 = default_matchbox_quark_masses2()
    lambda_squared = solve_matchbox_lambda_squared(args.input_alpha_s, args.input_scale)
    selected_flavors = resolve_flavors(args.flavor)

    max_disagreement = {
        "quark_collinear": (0.0, ""),
        "gluon_collinear": (0.0, ""),
        "qcdc_mapped": (0.0, ""),
        "qcdc_mapped_plus_c4diag": (0.0, ""),
        "bgf_mapped": (0.0, ""),
    }
    max_absolute_disagreement = {
        "virtual_abs": (0.0, ""),
        "collq_abs": (0.0, ""),
        "quark_collinear_abs": (0.0, ""),
        "gluon_collinear_abs": (0.0, ""),
        "qcdc_mapped_abs": (0.0, ""),
        "qcdc_mapped_abs_plus_c4diag": (0.0, ""),
        "qcdc_c4_minus_c2_abs": (0.0, ""),
        "bgf_mapped_abs": (0.0, ""),
        "nlo_sum_abs": (0.0, ""),
        "collq_ref_remainder_abs": (0.0, ""),
        "remainder_minus_odd": (0.0, ""),
        "qcdc_ref_remainder_abs": (0.0, ""),
        "remainder_minus_aq": (0.0, ""),
    }

    print(
        "Polarized NLO kernel audit\n"
        f"  channel = {args.channel}\n"
        f"  sum PDF = {args.sum_pdf}\n"
        f"  diff/POLDIS PDF = {args.diff_pdf}\n"
        f"  sin^2(theta_W) = {args.sin2_theta_w:.18f}\n"
        f"  sqrt(s) = {math.sqrt(args.s_hadronic):.6f} GeV"
    )

    for q2 in args.q2_values:
        herwig_alpha_s = matchbox_running_alpha(q2, lambda_squared, quark_masses2)
        poldis_alpha_s = diff_pdf.alphasQ2(q2)
        print(
            f"\nQ^2 = {q2:.3f} GeV^2"
            f"  alpha_s^Herwig = {herwig_alpha_s:.8f}"
            f"  alpha_s^POLDIS = {poldis_alpha_s:.8f}"
        )
        for y in args.y_values:
            x_b = q2 / (y * args.s_hadronic)
            if not (0.0 < x_b < 1.0):
                print(f"\n  Skipping y={y:.3f}: x_B={x_b:.6f} is outside (0,1).")
                continue
            xp_values = [x_b + fraction * (1.0 - x_b) for fraction in args.xp_fractions]
            for xp in xp_values:
                rows: List[List[str]] = []
                absolute_rows: List[List[str]] = []
                bgf_debug_rows: List[List[str]] = []
                collq_debug_rows: List[List[str]] = []
                qcdc_debug_rows: List[List[str]] = []
                qcdc_remainder_rows: List[List[str]] = []
                qcdc_c4_diag_rows: List[List[str]] = []
                for flavor in selected_flavors:
                    herwig = aggregate_herwig_terms(
                        flavor,
                        q2,
                        y,
                        xp,
                        args.s_hadronic,
                        args.channel,
                        ew_inputs,
                        sum_pdf,
                        diff_pdf,
                        herwig_alpha_s,
                    )
                    poldis = poldis_coefficients_single_flavor(
                        q2, x_b, y, xp, args.channel, flavor, ew_inputs, diff_pdf, poldis_alpha_s
                    )
                    pp_terms = None
                    if args.bgf_debug:
                        pp_terms = evaluate_herwig_helicity_terms(
                            flavor,
                            q2,
                            y,
                            xp,
                            args.s_hadronic,
                            args.channel,
                            ew_inputs,
                            sum_pdf,
                            diff_pdf,
                            herwig_alpha_s,
                        )["PP"]
                    rows.extend(
                        [
                            [
                                flavor.label,
                                f"{x_b:.6f}",
                                f"{xp:.6f}",
                                "Pq",
                                format_float(herwig.inputs["Pq"]),
                                format_float(herwig.inputs["Pq"]),
                                "1.00000000e+00",
                                "0.00000000e+00",
                            ],
                            [
                                flavor.label,
                                f"{x_b:.6f}",
                                f"{xp:.6f}",
                                "Pq_m",
                                format_float(herwig.inputs["Pq_m"]),
                                format_float(herwig.inputs["Pq_m"]),
                                "1.00000000e+00",
                                "0.00000000e+00",
                            ],
                            [
                                flavor.label,
                                f"{x_b:.6f}",
                                f"{xp:.6f}",
                                "Pg_m",
                                format_float(herwig.inputs["Pg_m"]),
                                format_float(herwig.inputs["Pg_m"]),
                                "1.00000000e+00",
                                "0.00000000e+00",
                            ],
                            [
                                flavor.label,
                                f"{x_b:.6f}",
                                f"{xp:.6f}",
                                "dqRatio",
                                format_float(herwig.inputs["dqRatio"]),
                                format_float(herwig.inputs["dqRatio"]),
                                "1.00000000e+00",
                                "0.00000000e+00",
                            ],
                            [
                                flavor.label,
                                f"{x_b:.6f}",
                                f"{xp:.6f}",
                                "dgRatio",
                                format_float(herwig.inputs["dgRatio"]),
                                format_float(herwig.inputs["dgRatio"]),
                                "1.00000000e+00",
                                "0.00000000e+00",
                            ],
                            [
                                flavor.label,
                                f"{x_b:.6f}",
                                f"{xp:.6f}",
                                "deltaqOverLo",
                                format_float(herwig.inputs["deltaqOverLo"]),
                                format_float(herwig.inputs["deltaqOverLo"]),
                                "1.00000000e+00",
                                "0.00000000e+00",
                            ],
                            [
                                flavor.label,
                                f"{x_b:.6f}",
                                f"{xp:.6f}",
                                "deltagOverLo",
                                format_float(herwig.inputs["deltagOverLo"]),
                                format_float(herwig.inputs["deltagOverLo"]),
                                "1.00000000e+00",
                                "0.00000000e+00",
                            ],
                            [
                                flavor.label,
                                f"{x_b:.6f}",
                                f"{xp:.6f}",
                                "qOddResponse",
                                format_float(herwig.inputs["qOddResponse"]),
                                format_float(herwig.inputs["qOddResponse"]),
                                "1.00000000e+00",
                                "0.00000000e+00",
                            ],
                            [
                                flavor.label,
                                f"{x_b:.6f}",
                                f"{xp:.6f}",
                                "gOddResponse",
                                format_float(herwig.inputs["gOddResponse"]),
                                format_float(herwig.inputs["gOddResponse"]),
                                "1.00000000e+00",
                                "0.00000000e+00",
                            ],
                            build_comparison_row(
                                flavor.label,
                                x_b,
                                xp,
                                "quark_collinear",
                                herwig.quark_collinear,
                                poldis.quark_collinear,
                            ),
                            build_comparison_row(
                                flavor.label,
                                x_b,
                                xp,
                                "gluon_collinear",
                                herwig.gluon_collinear,
                                poldis.gluon_collinear,
                            ),
                            build_comparison_row(
                                flavor.label,
                                x_b,
                                xp,
                                "qcdc_mapped",
                                herwig.qcdc_mapped,
                                poldis.qcdc_mapped,
                            ),
                            build_comparison_row(
                                flavor.label,
                                x_b,
                                xp,
                                "bgf_mapped",
                                herwig.bgf_mapped,
                                poldis.bgf_mapped,
                            ),
                        ]
                    )
                    qcdc_c4_diag = poldis.inputs.get("qcdc_c4_minus_c2", 0.0)
                    qcdc_c4_diag_abs = poldis.inputs.get("qcdc_c4_minus_c2_abs", 0.0)
                    herwig_qcdc_plus_c4diag = herwig.qcdc_mapped + qcdc_c4_diag
                    herwig_qcdc_plus_c4diag_abs = herwig.qcdc_mapped_abs + qcdc_c4_diag_abs
                    if args.c4_diagnostic:
                        rows.append(
                            build_comparison_row(
                                flavor.label,
                                x_b,
                                xp,
                                "qcdc_mapped_plus_c4diag",
                                herwig_qcdc_plus_c4diag,
                                poldis.qcdc_mapped,
                            )
                        )
                    if args.absolute_debug:
                        absolute_rows.extend(
                            [
                                build_comparison_row(
                                    flavor.label, x_b, xp, "born_abs", herwig.born_abs, poldis.born_abs
                                ),
                                build_comparison_row(
                                    flavor.label,
                                    x_b,
                                    xp,
                                    "virtual_abs",
                                    herwig.virtual_abs,
                                    poldis.virtual_abs,
                                ),
                                build_comparison_row(
                                    flavor.label,
                                    x_b,
                                    xp,
                                    "collq_abs",
                                    herwig.collq_abs,
                                    poldis.collq_abs,
                                ),
                                build_comparison_row(
                                    flavor.label,
                                    x_b,
                                    xp,
                                    "quark_collinear_abs",
                                    herwig.quark_collinear_abs,
                                    poldis.quark_collinear_abs,
                                ),
                                build_comparison_row(
                                    flavor.label,
                                    x_b,
                                    xp,
                                    "gluon_collinear_abs",
                                    herwig.gluon_collinear_abs,
                                    poldis.gluon_collinear_abs,
                                ),
                                build_comparison_row(
                                    flavor.label,
                                    x_b,
                                    xp,
                                    "qcdc_mapped_abs",
                                    herwig.qcdc_mapped_abs,
                                    poldis.qcdc_mapped_abs,
                                ),
                                *(
                                    [
                                        build_comparison_row(
                                            flavor.label,
                                            x_b,
                                            xp,
                                            "qcdc_mapped_abs_plus_c4diag",
                                            herwig_qcdc_plus_c4diag_abs,
                                            poldis.qcdc_mapped_abs,
                                        ),
                                        build_comparison_row(
                                            flavor.label,
                                            x_b,
                                            xp,
                                            "qcdc_c4_minus_c2_abs",
                                            qcdc_c4_diag_abs,
                                            qcdc_c4_diag_abs,
                                        ),
                                    ]
                                    if args.c4_diagnostic
                                    else []
                                ),
                                build_comparison_row(
                                    flavor.label,
                                    x_b,
                                    xp,
                                    "bgf_mapped_abs",
                                    herwig.bgf_mapped_abs,
                                    poldis.bgf_mapped_abs,
                                ),
                                build_comparison_row(
                                    flavor.label,
                                    x_b,
                                    xp,
                                    "nlo_sum_abs",
                                    herwig.nlo_sum_abs,
                                    poldis.nlo_sum_abs,
                                ),
                            ]
                        )
                        collq_even_abs = herwig.collq_k1_abs + herwig.collq_k2_abs
                        collq_split_minus_total = collq_even_abs + herwig.collq_odd_abs - herwig.collq_abs
                        collq_ref_remainder_abs = poldis.collq_abs - collq_even_abs
                        remainder_minus_odd = collq_ref_remainder_abs - herwig.collq_odd_abs
                        collq_debug_rows.append(
                            [
                                flavor.label,
                                f"{x_b:.6f}",
                                f"{xp:.6f}",
                                format_float(collq_even_abs),
                                format_float(herwig.collq_odd_abs),
                                format_float(poldis.collq_abs),
                                format_float(collq_ref_remainder_abs),
                                format_float(remainder_minus_odd),
                                format_float(collq_split_minus_total),
                            ]
                        )
                        qcdc_debug_rows.append(
                            [
                                flavor.label,
                                format_float(herwig.qcdc_even_abs),
                                format_float(herwig.qcdc_aq_abs),
                                format_float(herwig.qcdc_mapped_abs),
                                format_float(herwig.qcdc_even_abs + herwig.qcdc_aq_abs - herwig.qcdc_mapped_abs),
                            ]
                        )
                        qcdc_ref_remainder_abs = poldis.qcdc_mapped_abs - herwig.qcdc_even_abs
                        remainder_minus_aq = qcdc_ref_remainder_abs - herwig.qcdc_aq_abs
                        qcdc_remainder_rows.append(
                            [
                                flavor.label,
                                f"{x_b:.6f}",
                                f"{xp:.6f}",
                                format_float(herwig.qcdc_even_abs),
                                format_float(herwig.qcdc_aq_abs),
                                format_float(poldis.qcdc_mapped_abs),
                                format_float(qcdc_ref_remainder_abs),
                                format_float(remainder_minus_aq),
                            ]
                        )
                        if args.c4_diagnostic:
                            residual_before = herwig.qcdc_mapped_abs - poldis.qcdc_mapped_abs
                            residual_after = herwig_qcdc_plus_c4diag_abs - poldis.qcdc_mapped_abs
                            qcdc_c4_diag_rows.append(
                                [
                                    flavor.label,
                                    f"{x_b:.6f}",
                                    f"{xp:.6f}",
                                    format_float(herwig.qcdc_mapped_abs),
                                    format_float(qcdc_c4_diag_abs),
                                    format_float(herwig_qcdc_plus_c4diag_abs),
                                    format_float(poldis.qcdc_mapped_abs),
                                    format_float(residual_before),
                                    format_float(residual_after),
                                ]
                            )
                    for key, herwig_value, poldis_value in (
                        ("quark_collinear", herwig.quark_collinear, poldis.quark_collinear),
                        ("gluon_collinear", herwig.gluon_collinear, poldis.gluon_collinear),
                        ("qcdc_mapped", herwig.qcdc_mapped, poldis.qcdc_mapped),
                        ("qcdc_mapped_plus_c4diag", herwig_qcdc_plus_c4diag, poldis.qcdc_mapped),
                        ("bgf_mapped", herwig.bgf_mapped, poldis.bgf_mapped),
                    ):
                        diff = abs(herwig_value - poldis_value)
                        if diff > max_disagreement[key][0]:
                            max_disagreement[key] = (
                                diff,
                                f"flavor={flavor.label}, Q2={q2:.3f}, y={y:.3f}, xp={xp:.6f}",
                            )
                    for key, herwig_value, poldis_value in (
                        ("virtual_abs", herwig.virtual_abs, poldis.virtual_abs),
                        ("collq_abs", herwig.collq_abs, poldis.collq_abs),
                        ("quark_collinear_abs", herwig.quark_collinear_abs, poldis.quark_collinear_abs),
                        ("gluon_collinear_abs", herwig.gluon_collinear_abs, poldis.gluon_collinear_abs),
                        ("qcdc_mapped_abs", herwig.qcdc_mapped_abs, poldis.qcdc_mapped_abs),
                        ("qcdc_mapped_abs_plus_c4diag", herwig_qcdc_plus_c4diag_abs, poldis.qcdc_mapped_abs),
                        ("qcdc_c4_minus_c2_abs", qcdc_c4_diag_abs, 0.0),
                        ("bgf_mapped_abs", herwig.bgf_mapped_abs, poldis.bgf_mapped_abs),
                        ("nlo_sum_abs", herwig.nlo_sum_abs, poldis.nlo_sum_abs),
                    ):
                        diff = abs(herwig_value - poldis_value)
                        if diff > max_absolute_disagreement[key][0]:
                            max_absolute_disagreement[key] = (
                                diff,
                                f"flavor={flavor.label}, Q2={q2:.3f}, y={y:.3f}, xp={xp:.6f}",
                            )
                    if args.absolute_debug:
                        for key, value in (
                            ("collq_ref_remainder_abs", collq_ref_remainder_abs),
                            ("remainder_minus_odd", remainder_minus_odd),
                            ("qcdc_ref_remainder_abs", qcdc_ref_remainder_abs),
                            ("remainder_minus_aq", remainder_minus_aq),
                        ):
                            diff = abs(value)
                            if diff > max_absolute_disagreement[key][0]:
                                max_absolute_disagreement[key] = (
                                    diff,
                                    f"flavor={flavor.label}, Q2={q2:.3f}, y={y:.3f}, xp={xp:.6f}",
                                )
                    if pp_terms is not None:
                        bgf_debug_rows.append(
                            [
                                flavor.label,
                                format_float(pp_terms.a_born),
                                format_float(pp_terms.a_q_mapped),
                                format_float(pp_terms.a_g_r2),
                                format_float(pp_terms.a_g_r3),
                                format_float(pp_terms.bgf_split_unpol),
                                format_float(pp_terms.bgf_split_odd),
                                format_float(pp_terms.bgf_mapped),
                                format_float(pp_terms.bgf_legacy_unpol),
                                format_float(pp_terms.bgf_legacy_spin),
                                format_float(pp_terms.bgf_legacy_total),
                                format_float(pp_terms.bgf_mapped - pp_terms.bgf_legacy_total),
                            ]
                        )

                print(f"\n  y = {y:.3f}, xp = {xp:.6f}, x_B = {x_b:.6f}")
                print(
                    render_table(
                        ["flavor", "x_B", "xp", "term", "Herwig", "POLDIS", "ratio", "abs_diff"],
                        rows,
                    )
                )
                if absolute_rows:
                    print("\n  Absolute sigma_LL contributions")
                    print(
                        render_table(
                            ["flavor", "x_B", "xp", "term", "Herwig", "POLDIS", "ratio", "abs_diff"],
                            absolute_rows,
                        )
                    )
                if collq_debug_rows:
                    print("\n  Collq remainder diagnostic")
                    print(
                        render_table(
                            [
                                "flavor",
                                "x_B",
                                "xp",
                                "collqEvenAbs",
                                "collqOddAbs",
                                "poldisCollqTotalAbs",
                                "collqRefRemainderAbs",
                                "remainderMinusOdd",
                                "splitMinusTotal",
                            ],
                            collq_debug_rows,
                        )
                    )
                if qcdc_debug_rows:
                    print("\n  Herwig QCDC diagnostic (LL combined)")
                    print(
                        render_table(
                            ["flavor", "qcdcEvenAbs", "qcdcAqAbs", "qcdcMappedAbs", "splitMinusTotal"],
                            qcdc_debug_rows,
                        )
                    )
                if qcdc_remainder_rows:
                    print("\n  QCDC remainder diagnostic")
                    print(
                        render_table(
                            [
                                "flavor",
                                "x_B",
                                "xp",
                                "qcdcEvenAbs",
                                "qcdcAqAbs",
                                "poldisQcdcTotalAbs",
                                "qcdcRefRemainderAbs",
                                "remainderMinusAq",
                            ],
                            qcdc_remainder_rows,
                        )
                    )
                if qcdc_c4_diag_rows:
                    print("\n  QCDC c4-c2 diagnostic")
                    print(
                        render_table(
                            [
                                "flavor",
                                "x_B",
                                "xp",
                                "herwigQcdcAbs",
                                "c4MinusC2DiagAbs",
                                "herwigPlusDiagAbs",
                                "poldisQcdcAbs",
                                "residualBefore",
                                "residualAfter",
                            ],
                            qcdc_c4_diag_rows,
                        )
                    )
                if bgf_debug_rows:
                    print("\n  Herwig BGF diagnostic (PP helicity)")
                    print(
                        render_table(
                            [
                                "flavor",
                                "aBorn",
                                "aQMapped",
                                "aGR2",
                                "aGR3",
                                "bgfSplitUnpol",
                                "bgfSplitOdd",
                                "bgfSplitTotal",
                                "bgfLegacyUnpol",
                                "bgfLegacySpin",
                                "bgfLegacyTotal",
                                "splitMinusLegacy",
                            ],
                            bgf_debug_rows,
                        )
                    )

    print("\nLargest absolute disagreements by category")
    summary_rows = [
        [name, f"{value:.8e}", where or "n/a"]
        for name, (value, where) in max_disagreement.items()
    ]
    print(render_table(["category", "max_abs_diff", "location"], summary_rows))
    if args.absolute_debug:
        print("\nLargest absolute sigma_LL disagreements by category")
        absolute_summary_rows = [
            [name, f"{value:.8e}", where or "n/a"]
            for name, (value, where) in max_absolute_disagreement.items()
        ]
        print(render_table(["category", "max_abs_diff", "location"], absolute_summary_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
