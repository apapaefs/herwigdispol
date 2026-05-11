from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

from . import ew
from .born import DISPoint
from .pdfs import PDFPair, ratios_for_flavor


CF = 4.0 / 3.0
TR = 0.5
TWO_PI = 2.0 * math.pi
HERWIG_XP_POWER = 0.1


@dataclass(frozen=True)
class NLOTerms:
    virtual: float
    collinear_quark: float
    collinear_gluon: float
    real_qcdc: float
    real_bgf: float
    total: float


def herwig_xp_power_map(random_unit: float, x_b: float, power: float = HERWIG_XP_POWER) -> tuple[float, float]:
    if not (0.0 < x_b < 1.0):
        raise ValueError("Require 0 < x_b < 1 for the Herwig NLO xp map.")
    if not (0.0 < power < 1.0):
        raise ValueError("Require 0 < power < 1 for the Herwig NLO xp map.")
    r = min(max(random_unit, 1.0e-15), 1.0 - 1.0e-15)
    rhomin = (1.0 - x_b) ** (1.0 - power)
    x_p = 1.0 - (r * rhomin) ** (1.0 / (1.0 - power))
    jacobian = rhomin / (1.0 - power) * (1.0 - x_p) ** power
    return x_p, jacobian


def nlo_terms(
    point: DISPoint,
    x_p: float,
    pdfs: PDFPair,
    alpha_s: Callable[[float], float],
    jacobian: float = 1.0,
) -> NLOTerms:
    if not (0.0 < point.x_b < x_p < 1.0):
        raise ValueError("Require 0 < x_b < x_p < 1 for the Herwig NLO map.")
    mu2 = point.scale2
    a_s = alpha_s(mu2)
    cf_fact = CF * a_s / TWO_PI
    tr_fact = TR * a_s / TWO_PI

    x_b = point.x_b
    one_minus_xb = 1.0 - x_b
    virt = 1.0 + cf_fact * (
        -4.5
        - math.pi * math.pi / 3.0
        + 1.5 * math.log(point.q2 / mu2 / one_minus_xb)
        + 2.0 * math.log(one_minus_xb) * math.log(point.q2 / mu2)
        + math.log(one_minus_xb) ** 2
    )
    virt /= jacobian

    ratios = ratios_for_flavor(pdfs, point.flavor, point.x_b, x_p, mu2, point.hadron_pol)
    ell = point.ell
    a_born = ew.a_pol(point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, ratios.pq_born)
    blend = ew.collinear_blend_weights(
        point.setup,
        point.lepton_id,
        point.flavor,
        point.q2,
        point.lepton_pol,
        ratios.pq_born,
        ell,
    )
    response = ew.neutral_current_response(
        point.setup,
        point.lepton_id,
        point.flavor,
        point.q2,
        point.lepton_pol,
        ratios.pq_born,
        ell,
    )
    ratio_floor = 1.0e-12
    if response is not None:
        q_odd_weight = response.q_odd_response
        g_odd_weight = response.g_odd_response
    elif abs(ratios.pq_born) > ratio_floor:
        q_odd_weight = blend.q_polarized / ratios.pq_born
        g_odd_weight = blend.g_polarized / ratios.pq_born
    else:
        q_odd_weight = 0.0
        g_odd_weight = 0.0

    log_ratio = math.log((1.0 - x_p) * point.q2 / x_p / mu2)
    collg_over_born_unpol = tr_fact / x_p * ratios.g_ratio * (
        2.0 * x_p * (1.0 - x_p) + (x_p * x_p + (1.0 - x_p) ** 2) * log_ratio
    )
    collg_even = blend.g_unpolarized * collg_over_born_unpol
    collg_odd = (
        g_odd_weight
        * tr_fact
        / x_p
        * ratios.deltag_over_lo_eff
        * (2.0 * (1.0 - x_p) + (2.0 * x_p - 1.0) * log_ratio)
        if response is not None or ratios.has_stable_difference_ratio
        else 0.0
    )

    collq_k1 = 1.0 - x_p - 2.0 / (1.0 - x_p) * math.log(x_p) - (1.0 + x_p) * math.log(
        (1.0 - x_p) / x_p * point.q2 / mu2
    )
    collq_k2 = 2.0 / (1.0 - x_p) * math.log(point.q2 * (1.0 - x_p) / mu2) - 1.5 / (1.0 - x_p)
    collq_over_born_unpol = (
        cf_fact / x_p * ratios.q_ratio * collq_k1
        + cf_fact / x_p * (ratios.q_ratio - x_p) * collq_k2
    )
    collq_even = blend.q_unpolarized * collq_over_born_unpol
    collq_odd = (
        q_odd_weight
        * (
            cf_fact / x_p * ratios.deltaq_over_lo_eff * collq_k1
            + cf_fact / x_p * (ratios.deltaq_over_lo_eff - ratios.pq_born * x_p) * collq_k2
        )
        if ratios.has_stable_difference_ratio
        else 0.0
    )

    a_q_mapped = ew.a_pol(point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, ratios.pq_mapped)
    qcdc_den_ratio = ew.qcdc_mapped_denominator_ratio(
        point.setup,
        point.lepton_id,
        point.flavor,
        point.q2,
        point.lepton_pol,
        ratios.pq_born,
        ratios.pq_mapped,
    )
    born_factor = 1.0 + a_born * ell + ell * ell
    realq_prefactor = qcdc_den_ratio * cf_fact / x_p / born_factor * ratios.q_ratio
    realq_even = realq_prefactor * (2.0 + 2.0 * ell * ell - x_p + 3.0 * x_p * ell * ell)
    realq_odd = realq_prefactor * a_q_mapped * ell * (2.0 * x_p + 1.0)
    realq = realq_even + realq_odd

    realg_over_born_unpol = (
        -tr_fact
        / x_p
        * ratios.g_ratio
        * ((1.0 + ell * ell + 2.0 * (1.0 - 3.0 * ell * ell) * x_p * (1.0 - x_p)) / (1.0 + ell * ell))
    )
    realg_even = blend.g_unpolarized * realg_over_born_unpol
    realg_odd = (
        g_odd_weight * (-tr_fact / x_p * ratios.deltag_over_lo_eff * (2.0 * x_p - 1.0))
        if response is not None or ratios.has_stable_difference_ratio
        else 0.0
    )
    realg = realg_even + realg_odd

    collq = collq_even + collq_odd
    collg = collg_even + collg_odd
    total = virt + collq + collg + realq + realg
    return NLOTerms(virt, collq, collg, realq, realg, total)
