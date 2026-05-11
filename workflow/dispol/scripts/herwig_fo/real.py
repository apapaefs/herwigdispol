from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

from . import ew
from .born import DISPoint, lab_kinematics
from .config import BeamConfig
from .hepmc import EventParticle, HepMCEvent, ME_STATUS_ROLE2, ME_STATUS_ROLE3
from .nlo import CF, TR, TWO_PI
from .pdfs import PDFPair, PDFRatios, ratios_for_flavor
from .vectors import FourVector


@dataclass(frozen=True)
class PartonMomentum:
    pid: int
    p4: FourVector
    status: int
    role: int

    @property
    def px(self) -> float:
        return self.p4.px

    @property
    def py(self) -> float:
        return self.p4.py

    @property
    def pz(self) -> float:
        return self.p4.pz

    @property
    def e(self) -> float:
        return self.p4.e

    def pt(self) -> float:
        return self.p4.pt()


@dataclass(frozen=True)
class RealMomenta:
    incoming: PartonMomentum
    outgoing: tuple[PartonMomentum, PartonMomentum]


@dataclass(frozen=True)
class AzimuthalKernelCoefficients:
    c0: float
    c1: float
    c2: float

    def average(self) -> float:
        return self.c0 + 0.5 * self.c2

    def value(self, phi: float) -> float:
        cphi = math.cos(phi)
        return self.c0 + self.c1 * cphi + self.c2 * cphi * cphi


@dataclass(frozen=True)
class LocalRealWeightRatios:
    qcdc: float
    bgf: float


def _breit_two_body_momenta(q2: float, x_p: float, z_p: float, phi: float) -> tuple[FourVector, FourVector, FourVector]:
    if not (0.0 < x_p < 1.0):
        raise ValueError("x_p must be in (0, 1).")
    if not (0.0 < z_p < 1.0):
        raise ValueError("z_p must be in (0, 1).")
    q = math.sqrt(q2)
    x_perp = math.sqrt(max(0.0, 4.0 * (1.0 - x_p) * (1.0 - z_p) * z_p / x_p))
    x1 = -1.0 / x_p
    x2 = 1.0 - (1.0 - z_p) / x_p
    x3 = 2.0 + x1 - x2
    cphi = math.cos(phi)
    sphi = math.sin(phi)
    incoming = FourVector(0.0, 0.0, -0.5 * x1 * q, -0.5 * x1 * q)
    p1 = FourVector(0.5 * q * x_perp * cphi, 0.5 * q * x_perp * sphi, -0.5 * q * x2, 0.5 * q * math.hypot(x_perp, x2))
    p2 = FourVector(-p1.px, -p1.py, -0.5 * q * x3, 0.5 * q * math.hypot(x_perp, x3))
    return incoming, p1, p2


def breit_qcdc_momenta(q2: float, x_p: float, z_p: float, phi: float, flavor: int = 2) -> RealMomenta:
    incoming, p1, p2 = _breit_two_body_momenta(q2, x_p, z_p, phi)
    return RealMomenta(
        PartonMomentum(flavor, incoming, 990071, 1),
        (
            PartonMomentum(flavor, p1, ME_STATUS_ROLE2, 2),
            PartonMomentum(21, p2, ME_STATUS_ROLE3, 3),
        ),
    )


def breit_bgf_momenta(q2: float, x_p: float, z_p: float, phi: float, flavor: int = 2) -> RealMomenta:
    incoming, p1, p2 = _breit_two_body_momenta(q2, x_p, z_p, phi)
    return RealMomenta(
        PartonMomentum(21, incoming, 990071, 1),
        (
            PartonMomentum(flavor, p1, ME_STATUS_ROLE2, 2),
            PartonMomentum(-flavor, p2, ME_STATUS_ROLE3, 3),
        ),
    )


def _safe_den_sqrt(value: float) -> float:
    return math.sqrt(max(value, 1.0e-300))


def _is_safe_analyzing_power(value: float) -> bool:
    return math.isfinite(value) and abs(value) <= 2.0 + 1.0e-10


def _positive_ratio(numerator: float, denominator: float) -> tuple[bool, float]:
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        return False, 1.0
    floor = 1.0e-12 * max(1.0, abs(denominator))
    if abs(denominator) <= floor or numerator <= 0.0:
        return False, 1.0
    candidate = numerator / denominator
    if not math.isfinite(candidate) or candidate <= 0.0:
        return False, 1.0
    return True, candidate


def _real_variables(x_p: float, z_p: float) -> tuple[float, float, float, float]:
    if not (0.0 < x_p < 1.0):
        raise ValueError("x_p must be in (0, 1).")
    if not (0.0 < z_p < 1.0):
        raise ValueError("z_p must be in (0, 1).")
    x_perp2 = 4.0 * (1.0 - x_p) * (1.0 - z_p) * z_p / x_p
    x_perp = math.sqrt(max(0.0, x_perp2))
    x1 = -1.0 / x_p
    x2 = 1.0 - (1.0 - z_p) / x_p
    x3 = 2.0 + x1 - x2
    return x_perp, x1, x2, x3


def compton_azimuthal_coefficients(point: DISPoint, x_p: float, z_p: float, ratios: PDFRatios) -> AzimuthalKernelCoefficients:
    x_perp, _, x2, _ = _real_variables(x_p, z_p)
    den = _safe_den_sqrt(x2 * x2 + x_perp * x_perp)
    cos2 = x2 / den
    sin2 = x_perp / den
    ell = point.ell
    root = math.sqrt(max(0.0, ell * ell - 1.0))

    a_zero = ew.a_pol(point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, 0.0)
    a_born = ew.a_pol(point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, ratios.pq_born)
    if not _is_safe_analyzing_power(a_born):
        a_born = a_zero

    a_mapped = a_born
    mapped_den_ratio = 1.0
    d_born = ew.real_emission_denominator_factor(
        point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, ratios.pq_born
    )
    d_mapped = ew.real_emission_denominator_factor(
        point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, ratios.pq_mapped
    )
    a_candidate = ew.a_pol(
        point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, ratios.pq_mapped
    )
    ratio_ok, candidate_ratio = _positive_ratio(d_mapped, d_born)
    if _is_safe_analyzing_power(a_candidate) and ratio_ok:
        a_mapped = a_candidate
        mapped_den_ratio = candidate_ratio

    lo = 1.0 + a_born * ell + ell * ell
    if abs(lo) <= 1.0e-30:
        lo = math.copysign(1.0e-30, lo if lo != 0.0 else 1.0)
    mapped_lo = 1.0 + a_mapped * ell + ell * ell
    leading = mapped_den_ratio * mapped_lo / lo
    fact = mapped_den_ratio * x_p * x_p * (x2 * x2 + x_perp * x_perp) / lo
    return AzimuthalKernelCoefficients(
        leading + fact * (cos2 * cos2 + a_mapped * cos2 * ell + ell * ell),
        fact * (-a_mapped * cos2 * root * sin2 - 2.0 * ell * root * sin2),
        fact * (root * root * sin2 * sin2),
    )


def bgf_azimuthal_coefficients(point: DISPoint, x_p: float, z_p: float, ratios: PDFRatios) -> AzimuthalKernelCoefficients:
    x_perp, _, x2, x3 = _real_variables(x_p, z_p)
    den2 = _safe_den_sqrt(x2 * x2 + x_perp * x_perp)
    den3 = _safe_den_sqrt(x3 * x3 + x_perp * x_perp)
    cos2 = x2 / den2
    sin2 = x_perp / den2
    cos3 = x3 / den3
    sin3 = x_perp / den3
    fact2 = x_p * x_p * (x2 * x2 + x_perp * x_perp)
    fact3 = x_p * x_p * (x3 * x3 + x_perp * x_perp)
    ell = point.ell
    root = math.sqrt(max(0.0, ell * ell - 1.0))

    a_zero = ew.a_pol(point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, 0.0)
    a_born = ew.a_pol(point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, ratios.pq_born)
    if not _is_safe_analyzing_power(a_born):
        a_born = a_zero

    d_born = ew.real_emission_denominator_factor(
        point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, ratios.pq_born
    )

    a_r2 = a_born
    r2_den_ratio = 1.0
    d_r2 = ew.real_emission_denominator_factor(
        point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, ratios.pg_mapped
    )
    a_r2_candidate = ew.a_pol(
        point.setup, point.lepton_id, point.flavor, point.q2, point.lepton_pol, ratios.pg_mapped
    )
    r2_ok, r2_candidate_ratio = _positive_ratio(d_r2, d_born)
    if _is_safe_analyzing_power(a_r2_candidate) and r2_ok:
        a_r2 = a_r2_candidate
        r2_den_ratio = r2_candidate_ratio

    qbar_id = -point.flavor
    a_r3 = a_born
    r3_den_ratio = 1.0
    d_r3 = ew.real_emission_denominator_factor(
        point.setup, point.lepton_id, qbar_id, point.q2, point.lepton_pol, -ratios.pg_mapped
    )
    a_r3_candidate = ew.a_pol(
        point.setup, point.lepton_id, qbar_id, point.q2, point.lepton_pol, -ratios.pg_mapped
    )
    r3_ok, r3_candidate_ratio = _positive_ratio(d_r3, d_born)
    if _is_safe_analyzing_power(a_r3_candidate) and r3_ok:
        a_r3 = a_r3_candidate
        r3_den_ratio = r3_candidate_ratio

    lo = 1.0 + a_born * ell + ell * ell
    if abs(lo) <= 1.0e-30:
        lo = math.copysign(1.0e-30, lo if lo != 0.0 else 1.0)
    c0 = (
        fact2 * r2_den_ratio * (cos2 * cos2 + a_r2 * cos2 * ell + ell * ell)
        + fact3 * r3_den_ratio * (cos3 * cos3 - a_r3 * cos3 * ell + ell * ell)
    )
    c1 = (
        -fact2 * r2_den_ratio * (a_r2 * cos2 * root * sin2 + 2.0 * ell * root * sin2)
        - fact3 * r3_den_ratio * (a_r3 * cos3 * root * sin3 - 2.0 * ell * root * sin3)
    )
    c2 = (
        fact2 * r2_den_ratio * (root * root * sin2 * sin2)
        + fact3 * r3_den_ratio * (root * root * sin3 * sin3)
    )
    return AzimuthalKernelCoefficients(c0 / lo, c1 / lo, c2 / lo)


def qcdc_real_weight_ratio(
    point: DISPoint,
    x_p: float,
    z_p: float,
    phi: float,
    pdfs: PDFPair,
    alpha_s: Callable[[float], float],
) -> float:
    ratios = ratios_for_flavor(pdfs, point.flavor, point.x_b, x_p, point.scale2, point.hadron_pol)
    kernel = compton_azimuthal_coefficients(point, x_p, z_p, ratios).value(phi)
    return CF * alpha_s(point.scale2) / TWO_PI * ratios.q_ratio * kernel / (x_p * (1.0 - x_p) * (1.0 - z_p))


def bgf_real_weight_ratio(
    point: DISPoint,
    x_p: float,
    z_p: float,
    phi: float,
    pdfs: PDFPair,
    alpha_s: Callable[[float], float],
) -> float:
    ratios = ratios_for_flavor(pdfs, point.flavor, point.x_b, x_p, point.scale2, point.hadron_pol)
    kernel = bgf_azimuthal_coefficients(point, x_p, z_p, ratios).value(phi)
    # BGF is partitioned onto quark and antiquark underlying-Born branches.
    endpoint = (1.0 - z_p) if point.flavor > 0 else z_p
    return TR * alpha_s(point.scale2) / TWO_PI * ratios.g_ratio * kernel / (x_p * endpoint)


def local_real_weight_ratios(
    point: DISPoint,
    x_p: float,
    z_p: float,
    phi: float,
    pdfs: PDFPair,
    alpha_s: Callable[[float], float],
) -> LocalRealWeightRatios:
    return LocalRealWeightRatios(
        qcdc=qcdc_real_weight_ratio(point, x_p, z_p, phi, pdfs, alpha_s),
        bgf=bgf_real_weight_ratio(point, x_p, z_p, phi, pdfs, alpha_s),
    )


def _dot4(a: FourVector, b: FourVector) -> float:
    return a.e * b.e - a.px * b.px - a.py * b.py - a.pz * b.pz


def _project_orthogonal(seed: FourVector, basis: tuple[FourVector, ...]) -> FourVector:
    out = seed
    for axis in basis:
        axis2 = _dot4(axis, axis)
        if abs(axis2) <= 1.0e-14:
            continue
        out = out - axis.scale(_dot4(out, axis) / axis2)
    return out


def _normalize_timelike(v: FourVector) -> FourVector:
    norm2 = _dot4(v, v)
    if norm2 <= 0.0:
        raise ValueError("Expected a timelike Breit basis vector.")
    return v.scale(1.0 / math.sqrt(norm2))


def _normalize_spacelike(v: FourVector) -> FourVector:
    norm2 = -_dot4(v, v)
    if norm2 <= 0.0:
        raise ValueError("Expected a spacelike Breit basis vector.")
    return v.scale(1.0 / math.sqrt(norm2))


def _breit_basis_in_lab(point: DISPoint, beams: BeamConfig) -> tuple[FourVector, FourVector, FourVector, FourVector]:
    kin = lab_kinematics(point, beams)
    q = math.sqrt(point.q2)
    z_axis = _normalize_spacelike(kin.exchanged_boson.scale(-1.0 / q))
    t_axis = _normalize_timelike(kin.proton_in.scale(2.0 * point.x_b / q) - z_axis)

    x_seed = kin.electron_out
    x_axis = _project_orthogonal(x_seed, (t_axis, z_axis))
    if -_dot4(x_axis, x_axis) <= 1.0e-14:
        x_axis = _project_orthogonal(FourVector(1.0, 0.0, 0.0, 0.0), (t_axis, z_axis))
    x_axis = _normalize_spacelike(x_axis)

    y_axis = _project_orthogonal(FourVector(0.0, 1.0, 0.0, 0.0), (t_axis, z_axis, x_axis))
    if -_dot4(y_axis, y_axis) <= 1.0e-14:
        y_axis = _project_orthogonal(FourVector(0.0, 0.0, 1.0, 0.0), (t_axis, z_axis, x_axis))
    y_axis = _normalize_spacelike(y_axis)
    return t_axis, x_axis, y_axis, z_axis


def _breit_to_lab(point: DISPoint, p: FourVector, beams: BeamConfig) -> FourVector:
    t_axis, x_axis, y_axis, z_axis = _breit_basis_in_lab(point, beams)
    return (
        t_axis.scale(p.e)
        + x_axis.scale(p.px)
        + y_axis.scale(p.py)
        + z_axis.scale(p.pz)
    )


def real_event(
    point: DISPoint,
    channel: str,
    x_p: float,
    z_p: float,
    phi: float,
    weight: float,
    event_number: int,
    beams: BeamConfig = BeamConfig(),
) -> HepMCEvent:
    channel_key = channel.upper()
    if channel_key == "QCDC":
        momenta = breit_qcdc_momenta(point.q2, x_p, z_p, phi, point.flavor)
    elif channel_key == "BGF":
        momenta = breit_bgf_momenta(point.q2, x_p, z_p, phi, point.flavor)
    else:
        raise ValueError("channel must be QCDC or BGF.")

    kin = lab_kinematics(point, beams)
    incoming_parton_lab = _breit_to_lab(point, momenta.incoming.p4, beams)
    outgoing = [
        EventParticle(point.lepton_id, kin.electron_out, 1),
    ]
    for parton in momenta.outgoing:
        outgoing.append(
            EventParticle(
                parton.pid,
                _breit_to_lab(point, parton.p4, beams),
                parton.status,
                role=parton.role,
            )
        )
    outgoing.append(EventParticle(82, kin.proton_in - incoming_parton_lab, 1))

    return HepMCEvent(
        event_number=event_number,
        weight=weight,
        incoming=(
            EventParticle(beams.hadron_id, kin.proton_in, 4),
            EventParticle(point.lepton_id, kin.electron_in, 4),
        ),
        outgoing=tuple(outgoing),
        attributes={
            "HERWIG_FO_EVENT_CLASS": channel_key,
            "HERWIG_FO_SETUP": point.setup,
        },
    )
