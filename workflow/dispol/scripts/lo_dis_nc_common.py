#!/usr/bin/env python3.10
"""Common LO neutral-current DIS helpers for POLDIS-style reference checks."""

from __future__ import annotations

import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence

import numpy as np


DEFAULT_LHAPDF_DATA_PATH = "/opt/homebrew/Cellar/lhapdf/6.5.4/share/LHAPDF"
DEFAULT_SUM_PDF = "PDF4LHC15_nnlo_100_pdfas"
DEFAULT_DIFF_PDF = "BDSSV24-NNLO"
DEFAULT_BEAM_EA_GEV = 18.0
DEFAULT_BEAM_EB_GEV = 275.0
DEFAULT_Q2_POINT = 100.0
DEFAULT_Y_POINT = 0.4
DEFAULT_Q2_MIN = 49.0
DEFAULT_Q2_MAX = 2500.0
DEFAULT_Y_MIN = 0.2
DEFAULT_Y_MAX = 0.6
DEFAULT_ALPHA_EM = 0.00729927
DEFAULT_SIN2_THETA_W = 0.221639970483740179
DEFAULT_MW = 80.450
DEFAULT_MZ = 91.1876
DEFAULT_GAMMA_Z = 2.4952
DEFAULT_LEPTON_CHARGE = -1.0
POLDIS_NORM_PB = 389385.0 * 1000.0
THEPEG_HBARC_SQ_PB = (197.326968e-18 ** 2) / 1.0e-40

CHANNEL_GAMMA = "GAMMA"
CHANNEL_Z = "Z"
CHANNEL_INT = "INT"
CHANNEL_ALL = "ALL"
CHANNEL_ORDER = (CHANNEL_GAMMA, CHANNEL_Z, CHANNEL_INT, CHANNEL_ALL)
HERWIG_BASE_QUARK_MASSES_GEV: Dict[str, float] = {
    "u": 0.00216,
    "d": 0.00467,
    "s": 0.0934,
    "c": 1.27,
    "b": 4.18,
}


@dataclass(frozen=True)
class EWInputs:
    sin2_theta_w: float
    alpha_em: float
    mw: float
    mz: float
    gamma_z: float
    lepton_charge: float


@dataclass(frozen=True)
class LOCalculatorMode:
    name: str
    spacelike_z_width: bool
    enforce_ep_x_cuts: bool
    x1_min: float
    x2_min: float
    w2_min: float
    w2_max: float


@dataclass(frozen=True)
class QuarkMassModel:
    name: str
    rescale_pdf_x: bool
    apply_jacobian: bool
    use_massive_born_kernel: bool = False


DEFAULT_MODE_POLDIS = LOCalculatorMode(
    name="poldis",
    spacelike_z_width=True,
    enforce_ep_x_cuts=False,
    x1_min=0.0,
    x2_min=0.0,
    w2_min=0.0,
    w2_max=float("inf"),
)

DEFAULT_MODE_HERWIG_LO = LOCalculatorMode(
    name="herwig_lo",
    spacelike_z_width=False,
    enforce_ep_x_cuts=True,
    x1_min=1.0e-5,
    x2_min=1.0e-5,
    w2_min=0.0,
    w2_max=float("inf"),
)

DEFAULT_QUARK_MASS_MODEL = QuarkMassModel(
    name="massless",
    rescale_pdf_x=False,
    apply_jacobian=False,
)
DEFAULT_QUARK_MASS_MODEL_SLOWRESCALE = QuarkMassModel(
    name="slowrescale",
    rescale_pdf_x=True,
    apply_jacobian=False,
)
DEFAULT_QUARK_MASS_MODEL_SLOWRESCALE_JACOBIAN = QuarkMassModel(
    name="slowrescale_jacobian",
    rescale_pdf_x=True,
    apply_jacobian=True,
)
DEFAULT_QUARK_MASS_MODEL_MASSIVE_BORN = QuarkMassModel(
    name="massive_born",
    rescale_pdf_x=False,
    apply_jacobian=False,
    use_massive_born_kernel=True,
)


@dataclass(frozen=True)
class BeamInputs:
    beam_ea_gev: float
    beam_eb_gev: float

    @property
    def s_hadronic(self) -> float:
        return 4.0 * self.beam_ea_gev * self.beam_eb_gev


@dataclass(frozen=True)
class Flavor:
    name: str
    pdg_id: int
    poldis_index: int
    species_charge: float
    is_up_type: bool


@dataclass(frozen=True)
class ChannelPieces:
    vec_gamma: float
    vec_interference: float
    vec_z: float
    ax_interference: float
    ax_z: float


@dataclass(frozen=True)
class ChannelKernel:
    vec: float
    ax: float
    sigma_unpol_kernel: float
    sigma_ll_kernel: float


@dataclass(frozen=True)
class PointContribution:
    flavor: str
    channel: str
    x_b: float
    x_pdf: float
    mass_gev: float
    mass_weight: float
    xf_unpolarized: float
    xf_polarized: float
    vec: float
    ax: float
    sigma_unpol_kernel: float
    sigma_ll_kernel: float
    structure_unpol: float
    structure_ll: float
    d2sigma_unpol_pb_per_gev2: float
    d2sigma_ll_pb_per_gev2: float


@dataclass(frozen=True)
class IntegratedContribution:
    flavor: str
    channel: str
    sigma_unpol_pb: float
    sigma_ll_pb: float


FLAVORS: List[Flavor] = [
    Flavor("u", 2, 2, 2.0 / 3.0, True),
    Flavor("d", 1, 1, -1.0 / 3.0, False),
    Flavor("s", 3, 3, -1.0 / 3.0, False),
    Flavor("c", 4, 4, 2.0 / 3.0, True),
    Flavor("b", 5, 5, -1.0 / 3.0, False),
    Flavor("ubar", -2, -2, 2.0 / 3.0, True),
    Flavor("dbar", -1, -1, -1.0 / 3.0, False),
    Flavor("sbar", -3, -3, -1.0 / 3.0, False),
    Flavor("cbar", -4, -4, 2.0 / 3.0, True),
    Flavor("bbar", -5, -5, -1.0 / 3.0, False),
]
FLAVOR_LOOKUP: Dict[str, Flavor] = {flavor.name: flavor for flavor in FLAVORS}
PDG_FLAVOR_LOOKUP: Dict[int, Flavor] = {flavor.pdg_id: flavor for flavor in FLAVORS}
HELICITY_COMBINATIONS: tuple[tuple[int, int], ...] = ((+1, +1), (+1, -1), (-1, +1), (-1, -1))

IDENTITY4 = np.eye(4, dtype=np.complex128)
GAMMA0 = np.array(
    [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, -1.0, 0.0], [0.0, 0.0, 0.0, -1.0]],
    dtype=np.complex128,
)
GAMMA1 = np.array(
    [[0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0, 0.0], [0.0, -1.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0]],
    dtype=np.complex128,
)
GAMMA2 = np.array(
    [[0.0, 0.0, 0.0, -1.0j], [0.0, 0.0, 1.0j, 0.0], [0.0, 1.0j, 0.0, 0.0], [-1.0j, 0.0, 0.0, 0.0]],
    dtype=np.complex128,
)
GAMMA3 = np.array(
    [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, -1.0], [-1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]],
    dtype=np.complex128,
)
GAMMA5 = 1.0j * GAMMA0 @ GAMMA1 @ GAMMA2 @ GAMMA3
GAMMA_UP = (GAMMA0, GAMMA1, GAMMA2, GAMMA3)
GAMMA_DOWN = (GAMMA0, -GAMMA1, -GAMMA2, -GAMMA3)


def herwig_default_quark_mass_map() -> Dict[str, float]:
    masses: Dict[str, float] = {}
    for name, mass in HERWIG_BASE_QUARK_MASSES_GEV.items():
        masses[name] = mass
        masses[f"{name}bar"] = mass
    return masses


def _slash(momentum: Sequence[float]) -> np.ndarray:
    return sum(float(momentum[mu]) * GAMMA_DOWN[mu] for mu in range(4))


def _massless_incoming_density(momentum: Sequence[float], helicity: int) -> np.ndarray:
    return 0.5 * _slash(momentum) @ (IDENTITY4 + float(helicity) * GAMMA5)


def _massive_spin_vector(momentum: Sequence[float], mass_gev: float) -> np.ndarray:
    px, py, pz = float(momentum[1]), float(momentum[2]), float(momentum[3])
    p_abs = math.sqrt(px * px + py * py + pz * pz)
    if p_abs <= 0.0:
        return np.array([0.0, 0.0, 0.0, 1.0], dtype=float)
    energy = float(momentum[0])
    return np.array(
        [p_abs / mass_gev, energy * px / (mass_gev * p_abs), energy * py / (mass_gev * p_abs), energy * pz / (mass_gev * p_abs)],
        dtype=float,
    )


def _massive_incoming_density(momentum: Sequence[float], mass_gev: float, helicity: int) -> np.ndarray:
    spin_vec = _massive_spin_vector(momentum, mass_gev)
    return 0.5 * (_slash(momentum) + mass_gev * IDENTITY4) @ (IDENTITY4 - float(helicity) * GAMMA5 @ _slash(spin_vec))


def _fermion_sum_density(momentum: Sequence[float], mass_gev: float) -> np.ndarray:
    return _slash(momentum) + mass_gev * IDENTITY4


def _vertex_up(gv: float, ga: float, mu: int) -> np.ndarray:
    return GAMMA_UP[mu] @ (gv * IDENTITY4 + ga * GAMMA5)


def _vertex_down(gv: float, ga: float, mu: int) -> np.ndarray:
    return GAMMA_DOWN[mu] @ (gv * IDENTITY4 + ga * GAMMA5)


def _partonic_kinematics(q2: float, y: float, mass_gev: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    shat = q2 / y
    sqrt_shat = math.sqrt(shat)
    p_abs = (shat - mass_gev * mass_gev) / (2.0 * sqrt_shat)
    if p_abs <= 0.0:
        raise ValueError("Partonic COM momentum is non-positive.")
    cos_theta = 1.0 - q2 / (2.0 * p_abs * p_abs)
    if cos_theta < -1.0 - 1.0e-10 or cos_theta > 1.0 + 1.0e-10:
        raise ValueError("Requested (Q^2, y, m_q) lies outside the 2->2 kinematic range.")
    cos_theta = max(-1.0, min(1.0, cos_theta))
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    lepton_energy = p_abs
    quark_energy = (shat + mass_gev * mass_gev) / (2.0 * sqrt_shat)
    k_in = np.array([lepton_energy, 0.0, 0.0, lepton_energy], dtype=float)
    p_in = np.array([quark_energy, 0.0, 0.0, -p_abs], dtype=float)
    k_out = np.array([lepton_energy, lepton_energy * sin_theta, 0.0, lepton_energy * cos_theta], dtype=float)
    p_out = np.array([quark_energy, -lepton_energy * sin_theta, 0.0, -lepton_energy * cos_theta], dtype=float)
    return k_in, p_in, k_out, p_out


def _numerical_boson_data(
    flavor: Flavor, q2: float, inputs: EWInputs, mode: LOCalculatorMode
) -> dict[str, dict[str, object]]:
    eq, cvz, caz = poldis_quark_couplings(flavor, inputs)
    cvze, caze = poldis_lepton_couplings(inputs)
    width = inputs.mz * inputs.gamma_z if mode.spacelike_z_width else 0.0
    z_coeff = inputs.lepton_charge * q2 / complex(q2 + inputs.mz * inputs.mz, width)
    return {
        CHANNEL_GAMMA: {
            "coeff": complex(1.0, 0.0),
            "lepton": (1.0, 0.0),
            "quark": (eq, 0.0),
        },
        CHANNEL_Z: {
            "coeff": z_coeff,
            "lepton": (cvze, caze),
            "quark": (cvz, caz),
        },
    }


def _channel_helicity_weight(
    flavor: Flavor,
    q2: float,
    y: float,
    inputs: EWInputs,
    mass_gev: float,
    lepton_helicity: int,
    quark_helicity: int,
    mode: LOCalculatorMode,
    channel: str,
) -> float:
    k_in, p_in, k_out, p_out = _partonic_kinematics(q2, y, mass_gev)
    rho_l_in = _massless_incoming_density(k_in, lepton_helicity)
    rho_l_out = _fermion_sum_density(k_out, 0.0)
    if mass_gev > 0.0:
        rho_q_in = _massive_incoming_density(p_in, mass_gev, quark_helicity)
    else:
        rho_q_in = _massless_incoming_density(p_in, quark_helicity)
    rho_q_out = _fermion_sum_density(p_out, mass_gev)

    bosons = _numerical_boson_data(flavor, q2, inputs, mode)
    boson_terms = ("GAMMA", "Z") if channel == CHANNEL_ALL else (channel,)
    total = 0.0 + 0.0j
    for left in boson_terms:
        for right in boson_terms:
            coeff = bosons[left]["coeff"] * complex(bosons[right]["coeff"]).conjugate()
            lepton_left = bosons[left]["lepton"]
            lepton_right = bosons[right]["lepton"]
            quark_left = bosons[left]["quark"]
            quark_right = bosons[right]["quark"]
            block = 0.0 + 0.0j
            for mu in range(4):
                v_l_up = _vertex_up(float(lepton_left[0]), float(lepton_left[1]), mu)
                v_q_down = _vertex_down(float(quark_left[0]), float(quark_left[1]), mu)
                for nu in range(4):
                    v_l_up_r = _vertex_up(float(lepton_right[0]), float(lepton_right[1]), nu)
                    v_q_down_r = _vertex_down(float(quark_right[0]), float(quark_right[1]), nu)
                    l_trace = np.trace(rho_l_in @ v_l_up @ rho_l_out @ v_l_up_r)
                    q_trace = np.trace(rho_q_in @ v_q_down @ rho_q_out @ v_q_down_r)
                    block += l_trace * q_trace
            total += coeff * block
    return float(total.real)


def massive_born_kernel_ratios(
    flavor: Flavor,
    q2: float,
    y: float,
    inputs: EWInputs,
    mass_gev: float,
    mode: LOCalculatorMode,
) -> dict[str, tuple[float, float]]:
    if mass_gev <= 0.0:
        return {channel: (1.0, 1.0) for channel in CHANNEL_ORDER}

    ratios: dict[str, tuple[float, float]] = {}
    for channel in (CHANNEL_GAMMA, CHANNEL_Z, CHANNEL_ALL):
        weights_massive: list[float] = []
        weights_massless: list[float] = []
        for lepton_helicity, quark_helicity in HELICITY_COMBINATIONS:
            weights_massive.append(
                _channel_helicity_weight(
                    flavor, q2, y, inputs, mass_gev, lepton_helicity, quark_helicity, mode, channel
                )
            )
            weights_massless.append(
                _channel_helicity_weight(
                    flavor, q2, y, inputs, 0.0, lepton_helicity, quark_helicity, mode, channel
                )
            )

        sigma0_massive = 0.25 * sum(weights_massive)
        sigma0_massless = 0.25 * sum(weights_massless)
        sigma_ll_massive = 0.25 * (
            weights_massive[0] + weights_massive[3] - weights_massive[1] - weights_massive[2]
        )
        sigma_ll_massless = 0.25 * (
            weights_massless[0] + weights_massless[3] - weights_massless[1] - weights_massless[2]
        )
        r0 = sigma0_massive / sigma0_massless if abs(sigma0_massless) > 1.0e-30 else 1.0
        rll = sigma_ll_massive / sigma_ll_massless if abs(sigma_ll_massless) > 1.0e-30 else 1.0
        ratios[channel] = (r0, rll)

    all_r0, all_rll = ratios[CHANNEL_ALL]
    gamma_r0, gamma_rll = ratios[CHANNEL_GAMMA]
    z_r0, z_rll = ratios[CHANNEL_Z]
    ratios[CHANNEL_INT] = (
        all_r0 if abs(all_r0) > 0.0 else 1.0,
        all_rll if abs(all_rll) > 0.0 else 1.0,
    )
    return ratios


def resolve_sin2_theta_w(theta: float | None, sin2_theta_w: float | None) -> float:
    if theta is None and sin2_theta_w is None:
        return DEFAULT_SIN2_THETA_W
    if theta is None:
        assert sin2_theta_w is not None
        return sin2_theta_w
    derived = math.sin(theta) ** 2
    if sin2_theta_w is not None and abs(derived - sin2_theta_w) > 1.0e-12:
        raise SystemExit(
            f"Inconsistent electroweak inputs: theta implies sin^2(theta_W)={derived:.15f} "
            f"but --sin2-theta-w={sin2_theta_w:.15f}."
        )
    return derived


def theta_w(inputs: EWInputs) -> float:
    return math.asin(math.sqrt(inputs.sin2_theta_w))


def poldis_lepton_couplings(inputs: EWInputs) -> tuple[float, float]:
    theta = theta_w(inputs)
    cvze = (-0.5 / math.sin(2.0 * theta) + math.tan(theta)) * (-inputs.lepton_charge)
    caze = -0.5 / math.sin(2.0 * theta)
    return cvze, caze


def poldis_quark_couplings(flavor: Flavor, inputs: EWInputs) -> tuple[float, float, float]:
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


def propagator_factors(
    q2: float, inputs: EWInputs, mode: LOCalculatorMode = DEFAULT_MODE_POLDIS
) -> tuple[float, float]:
    width_term = (inputs.mz * inputs.gamma_z) ** 2 if mode.spacelike_z_width else 0.0
    denominator = (q2 + inputs.mz * inputs.mz) ** 2 + width_term
    propz = q2 * q2 / denominator
    intgz = inputs.lepton_charge * (q2 + inputs.mz * inputs.mz) / q2
    return propz, intgz


def y_plus(y: float) -> float:
    return 1.0 + (1.0 - y) ** 2


def y_minus(y: float) -> float:
    return 1.0 - (1.0 - y) ** 2


def x_bjorken(q2: float, y: float, s_hadronic: float) -> float:
    return q2 / (y * s_hadronic)


def differential_prefactor_pb(q2: float, y: float, alpha_em: float) -> float:
    return 2.0 * math.pi * alpha_em * alpha_em * POLDIS_NORM_PB / (y * q2 * q2)


def channel_pieces(
    flavor: Flavor,
    q2: float,
    inputs: EWInputs,
    mode: LOCalculatorMode = DEFAULT_MODE_POLDIS,
) -> ChannelPieces:
    eq, cvz, caz = poldis_quark_couplings(flavor, inputs)
    cvze, caze = poldis_lepton_couplings(inputs)
    propz, intgz = propagator_factors(q2, inputs, mode)
    return ChannelPieces(
        vec_gamma=eq * eq,
        vec_interference=2.0 * eq * propz * intgz * (cvz * cvze),
        vec_z=propz * (cvz * cvz + caz * caz) * (cvze * cvze + caze * caze),
        ax_interference=2.0 * eq * propz * intgz * (caz * caze),
        ax_z=propz * (4.0 * cvz * caz * cvze * caze),
    )


def kernel_from_pieces(pieces: ChannelPieces, y: float) -> Dict[str, ChannelKernel]:
    yp = y_plus(y)
    ym = y_minus(y)

    def build(vec: float, ax: float) -> ChannelKernel:
        return ChannelKernel(
            vec=vec,
            ax=ax,
            sigma_unpol_kernel=yp * vec + ym * ax,
            sigma_ll_kernel=ym * vec + yp * ax,
        )

    gamma = build(pieces.vec_gamma, 0.0)
    z_only = build(pieces.vec_z, pieces.ax_z)
    interference = build(pieces.vec_interference, pieces.ax_interference)
    all_nc = build(
        pieces.vec_gamma + pieces.vec_interference + pieces.vec_z,
        pieces.ax_interference + pieces.ax_z,
    )
    return {
        CHANNEL_GAMMA: gamma,
        CHANNEL_Z: z_only,
        CHANNEL_INT: interference,
        CHANNEL_ALL: all_nc,
    }


def resolve_channels(channel: str) -> List[str]:
    upper = channel.upper()
    if upper == "ALL":
        return [CHANNEL_GAMMA, CHANNEL_Z, CHANNEL_INT, CHANNEL_ALL]
    if upper in CHANNEL_ORDER:
        return [upper]
    raise SystemExit(f"Unknown channel selection '{channel}'.")


def import_lhapdf():
    if "LHAPDF_DATA_PATH" not in os.environ:
        default_path = Path(DEFAULT_LHAPDF_DATA_PATH)
        if default_path.exists():
            os.environ["LHAPDF_DATA_PATH"] = str(default_path)
    try:
        import lhapdf  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "Could not import the LHAPDF Python bindings. Run this tool with a Python "
            "environment that has `lhapdf` available, such as python3.10 in this workspace."
        ) from exc
    if hasattr(lhapdf, "setVerbosity"):
        lhapdf.setVerbosity(0)
    return lhapdf


def pdf_xf(pdf, pid: int, x: float, q2: float) -> float:
    return float(pdf.xfxQ2(pid, x, q2))


def pdf_density(pdf, pid: int, x: float, q2: float) -> float:
    xf = pdf_xf(pdf, pid, x, q2)
    return xf / x


def lookup_flavor_by_pdg_id(pid: int) -> Flavor:
    try:
        return PDG_FLAVOR_LOOKUP[int(pid)]
    except KeyError as exc:
        raise SystemExit(f"Unsupported DIS flavor PDG id for LO audit: {pid}") from exc


def born_me2_from_unpolarized_kernel(y: float, alpha_em: float, sigma_unpol_kernel: float) -> float:
    if y <= 0.0:
        raise SystemExit("Require y > 0 to reconstruct the Born matrix element.")
    return 32.0 * math.pi * math.pi * alpha_em * alpha_em * sigma_unpol_kernel / (y * y)


def sigma_hat_pb_from_unpolarized_kernel(
    y: float,
    alpha_em: float,
    sigma_unpol_kernel: float,
    jacobian: float,
    s_hat: float,
    conversion_pb: float = POLDIS_NORM_PB,
) -> float:
    if s_hat <= 0.0:
        raise SystemExit("Require s_hat > 0 to reconstruct sigma_hat.")
    me2 = born_me2_from_unpolarized_kernel(y, alpha_em, sigma_unpol_kernel)
    return me2 * jacobian * conversion_pb / (16.0 * math.pi * math.pi * s_hat)


def plain_lo_gamma_hwmebase_jacobian(
    q2: float,
    s_hat: float,
    q2_min_cut: float,
) -> Dict[str, float | str]:
    if q2 <= 0.0:
        raise SystemExit("Require q2 > 0 to reconstruct the LO Gamma Jacobian.")
    if s_hat <= 0.0:
        raise SystemExit("Require s_hat > 0 to reconstruct the LO Gamma Jacobian.")
    if q2 >= s_hat:
        raise SystemExit("Require q2 < s_hat for the LO Gamma Jacobian audit.")

    cth = 1.0 - 2.0 * q2 / s_hat
    if q2_min_cut <= 0.0:
        return {
            "branch": "flat",
            "ctmin": -1.0,
            "ctmax": 1.0,
            "cth": cth,
            "jacobian": math.pi,
        }

    ctmax = 1.0 - 2.0 * q2_min_cut / s_hat
    if ctmax <= -1.0:
        raise SystemExit(
            f"MinQ2={q2_min_cut:.8f} GeV^2 is too large for sHat={s_hat:.8f} GeV^2."
        )
    if cth > ctmax + 1.0e-12 or cth < -1.0 - 1.0e-12:
        raise SystemExit(
            f"Point with q2={q2:.8f} GeV^2 and sHat={s_hat:.8f} GeV^2 lies outside the "
            f"plain LO DIS angular map: cth={cth:.16e}, ctmax={ctmax:.16e}."
        )

    jacobian = math.pi * (q2 / s_hat) * math.log(s_hat / q2_min_cut)
    return {
        "branch": "one_sided_q2min",
        "ctmin": -1.0,
        "ctmax": ctmax,
        "cth": cth,
        "jacobian": jacobian,
    }


def point_audit_observables(
    q2: float,
    y: float,
    qid: int,
    ew_inputs: EWInputs,
    beams: BeamInputs,
    sum_pdf,
    diff_pdf,
    channel: str = CHANNEL_GAMMA,
    mode: LOCalculatorMode = DEFAULT_MODE_POLDIS,
    pdf_q2: float | None = None,
    jacobian: float | None = None,
    s_hat: float | None = None,
    x_b: float | None = None,
) -> Dict[str, float | str]:
    flavor = lookup_flavor_by_pdg_id(qid)
    rows = point_contributions(
        q2,
        y,
        ew_inputs,
        beams,
        sum_pdf,
        diff_pdf,
        flavors=[flavor],
        mode=mode,
        quark_masses_gev={},
        mass_model=DEFAULT_QUARK_MASS_MODEL,
    )
    matching = [row for row in rows if row.channel == channel]
    if len(matching) != 1:
        raise SystemExit(
            f"Expected exactly one point contribution for flavor {flavor.name} and channel {channel}, found {len(matching)}."
        )
    row = matching[0]
    x_point = row.x_b if x_b is None else x_b
    if x_point <= 0.0 or x_point >= 1.0:
        raise SystemExit(f"Audit x_B={x_point:.8f} lies outside the physical interval (0,1).")
    pdf_scale_q2 = q2 if pdf_q2 is None else pdf_q2
    pdf_sum = pdf_density(sum_pdf, flavor.pdg_id, x_point, pdf_scale_q2)
    hard_y = q2 / s_hat if s_hat is not None and s_hat > 0.0 else y
    pieces = channel_pieces(flavor, q2, ew_inputs, mode)
    hard_kernel = kernel_from_pieces(pieces, hard_y)[channel].sigma_unpol_kernel
    result: Dict[str, float | str] = {
        "flavor": flavor.name,
        "channel": channel,
        "qid": int(qid),
        "x_b": x_point,
        "x_pdf": x_point,
        "pdf_q2": pdf_scale_q2,
        "pdf_sum": pdf_sum,
        "sigma_unpol_kernel": hard_kernel,
        "d2sigma_unpol_pb_per_gev2": differential_prefactor_pb(q2, y, ew_inputs.alpha_em)
        * pdf_xf(sum_pdf, flavor.pdg_id, x_point, pdf_scale_q2)
        * row.sigma_unpol_kernel,
    }
    if jacobian is not None and s_hat is not None:
        sigma_hat_pb = sigma_hat_pb_from_unpolarized_kernel(
            hard_y,
            ew_inputs.alpha_em,
            hard_kernel,
            jacobian,
            s_hat,
            conversion_pb=THEPEG_HBARC_SQ_PB,
        )
        result["sigma_hat_pb"] = sigma_hat_pb
        result["sigma_hat_nb"] = sigma_hat_pb / 1000.0
    return result


def flavor_pdf_x_and_weight(
    flavor: Flavor,
    x_b: float,
    q2: float,
    quark_masses_gev: Mapping[str, float],
    mass_model: QuarkMassModel,
) -> tuple[float, float, float]:
    mass_gev = float(quark_masses_gev.get(flavor.name, 0.0))
    if mass_gev <= 0.0 or not mass_model.rescale_pdf_x:
        return mass_gev, x_b, 1.0
    mass2 = mass_gev * mass_gev
    x_pdf = x_b * (1.0 + mass2 / q2)
    weight = q2 / (q2 + mass2) if mass_model.apply_jacobian else 1.0
    return mass_gev, x_pdf, weight


def point_contributions(
    q2: float,
    y: float,
    inputs: EWInputs,
    beams: BeamInputs,
    sum_pdf,
    diff_pdf,
    flavors: Sequence[Flavor] | None = None,
    mode: LOCalculatorMode = DEFAULT_MODE_POLDIS,
    quark_masses_gev: Mapping[str, float] | None = None,
    mass_model: QuarkMassModel = DEFAULT_QUARK_MASS_MODEL,
) -> List[PointContribution]:
    if q2 <= 0.0:
        raise SystemExit("Require q2 > 0.")
    if y <= 0.0 or y >= 1.0:
        raise SystemExit("Require 0 < y < 1.")

    x_b = x_bjorken(q2, y, beams.s_hadronic)
    if x_b <= 0.0 or x_b >= 1.0:
        raise SystemExit(
            f"x_B={x_b:.8f} lies outside the physical interval (0,1) for the chosen point."
        )
    if not point_passes_lo_mode_cuts(q2, x_b, mode):
        raise SystemExit(
            f"The chosen point fails the active LO convention cuts: x_B={x_b:.8f}, "
            f"W2={w2_from_x_q2(x_b, q2):.8f} GeV^2."
        )

    prefactor = differential_prefactor_pb(q2, y, inputs.alpha_em)
    rows: List[PointContribution] = []
    masses = quark_masses_gev or {}

    for flavor in list(flavors) if flavors is not None else FLAVORS:
        pieces = channel_pieces(flavor, q2, inputs, mode)
        kernels = kernel_from_pieces(pieces, y)
        mass_gev, x_pdf, mass_weight = flavor_pdf_x_and_weight(
            flavor, x_b, q2, masses, mass_model
        )
        kernel_ratios = (
            massive_born_kernel_ratios(flavor, q2, y, inputs, mass_gev, mode)
            if mass_model.use_massive_born_kernel
            else {channel: (1.0, 1.0) for channel in CHANNEL_ORDER}
        )
        if x_pdf <= 0.0 or x_pdf >= 1.0:
            continue
        xf_unpol = pdf_xf(sum_pdf, flavor.pdg_id, x_pdf, q2)
        xf_pol = pdf_xf(diff_pdf, flavor.pdg_id, x_pdf, q2)
        sigma_unpol_kernels = {
            CHANNEL_GAMMA: kernels[CHANNEL_GAMMA].sigma_unpol_kernel * kernel_ratios[CHANNEL_GAMMA][0],
            CHANNEL_Z: kernels[CHANNEL_Z].sigma_unpol_kernel * kernel_ratios[CHANNEL_Z][0],
            CHANNEL_ALL: kernels[CHANNEL_ALL].sigma_unpol_kernel * kernel_ratios[CHANNEL_ALL][0],
        }
        sigma_ll_kernels = {
            CHANNEL_GAMMA: kernels[CHANNEL_GAMMA].sigma_ll_kernel * kernel_ratios[CHANNEL_GAMMA][1],
            CHANNEL_Z: kernels[CHANNEL_Z].sigma_ll_kernel * kernel_ratios[CHANNEL_Z][1],
            CHANNEL_ALL: kernels[CHANNEL_ALL].sigma_ll_kernel * kernel_ratios[CHANNEL_ALL][1],
        }
        sigma_unpol_kernels[CHANNEL_INT] = (
            sigma_unpol_kernels[CHANNEL_ALL]
            - sigma_unpol_kernels[CHANNEL_GAMMA]
            - sigma_unpol_kernels[CHANNEL_Z]
        )
        sigma_ll_kernels[CHANNEL_INT] = (
            sigma_ll_kernels[CHANNEL_ALL]
            - sigma_ll_kernels[CHANNEL_GAMMA]
            - sigma_ll_kernels[CHANNEL_Z]
        )
        for channel in CHANNEL_ORDER:
            kernel = kernels[channel]
            structure_unpol = mass_weight * xf_unpol * sigma_unpol_kernels[channel]
            structure_ll = mass_weight * xf_pol * sigma_ll_kernels[channel]
            rows.append(
                PointContribution(
                    flavor=flavor.name,
                    channel=channel,
                    x_b=x_b,
                    x_pdf=x_pdf,
                    mass_gev=mass_gev,
                    mass_weight=mass_weight,
                    xf_unpolarized=xf_unpol,
                    xf_polarized=xf_pol,
                    vec=kernel.vec,
                    ax=kernel.ax,
                    sigma_unpol_kernel=sigma_unpol_kernels[channel],
                    sigma_ll_kernel=sigma_ll_kernels[channel],
                    structure_unpol=structure_unpol,
                    structure_ll=structure_ll,
                    d2sigma_unpol_pb_per_gev2=prefactor * structure_unpol,
                    d2sigma_ll_pb_per_gev2=prefactor * structure_ll,
                )
            )
    return rows


def sum_point_contributions(rows: Iterable[PointContribution]) -> Dict[str, PointContribution]:
    totals: Dict[str, PointContribution] = {}
    grouped: Dict[str, List[PointContribution]] = {}
    for row in rows:
        grouped.setdefault(row.channel, []).append(row)

    for channel, members in grouped.items():
        prototype = members[0]
        totals[channel] = PointContribution(
            flavor="TOTAL",
            channel=channel,
            x_b=prototype.x_b,
            x_pdf=prototype.x_b,
            mass_gev=0.0,
            mass_weight=1.0,
            xf_unpolarized=sum(item.xf_unpolarized for item in members),
            xf_polarized=sum(item.xf_polarized for item in members),
            vec=sum(item.vec for item in members),
            ax=sum(item.ax for item in members),
            sigma_unpol_kernel=sum(item.structure_unpol for item in members),
            sigma_ll_kernel=sum(item.structure_ll for item in members),
            structure_unpol=sum(item.structure_unpol for item in members),
            structure_ll=sum(item.structure_ll for item in members),
            d2sigma_unpol_pb_per_gev2=sum(item.d2sigma_unpol_pb_per_gev2 for item in members),
            d2sigma_ll_pb_per_gev2=sum(item.d2sigma_ll_pb_per_gev2 for item in members),
        )
    return totals


def integrate_contributions(
    q2_min: float,
    q2_max: float,
    y_min: float,
    y_max: float,
    nq2: int,
    ny: int,
    inputs: EWInputs,
    beams: BeamInputs,
    sum_pdf,
    diff_pdf,
    flavors: Sequence[Flavor] | None = None,
    mode: LOCalculatorMode = DEFAULT_MODE_POLDIS,
    quark_masses_gev: Mapping[str, float] | None = None,
    mass_model: QuarkMassModel = DEFAULT_QUARK_MASS_MODEL,
) -> List[IntegratedContribution]:
    if q2_min <= 0.0 or q2_max <= 0.0 or q2_min >= q2_max:
        raise SystemExit("Require 0 < q2-min < q2-max.")
    if y_min <= 0.0 or y_max >= 1.0 or y_min >= y_max:
        raise SystemExit("Require 0 < y-min < y-max < 1.")
    if nq2 < 2 or ny < 2:
        raise SystemExit("Require nq2 >= 2 and ny >= 2.")

    import numpy as np

    nodes_q, weights_q = np.polynomial.legendre.leggauss(nq2)
    nodes_y, weights_y = np.polynomial.legendre.leggauss(ny)

    log_q2_min = math.log(q2_min)
    log_q2_max = math.log(q2_max)
    half_log_span = 0.5 * (log_q2_max - log_q2_min)
    log_mid = 0.5 * (log_q2_max + log_q2_min)
    y_mid = 0.5 * (y_max + y_min)
    y_half_span = 0.5 * (y_max - y_min)

    totals: Dict[tuple[str, str], List[float]] = {}
    active_flavors = list(flavors) if flavors is not None else FLAVORS
    masses = quark_masses_gev or {}

    for uq, wq in zip(nodes_q, weights_q):
        log_q2 = log_mid + half_log_span * float(uq)
        q2 = math.exp(log_q2)
        jac_q2 = half_log_span * q2
        for uy, wy in zip(nodes_y, weights_y):
            y = y_mid + y_half_span * float(uy)
            jac_y = y_half_span
            x_b = x_bjorken(q2, y, beams.s_hadronic)
            if x_b <= 0.0 or x_b >= 1.0:
                continue
            if not point_passes_lo_mode_cuts(q2, x_b, mode):
                continue
            prefactor = differential_prefactor_pb(q2, y, inputs.alpha_em)
            weight = float(wq) * float(wy) * jac_q2 * jac_y
            for flavor in active_flavors:
                pieces = channel_pieces(flavor, q2, inputs, mode)
                kernels = kernel_from_pieces(pieces, y)
                mass_gev, x_pdf, mass_weight = flavor_pdf_x_and_weight(
                    flavor, x_b, q2, masses, mass_model
                )
                kernel_ratios = (
                    massive_born_kernel_ratios(flavor, q2, y, inputs, mass_gev, mode)
                    if mass_model.use_massive_born_kernel
                    else {channel: (1.0, 1.0) for channel in CHANNEL_ORDER}
                )
                if x_pdf <= 0.0 or x_pdf >= 1.0:
                    continue
                xf_unpol = pdf_xf(sum_pdf, flavor.pdg_id, x_pdf, q2)
                xf_pol = pdf_xf(diff_pdf, flavor.pdg_id, x_pdf, q2)
                sigma_unpol_kernels = {
                    CHANNEL_GAMMA: kernels[CHANNEL_GAMMA].sigma_unpol_kernel * kernel_ratios[CHANNEL_GAMMA][0],
                    CHANNEL_Z: kernels[CHANNEL_Z].sigma_unpol_kernel * kernel_ratios[CHANNEL_Z][0],
                    CHANNEL_ALL: kernels[CHANNEL_ALL].sigma_unpol_kernel * kernel_ratios[CHANNEL_ALL][0],
                }
                sigma_ll_kernels = {
                    CHANNEL_GAMMA: kernels[CHANNEL_GAMMA].sigma_ll_kernel * kernel_ratios[CHANNEL_GAMMA][1],
                    CHANNEL_Z: kernels[CHANNEL_Z].sigma_ll_kernel * kernel_ratios[CHANNEL_Z][1],
                    CHANNEL_ALL: kernels[CHANNEL_ALL].sigma_ll_kernel * kernel_ratios[CHANNEL_ALL][1],
                }
                sigma_unpol_kernels[CHANNEL_INT] = (
                    sigma_unpol_kernels[CHANNEL_ALL]
                    - sigma_unpol_kernels[CHANNEL_GAMMA]
                    - sigma_unpol_kernels[CHANNEL_Z]
                )
                sigma_ll_kernels[CHANNEL_INT] = (
                    sigma_ll_kernels[CHANNEL_ALL]
                    - sigma_ll_kernels[CHANNEL_GAMMA]
                    - sigma_ll_kernels[CHANNEL_Z]
                )
                for channel in CHANNEL_ORDER:
                    sigma_unpol = prefactor * mass_weight * xf_unpol * sigma_unpol_kernels[channel]
                    sigma_ll = prefactor * mass_weight * xf_pol * sigma_ll_kernels[channel]
                    key = (flavor.name, channel)
                    total = totals.setdefault(key, [0.0, 0.0])
                    total[0] += weight * sigma_unpol
                    total[1] += weight * sigma_ll

    rows = [
        IntegratedContribution(
            flavor=flavor_name,
            channel=channel,
            sigma_unpol_pb=values[0],
            sigma_ll_pb=values[1],
        )
        for (flavor_name, channel), values in sorted(totals.items())
    ]
    return rows


def sum_integrated_contributions(
    rows: Iterable[IntegratedContribution],
) -> Dict[str, IntegratedContribution]:
    grouped: Dict[str, List[IntegratedContribution]] = {}
    for row in rows:
        grouped.setdefault(row.channel, []).append(row)
    return {
        channel: IntegratedContribution(
            flavor="TOTAL",
            channel=channel,
            sigma_unpol_pb=sum(item.sigma_unpol_pb for item in items),
            sigma_ll_pb=sum(item.sigma_ll_pb for item in items),
        )
        for channel, items in grouped.items()
    }


def serialize_dataclass_rows(rows: Iterable[object]) -> List[Dict[str, object]]:
    return [asdict(row) for row in rows]


def filter_rows_by_channel(rows: Iterable[object], channels: Sequence[str]) -> List[object]:
    allowed = set(channels)
    return [row for row in rows if getattr(row, "channel") in allowed]


def flavor_names(flavors: Sequence[Flavor] | None = None) -> List[str]:
    return [flavor.name for flavor in (list(flavors) if flavors is not None else FLAVORS)]


def lookup_flavors(names: Sequence[str]) -> List[Flavor]:
    missing = [name for name in names if name not in FLAVOR_LOOKUP]
    if missing:
        raise SystemExit(f"Unknown flavor names: {', '.join(missing)}")
    return [FLAVOR_LOOKUP[name] for name in names]


def w2_from_x_q2(x_b: float, q2: float) -> float:
    return (1.0 - x_b) * q2 / x_b


def point_passes_lo_mode_cuts(q2: float, x_b: float, mode: LOCalculatorMode) -> bool:
    if not mode.enforce_ep_x_cuts:
        return True
    x1 = 1.0
    x2 = x_b
    w2 = w2_from_x_q2(x_b, q2)
    if x1 <= mode.x1_min or x2 <= mode.x2_min:
        return False
    if w2 <= mode.w2_min:
        return False
    if math.isfinite(mode.w2_max) and w2 >= mode.w2_max:
        return False
    return True
