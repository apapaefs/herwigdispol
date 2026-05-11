from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class EWInputs:
    mz: float = 91.1876
    width_z: float = 2.4952
    sin2_theta_w: float = 1.0 - (80.450 / 91.1876) ** 2
    alpha_em: float = 0.00729927

    @property
    def sin_theta_w(self) -> float:
        return math.sqrt(self.sin2_theta_w)

    @property
    def cos_theta_w(self) -> float:
        return math.sqrt(1.0 - self.sin2_theta_w)


@dataclass(frozen=True)
class NCCoefficients:
    N0: float = 0.0
    Nl: float = 0.0
    Nq: float = 0.0
    Nlq: float = 0.0
    D0: float = 0.0
    Dl: float = 0.0
    Dq: float = 0.0
    Dlq: float = 0.0


@dataclass(frozen=True)
class CollinearBlendWeights:
    q_unpolarized: float
    q_polarized: float
    g_unpolarized: float
    g_polarized: float


@dataclass(frozen=True)
class NeutralCurrentResponse:
    q_odd_response: float
    g_odd_response: float


DEFAULT_EW_INPUTS = EWInputs()


def _setup_code(setup: str) -> int:
    value = setup.upper()
    if value == "ALL":
        return 0
    if value == "GAMMA":
        return 1
    if value == "Z":
        return 2
    raise ValueError(f"Unsupported NC setup '{setup}'.")


def _fermion_charge(pid: int) -> float:
    apid = abs(pid)
    if apid in (11, 13, 15):
        return -1.0
    if apid in (12, 14, 16):
        return 0.0
    if apid in (2, 4, 6):
        return 2.0 / 3.0
    if apid in (1, 3, 5):
        return -1.0 / 3.0
    raise ValueError(f"Unsupported fermion id {pid}.")


def _fermion_t3(pid: int) -> float:
    apid = abs(pid)
    if apid in (12, 14, 16, 2, 4, 6):
        return 0.5
    if apid in (11, 13, 15, 1, 3, 5):
        return -0.5
    raise ValueError(f"Unsupported fermion id {pid}.")


def _herwig_vector_axial(pid: int, inputs: EWInputs) -> tuple[float, float]:
    charge = _fermion_charge(pid)
    t3 = _fermion_t3(pid)
    vector = 0.5 * (t3 - 2.0 * charge * inputs.sin2_theta_w)
    axial = 0.5 * t3
    return vector, axial


def nc_coefficients(setup: str, lepton_id: int, quark_id: int, q2: float, inputs: EWInputs = DEFAULT_EW_INPUTS) -> NCCoefficients:
    gamma_z = _setup_code(setup)
    eta_l = -1.0 if lepton_id < 0 else 1.0
    eta_q = -1.0 if quark_id < 0 else 1.0

    ql = _fermion_charge(lepton_id) if gamma_z == 0 or gamma_z == 1 else 0.0
    qq = _fermion_charge(quark_id) if gamma_z == 0 or gamma_z == 1 else 0.0
    cv_l, ca_l = _herwig_vector_axial(lepton_id, inputs)
    cv_q, ca_q = _herwig_vector_axial(quark_id, inputs)

    ql2 = ql * ql
    qq2 = qq * qq
    cvl2 = cv_l * cv_l
    cal2 = ca_l * ca_l
    cvq2 = cv_q * cv_q
    caq2 = ca_q * ca_q
    cvl_cvq = cv_l * cv_q
    cal_caq = ca_l * ca_q

    k = 1.0 / (inputs.sin_theta_w * inputs.cos_theta_w) ** 2
    z_int = 0.0
    z_sq = 0.0
    if gamma_z == 0 or gamma_z == 2:
        mz2 = inputs.mz * inputs.mz
        gamma_width = inputs.mz * inputs.width_z
        den = (q2 + mz2) * (q2 + mz2) + gamma_width * gamma_width
        r_int = q2 * (q2 + mz2) / den
        r_sq = q2 * q2 / den
        z_int = k * r_int
        z_sq = k * k * r_sq

    return NCCoefficients(
        N0=eta_l * eta_q * 4.0 * cal_caq * (z_int * ql * qq + 2.0 * z_sq * cvl_cvq),
        Nl=eta_q * -4.0 * ca_q * (z_sq * cv_q * (cal2 + cvl2) + z_int * cv_l * ql * qq),
        Nq=eta_l * -4.0 * ca_l * (z_sq * cv_l * (caq2 + cvq2) + z_int * cv_q * ql * qq),
        Nlq=2.0 * (z_sq * (caq2 + cvq2) * (cal2 + cvl2) + 2.0 * z_int * cvl_cvq * ql * qq + ql2 * qq2),
        D0=ql2 * qq2 + 2.0 * z_int * ql * qq * cvl_cvq + z_sq * (cvl2 + cal2) * (cvq2 + caq2),
        Dl=eta_l * -2.0 * ca_l * (z_sq * cv_l * (caq2 + cvq2) + z_int * cv_q * ql * qq),
        Dq=eta_q * -2.0 * ca_q * (z_sq * cv_q * (cal2 + cvl2) + z_int * cv_l * ql * qq),
        Dlq=eta_l * eta_q * 2.0 * cal_caq * (2.0 * z_sq * cvl_cvq + z_int * ql * qq),
    )


def a_pol(setup: str, lepton_id: int, quark_id: int, q2: float, lepton_pol: float, quark_pol: float) -> float:
    if _setup_code(setup) == 1:
        return 2.0 * lepton_pol * quark_pol
    coeff = nc_coefficients(setup, lepton_id, quark_id, q2)
    num = coeff.N0 + lepton_pol * coeff.Nl + quark_pol * coeff.Nq + lepton_pol * quark_pol * coeff.Nlq
    den = coeff.D0 + lepton_pol * coeff.Dl + quark_pol * coeff.Dq + lepton_pol * quark_pol * coeff.Dlq
    if den == 0.0:
        return 0.0
    return num / den


def sigma_born_factor(
    setup: str,
    lepton_id: int,
    quark_id: int,
    q2: float,
    lepton_pol: float,
    quark_pol: float,
    ell: float,
) -> float:
    if _setup_code(setup) == 1:
        return 1.0 + 2.0 * lepton_pol * quark_pol * ell + ell * ell
    coeff = nc_coefficients(setup, lepton_id, quark_id, q2)
    d_even = coeff.D0 + lepton_pol * coeff.Dl
    d_spin = quark_pol * (coeff.Dq + lepton_pol * coeff.Dlq)
    n_even = coeff.N0 + lepton_pol * coeff.Nl
    n_spin = quark_pol * (coeff.Nq + lepton_pol * coeff.Nlq)
    return (1.0 + ell * ell) * (d_even + d_spin) + ell * (n_even + n_spin)


def born_prefactor_charge(setup: str, lepton_id: int, quark_id: int) -> float:
    if _setup_code(setup) != 1:
        return 1.0
    ql = _fermion_charge(lepton_id)
    qq = _fermion_charge(quark_id)
    return ql * ql * qq * qq


def collinear_blend_weights(
    setup: str,
    lepton_id: int,
    quark_id: int,
    q2: float,
    lepton_pol: float,
    quark_pol: float,
    ell: float,
) -> CollinearBlendWeights:
    if _setup_code(setup) == 1:
        a = a_pol(setup, lepton_id, quark_id, q2, lepton_pol, quark_pol)
        denom = 1.0 + a * ell + ell * ell
        f = a * ell / denom if abs(denom) > 1e-30 else 0.0
        return CollinearBlendWeights(1.0 - f, f, 1.0 - f, f)

    coeff = nc_coefficients(setup, lepton_id, quark_id, q2)
    d_even = coeff.D0 + lepton_pol * coeff.Dl
    d_spin = quark_pol * (coeff.Dq + lepton_pol * coeff.Dlq)
    n_even = coeff.N0 + lepton_pol * coeff.Nl
    n_spin = quark_pol * (coeff.Nq + lepton_pol * coeff.Nlq)
    sigma = (1.0 + ell * ell) * (d_even + d_spin) + ell * (n_even + n_spin)
    if abs(sigma) <= 1e-30:
        return CollinearBlendWeights(1.0, 0.0, 1.0, 0.0)
    return CollinearBlendWeights(
        ((1.0 + ell * ell) * d_even + ell * n_even) / sigma,
        ((1.0 + ell * ell) * d_spin + ell * n_spin) / sigma,
        ((1.0 + ell * ell) * d_even) / sigma,
        (ell * n_spin) / sigma,
    )


def qcdc_mapped_denominator_ratio(
    setup: str,
    lepton_id: int,
    quark_id: int,
    q2: float,
    lepton_pol: float,
    quark_pol_born: float,
    quark_pol_mapped: float,
) -> float:
    if _setup_code(setup) == 1:
        return 1.0
    coeff = nc_coefficients(setup, lepton_id, quark_id, q2)
    born = coeff.D0 + lepton_pol * coeff.Dl + quark_pol_born * coeff.Dq + lepton_pol * quark_pol_born * coeff.Dlq
    if abs(born) <= 1e-30:
        return 1.0
    mapped = coeff.D0 + lepton_pol * coeff.Dl + quark_pol_mapped * coeff.Dq + lepton_pol * quark_pol_mapped * coeff.Dlq
    return mapped / born


def real_emission_denominator_factor(
    setup: str,
    lepton_id: int,
    quark_id: int,
    q2: float,
    lepton_pol: float,
    mapped_parton_pol: float,
) -> float:
    if _setup_code(setup) == 1:
        return 1.0
    coeff = nc_coefficients(setup, lepton_id, quark_id, q2)
    return (
        coeff.D0
        + lepton_pol * coeff.Dl
        + mapped_parton_pol * coeff.Dq
        + lepton_pol * mapped_parton_pol * coeff.Dlq
    )


def neutral_current_response(
    setup: str,
    lepton_id: int,
    quark_id: int,
    q2: float,
    lepton_pol: float,
    quark_pol_born: float,
    ell: float,
) -> NeutralCurrentResponse | None:
    if _setup_code(setup) == 1:
        return None
    coeff = nc_coefficients(setup, lepton_id, quark_id, q2)
    d_even = coeff.D0 + lepton_pol * coeff.Dl
    d_spin = quark_pol_born * (coeff.Dq + lepton_pol * coeff.Dlq)
    n_even = coeff.N0 + lepton_pol * coeff.Nl
    n_spin = quark_pol_born * (coeff.Nq + lepton_pol * coeff.Nlq)
    born = (1.0 + ell * ell) * (d_even + d_spin) + ell * (n_even + n_spin)
    if not math.isfinite(born) or abs(born) <= 1e-30:
        return None
    q_polarized = ((1.0 + ell * ell) * d_spin + ell * n_spin) / born
    g_polarized = (ell * n_spin) / born
    if abs(quark_pol_born) <= 1e-30:
        return NeutralCurrentResponse(0.0, 0.0)
    return NeutralCurrentResponse(q_polarized / quark_pol_born, g_polarized / quark_pol_born)
