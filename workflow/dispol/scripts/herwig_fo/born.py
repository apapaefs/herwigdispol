from __future__ import annotations

from dataclasses import dataclass
import math

from .config import BeamConfig
from .ew import DEFAULT_EW_INPUTS, EWInputs, born_prefactor_charge, sigma_born_factor
from .hepmc import EventParticle, HepMCEvent, ME_STATUS_ROLE2
from .pdfs import PDFPair, ratios_for_flavor
from .vectors import FourVector


PB_PER_GEV2 = 0.389379338e9


@dataclass(frozen=True)
class DISPoint:
    q2: float
    y: float
    x_b: float
    flavor: int
    setup: str = "ALL"
    lepton_pol: float = 0.0
    hadron_pol: float = 0.0
    lepton_id: int = 11
    mu2: float | None = None

    @property
    def ell(self) -> float:
        return 2.0 / self.y - 1.0

    @property
    def scale2(self) -> float:
        return self.q2 if self.mu2 is None else self.mu2


@dataclass(frozen=True)
class LabKinematics:
    electron_in: FourVector
    proton_in: FourVector
    exchanged_boson: FourVector
    electron_out: FourVector
    parton_in: FourVector
    parton_out: FourVector
    remnant: FourVector


def xb_from_q2_y(q2: float, y: float, beams: BeamConfig = BeamConfig()) -> float:
    return q2 / (beams.s_ep * y)


def lab_kinematics(point: DISPoint, beams: BeamConfig = BeamConfig()) -> LabKinematics:
    ee = beams.lepton_energy
    ep = beams.hadron_energy
    q2 = point.q2
    qx = math.sqrt(max(0.0, q2 * (1.0 - point.y)))
    q0 = ee * point.y - q2 / (4.0 * ee)
    qz = -ee * point.y - q2 / (4.0 * ee)
    q = FourVector(qx, 0.0, qz, q0)

    electron_in = FourVector(0.0, 0.0, -ee, ee)
    proton_in = FourVector(0.0, 0.0, ep, ep)
    electron_out = electron_in - q
    parton_in = FourVector(0.0, 0.0, point.x_b * ep, point.x_b * ep)
    parton_out = parton_in + q
    remnant = proton_in - parton_in
    return LabKinematics(electron_in, proton_in, q, electron_out, parton_in, parton_out, remnant)


def born_weight(
    point: DISPoint,
    pdfs: PDFPair,
    ew_inputs: EWInputs = DEFAULT_EW_INPUTS,
) -> float:
    ratios = ratios_for_flavor(pdfs, point.flavor, point.x_b, 1.0, point.scale2, point.hadron_pol)
    angular = sigma_born_factor(
        point.setup,
        point.lepton_id,
        point.flavor,
        point.q2,
        point.lepton_pol,
        ratios.pq_born,
        point.ell,
    )
    prefactor = 2.0 * math.pi * ew_inputs.alpha_em * ew_inputs.alpha_em * PB_PER_GEV2
    xf = pdfs.unpolarized.xfx(point.flavor, point.x_b, point.scale2)
    charge = born_prefactor_charge(point.setup, point.lepton_id, point.flavor)
    return prefactor * xf * charge * (0.5 * point.y * point.y * angular) / (point.y * point.q2 * point.q2)


def born_event(
    point: DISPoint,
    weight: float,
    event_number: int = 0,
    beams: BeamConfig = BeamConfig(),
    role: int = 2,
) -> HepMCEvent:
    kin = lab_kinematics(point, beams)
    status = ME_STATUS_ROLE2 if role == 2 else 1
    outgoing = (
        EventParticle(point.lepton_id, kin.electron_out, 1),
        EventParticle(point.flavor, kin.parton_out, status, role=role),
        EventParticle(82, kin.remnant, 1),
    )
    return HepMCEvent(
        event_number=event_number,
        weight=weight,
        incoming=(
            EventParticle(beams.hadron_id, kin.proton_in, 4),
            EventParticle(point.lepton_id, kin.electron_in, 4),
        ),
        outgoing=outgoing,
        attributes={
            "HERWIG_FO_EVENT_CLASS": "BORN_PROJECTED",
            "HERWIG_FO_SETUP": point.setup,
        },
    )
