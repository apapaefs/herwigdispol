from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping


NC_SETUPS = ("GAMMA", "Z", "ALL")

HELICITIES: Mapping[str, tuple[float, float]] = {
    "00": (0.0, 0.0),
    "PP": (1.0, 1.0),
    "PM": (1.0, -1.0),
    "MP": (-1.0, 1.0),
    "MM": (-1.0, -1.0),
}

ACTIVE_QUARK_FLAVORS = (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5)


@dataclass(frozen=True)
class BeamConfig:
    lepton_energy: float = 18.0
    hadron_energy: float = 275.0
    lepton_id: int = 11
    hadron_id: int = 2212

    @property
    def s_ep(self) -> float:
        return 4.0 * self.lepton_energy * self.hadron_energy


@dataclass(frozen=True)
class DISCuts:
    q2_min: float = 100.0
    q2_max: float = 2500.0
    y_min: float = 0.2
    y_max: float = 0.6


CUT_WINDOWS: Mapping[str, DISCuts] = {
    "validation": DISCuts(),
    "plain66": DISCuts(q2_min=49.0),
}


@dataclass(frozen=True)
class PDFProfile:
    name: str
    unpolarized: str
    polarized_difference: str


PDF_PROFILES: Mapping[str, PDFProfile] = {
    "nnpdf_paired": PDFProfile(
        name="nnpdf_paired",
        unpolarized="NNPDF40_nlo_pch_as_01180",
        polarized_difference="NNPDFpol20_nlo_as_01180",
    ),
}


@dataclass(frozen=True)
class RunConfig:
    beams: BeamConfig = BeamConfig()
    cuts: DISCuts = DISCuts()
    pdf_profile: PDFProfile = PDF_PROFILES["nnpdf_paired"]
    setups: tuple[str, ...] = NC_SETUPS
    helicities: tuple[str, ...] = tuple(HELICITIES.keys())
    alpha_s_mz: float = 0.118
    mz: float = 91.188


def default_run_config() -> RunConfig:
    return RunConfig()


def cuts_for_window(name: str) -> DISCuts:
    key = name.lower()
    if key not in CUT_WINDOWS:
        raise ValueError(f"Unsupported DIS cut window '{name}'. Choose one of {', '.join(CUT_WINDOWS)}.")
    return CUT_WINDOWS[key]


def require_setup(setup: str) -> str:
    value = setup.upper()
    if value not in NC_SETUPS:
        raise ValueError(f"Unsupported neutral-current setup '{setup}'. Choose one of {', '.join(NC_SETUPS)}.")
    return value


def require_helicity(label: str) -> tuple[float, float]:
    key = label.upper()
    if key not in HELICITIES:
        raise ValueError(f"Unsupported helicity '{label}'. Choose one of {', '.join(HELICITIES)}.")
    return HELICITIES[key]


def require_pdf_profile(name: str) -> PDFProfile:
    key = name.lower()
    if key not in PDF_PROFILES:
        raise ValueError(f"Unsupported PDF profile '{name}'. Choose one of {', '.join(PDF_PROFILES)}.")
    return PDF_PROFILES[key]
