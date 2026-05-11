from __future__ import annotations

from dataclasses import dataclass
import math
import os
from pathlib import Path
from typing import Protocol

from .config import PDFProfile


DEFAULT_LHAPDF_DATA_PATHS = (
    "/Users/apapaefs/Projects/Herwig/Herwig-stable-gcc-full/share/LHAPDF",
    "/opt/homebrew/Cellar/lhapdf/6.5.4/share/LHAPDF",
)


class PDFProvider(Protocol):
    name: str

    def xfx(self, pid: int, x: float, q2: float) -> float:
        ...

    def alpha_s(self, q2: float) -> float | None:
        ...


@dataclass(frozen=True)
class PDFPair:
    unpolarized: PDFProvider
    difference: PDFProvider | None = None


@dataclass(frozen=True)
class PDFRatios:
    lo_pdf: float
    q_pdf: float
    g_pdf: float
    dq_pdf: float
    dg_pdf: float
    dlo_pdf: float
    q_ratio: float
    g_ratio: float
    deltaq_over_lo: float
    deltag_over_lo: float
    deltaq_over_lo_eff: float
    deltag_over_lo_eff: float
    pq_born: float
    pq_mapped: float
    pg_mapped: float
    has_difference_pdf: bool
    has_stable_difference_ratio: bool


@dataclass(frozen=True)
class LHAPDFProvider:
    name: str
    member: int = 0

    def __post_init__(self) -> None:
        if "LHAPDF_DATA_PATH" not in os.environ:
            for candidate in DEFAULT_LHAPDF_DATA_PATHS:
                if Path(candidate).exists():
                    os.environ["LHAPDF_DATA_PATH"] = candidate
                    break
        try:
            import lhapdf
        except ImportError as exc:
            raise RuntimeError(
                "Could not import LHAPDF Python bindings. Install/load LHAPDF or pass --toy-pdfs for tests."
            ) from exc
        try:
            pdf = lhapdf.mkPDF(self.name, self.member)
        except Exception as exc:
            raise RuntimeError(
                f"Could not load LHAPDF set '{self.name}' member {self.member}. "
                f"Check LHAPDF_DATA_PATH={os.environ.get('LHAPDF_DATA_PATH', '<unset>')}."
            ) from exc
        object.__setattr__(self, "_pdf", pdf)

    def xfx(self, pid: int, x: float, q2: float) -> float:
        return float(self._pdf.xfxQ2(pid, x, q2))

    def alpha_s(self, q2: float) -> float | None:
        return float(self._pdf.alphasQ2(q2))


@dataclass(frozen=True)
class ToyPDFProvider:
    name: str = "toy"

    def xfx(self, pid: int, x: float, q2: float) -> float:
        if x <= 0.0 or x >= 1.0:
            return 0.0
        apid = abs(pid)
        logq = math.log(max(q2, 1.0) / 100.0)
        if pid == 21:
            norm = 4.5
            eta = 4.0
            helicity = 0.20
        elif apid in (1, 2):
            norm = 1.2 if apid == 2 else 0.9
            eta = 3.2
            helicity = 0.35 if apid == 2 else -0.18
        elif apid in (3, 4, 5):
            norm = 0.18
            eta = 5.0
            helicity = 0.08 if apid % 2 == 0 else -0.06
        else:
            return 0.0
        sea_factor = 0.55 if pid < 0 and pid != 21 else 1.0
        base = norm * sea_factor * x ** (-0.22) * (1.0 - x) ** eta * (1.0 + 0.015 * logq)
        if "diff" in self.name.lower() or "pol" in self.name.lower() or "difference" in self.name.lower():
            anti_sign = -1.0 if pid < 0 and pid != 21 else 1.0
            return anti_sign * helicity * base
        return base

    def alpha_s(self, q2: float) -> float | None:
        return None


def load_pdf_pair(profile: PDFProfile, use_toy: bool = False) -> PDFPair:
    if use_toy:
        return PDFPair(ToyPDFProvider("toy-sum"), ToyPDFProvider("toy-difference"))
    return PDFPair(
        LHAPDFProvider(profile.unpolarized),
        LHAPDFProvider(profile.polarized_difference),
    )


def density(pdf: PDFProvider | None, pid: int, x: float, q2: float) -> float:
    if pdf is None or x <= 0.0 or x >= 1.0:
        return 0.0
    return pdf.xfx(pid, x, q2) / x


def clamp(value: float, low: float = -1.0, high: float = 1.0) -> float:
    return max(low, min(high, value))


def ratios_for_flavor(
    pdfs: PDFPair,
    flavor: int,
    x_b: float,
    x_p: float,
    q2: float,
    hadron_pol: float,
    ratio_floor: float = 1.0e-12,
) -> PDFRatios:
    if x_p <= 0.0:
        raise ValueError("x_p must be positive.")
    mapped_x = x_b / x_p
    lo_pdf = density(pdfs.unpolarized, flavor, x_b, q2)
    if abs(lo_pdf) <= ratio_floor:
        raise ZeroDivisionError(f"Zero unpolarized Born PDF for flavor {flavor} at x={x_b:.6g}, Q2={q2:.6g}.")

    q_pdf = density(pdfs.unpolarized, flavor, mapped_x, q2)
    g_pdf = density(pdfs.unpolarized, 21, mapped_x, q2)
    dq_pdf = 0.0
    dg_pdf = 0.0
    dlo_pdf = 0.0
    has_difference_pdf = pdfs.difference is not None and abs(hadron_pol) > 1e-12
    if has_difference_pdf:
        dq_pdf = density(pdfs.difference, flavor, mapped_x, q2)
        dg_pdf = density(pdfs.difference, 21, mapped_x, q2)
        dlo_pdf = density(pdfs.difference, flavor, x_b, q2)

    pq_born = clamp(hadron_pol * dlo_pdf / lo_pdf) if has_difference_pdf else 0.0
    pq_mapped = clamp(hadron_pol * dq_pdf / q_pdf) if has_difference_pdf and abs(q_pdf) > ratio_floor else pq_born
    pg_mapped = clamp(hadron_pol * dg_pdf / g_pdf) if has_difference_pdf and abs(g_pdf) > ratio_floor else 0.0

    min_dlo = max(ratio_floor, 1.0e-4 * abs(lo_pdf))
    has_stable_difference_ratio = has_difference_pdf and abs(dlo_pdf) > min_dlo
    deltaq_over_lo = hadron_pol * dq_pdf / lo_pdf if has_difference_pdf else 0.0
    deltag_over_lo = hadron_pol * dg_pdf / lo_pdf if has_difference_pdf else 0.0
    if has_stable_difference_ratio:
        deltaq_over_lo_eff = pq_born * dq_pdf / dlo_pdf
        deltag_over_lo_eff = pq_born * dg_pdf / dlo_pdf
    else:
        deltaq_over_lo_eff = 0.0
        deltag_over_lo_eff = 0.0

    return PDFRatios(
        lo_pdf=lo_pdf,
        q_pdf=q_pdf,
        g_pdf=g_pdf,
        dq_pdf=dq_pdf,
        dg_pdf=dg_pdf,
        dlo_pdf=dlo_pdf,
        q_ratio=q_pdf / lo_pdf,
        g_ratio=g_pdf / lo_pdf,
        deltaq_over_lo=deltaq_over_lo,
        deltag_over_lo=deltag_over_lo,
        deltaq_over_lo_eff=deltaq_over_lo_eff,
        deltag_over_lo_eff=deltag_over_lo_eff,
        pq_born=pq_born,
        pq_mapped=pq_mapped,
        pg_mapped=pg_mapped,
        has_difference_pdf=has_difference_pdf,
        has_stable_difference_ratio=has_stable_difference_ratio,
    )
