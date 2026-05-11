from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class FourVector:
    px: float
    py: float
    pz: float
    e: float

    def __add__(self, other: "FourVector") -> "FourVector":
        return FourVector(self.px + other.px, self.py + other.py, self.pz + other.pz, self.e + other.e)

    def __sub__(self, other: "FourVector") -> "FourVector":
        return FourVector(self.px - other.px, self.py - other.py, self.pz - other.pz, self.e - other.e)

    def scale(self, factor: float) -> "FourVector":
        return FourVector(self.px * factor, self.py * factor, self.pz * factor, self.e * factor)

    def pt(self) -> float:
        return math.hypot(self.px, self.py)

    def p_abs(self) -> float:
        return math.sqrt(max(0.0, self.px * self.px + self.py * self.py + self.pz * self.pz))

    def mass2(self) -> float:
        return self.e * self.e - self.px * self.px - self.py * self.py - self.pz * self.pz

    def mass(self) -> float:
        m2 = self.mass2()
        if abs(m2) < 1e-10:
            return 0.0
        return math.copysign(math.sqrt(abs(m2)), m2)


def massless(px: float, py: float, pz: float) -> FourVector:
    return FourVector(px, py, pz, math.sqrt(px * px + py * py + pz * pz))


def dot3(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross3(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm3(v: tuple[float, float, float]) -> float:
    return math.sqrt(max(0.0, dot3(v, v)))


def unit3(v: tuple[float, float, float], fallback: tuple[float, float, float]) -> tuple[float, float, float]:
    n = norm3(v)
    if n <= 0.0:
        return fallback
    return (v[0] / n, v[1] / n, v[2] / n)


def combine3(
    ax: float,
    a: tuple[float, float, float],
    bx: float,
    b: tuple[float, float, float],
    cx: float,
    c: tuple[float, float, float],
) -> tuple[float, float, float]:
    return (
        ax * a[0] + bx * b[0] + cx * c[0],
        ax * a[1] + bx * b[1] + cx * c[1],
        ax * a[2] + bx * b[2] + cx * c[2],
    )
