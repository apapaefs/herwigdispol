from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Callable, Sequence


def matchbox_beta_coefficients(nf: int) -> tuple[float, float]:
    beta0 = (33.0 - 2.0 * nf) / (12.0 * math.pi)
    beta1 = (153.0 - 19.0 * nf) / (24.0 * math.pi * math.pi)
    return beta0, beta1


def bisect_root(function: Callable[[float], float], low: float, high: float, label: str) -> float:
    f_low = function(low)
    f_high = function(high)
    if f_low == 0.0:
        return low
    if f_high == 0.0:
        return high
    if f_low * f_high > 0.0:
        raise ValueError(
            f"Could not bracket root for {label}: f({low:.6g})={f_low:.6g}, f({high:.6g})={f_high:.6g}"
        )
    for _ in range(200):
        mid = 0.5 * (low + high)
        f_mid = function(mid)
        if f_mid == 0.0:
            return mid
        if f_low * f_mid <= 0.0:
            high = mid
            f_high = f_mid
        else:
            low = mid
            f_low = f_mid
    return 0.5 * (low + high)


def matchbox_large_scale_alpha(scale_q2: float, lambda2: float, nf: int, two_largeq_terms: bool = True) -> float:
    if scale_q2 <= lambda2:
        raise ValueError(
            f"Matchbox NLOAlphaS requires Q2 > lambda2. Got Q2={scale_q2:.6g}, lambda2={lambda2:.6g}."
        )
    beta0, beta1 = matchbox_beta_coefficients(nf)
    slog = math.log(scale_q2 / lambda2)
    leading = 1.0 / (beta0 * slog)
    result = leading * (1.0 - (beta1 / (beta0 * beta0)) * math.log(slog) / slog)
    if two_largeq_terms:
        result += leading * (beta1 / (beta0 * beta0 * slog)) ** 2 * ((math.log(slog) - 0.5) ** 2 - 1.25)
    return result


def default_matchbox_quark_masses2() -> list[float]:
    masses = [0.0, 0.005, 0.0023, 0.095, 1.25, 4.2, 174.2]
    masses2 = [mass * mass for mass in masses]
    if masses2[1] > masses2[2]:
        masses2[1], masses2[2] = masses2[2], masses2[1]
    return masses2


def active_flavours(scale_q2: float, quark_masses2: Sequence[float], max_active_flavours: int) -> int:
    active = 0
    if scale_q2 > 0.0:
        while quark_masses2[active] < scale_q2:
            active += 1
            if active == max_active_flavours + 1:
                break
        active -= 1
    return active


def solve_lambda2_for_input_alpha(input_alpha_s: float, input_scale: float, nf: int) -> float:
    input_q2 = float(input_scale) * float(input_scale)

    def equation(lambda2: float) -> float:
        return matchbox_large_scale_alpha(input_q2, lambda2, nf) - input_alpha_s

    return bisect_root(equation, 1.0e-6, 1.0, f"lambda2 for nf={nf} at input scale")


def solve_matched_lambda2(target_alpha: float, threshold_q2: float, nf: int) -> float:
    def equation(lambda2: float) -> float:
        return matchbox_large_scale_alpha(threshold_q2, lambda2, nf) - target_alpha

    return bisect_root(equation, 1.0e-6, 1.0, f"matched lambda2 for nf={nf}")


def solve_matchbox_lambda_squared(
    input_alpha_s: float,
    input_scale: float,
    min_active_flavours: int = 3,
    max_active_flavours: int = 6,
) -> list[float]:
    quark_masses2 = default_matchbox_quark_masses2()
    lambdas2 = [0.0] * 7
    input_q2 = float(input_scale) * float(input_scale)
    active_at_input = active_flavours(input_q2, quark_masses2, max_active_flavours)
    lambdas2[active_at_input] = solve_lambda2_for_input_alpha(input_alpha_s, input_scale, active_at_input)

    below = active_at_input
    while below > min_active_flavours:
        threshold_q2 = quark_masses2[below]
        target_alpha = matchbox_large_scale_alpha(threshold_q2, lambdas2[below], below)
        lambdas2[below - 1] = solve_matched_lambda2(target_alpha, threshold_q2, below - 1)
        below -= 1

    above = active_at_input
    while above < max_active_flavours:
        threshold_q2 = quark_masses2[above + 1]
        target_alpha = matchbox_large_scale_alpha(threshold_q2, lambdas2[above], above)
        lambdas2[above + 1] = solve_matched_lambda2(target_alpha, threshold_q2, above + 1)
        above += 1

    for flavor in range(min_active_flavours):
        lambdas2[flavor] = lambdas2[min_active_flavours]
    for flavor in range(max_active_flavours + 1, 7):
        lambdas2[flavor] = lambdas2[max_active_flavours]
    return lambdas2


def matchbox_running_alpha(
    q2: float,
    lambda_squared: Sequence[float],
    quark_masses2: Sequence[float],
    min_active_flavours: int = 3,
    max_active_flavours: int = 6,
) -> float:
    nf = active_flavours(float(q2), quark_masses2, max_active_flavours)
    nf = max(min_active_flavours, min(max_active_flavours, nf))
    return matchbox_large_scale_alpha(float(q2), lambda_squared[nf], nf)


@dataclass(frozen=True)
class MatchboxAlphaS:
    alpha_s_mz: float = 0.118
    mz: float = 91.188
    min_active_flavours: int = 3
    max_active_flavours: int = 6
    quark_masses2: tuple[float, ...] = field(default_factory=lambda: tuple(default_matchbox_quark_masses2()))
    lambda_squared: tuple[float, ...] = field(init=False)

    def __post_init__(self) -> None:
        lambdas = solve_matchbox_lambda_squared(
            self.alpha_s_mz,
            self.mz,
            self.min_active_flavours,
            self.max_active_flavours,
        )
        object.__setattr__(self, "lambda_squared", tuple(lambdas))

    def alpha_s(self, q2: float) -> float:
        return matchbox_running_alpha(
            q2,
            self.lambda_squared,
            self.quark_masses2,
            self.min_active_flavours,
            self.max_active_flavours,
        )
