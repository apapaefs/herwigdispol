from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math


DEFAULT_OBSERVABLE_TOKENS = ("Q2", "XBj", "pT1", "pT2", "Mjj", "Eta", "Zeta", "D", "ALL")


@dataclass(frozen=True)
class BinComparison:
    path: str
    index: int
    herwig: float
    reference: float
    herwig_error: float
    reference_error: float
    pull: float


@dataclass(frozen=True)
class ComparisonSummary:
    ok: bool
    max_abs_pull: float
    bins: tuple[BinComparison, ...]
    missing_in_herwig: tuple[str, ...]
    missing_in_reference: tuple[str, ...]


def _read_yoda(path: Path):
    try:
        import yoda
    except ImportError as exc:
        raise RuntimeError("Could not import yoda. Load a Rivet/YODA Python environment first.") from exc
    if hasattr(yoda, "read"):
        return yoda.read(str(path))
    if hasattr(yoda, "readYODA"):
        return yoda.readYODA(str(path))
    raise RuntimeError("The imported yoda module does not expose read/readYODA.")


def _objects_by_path(path: Path) -> dict[str, object]:
    raw = _read_yoda(path)
    if isinstance(raw, dict):
        return {str(k): v for k, v in raw.items()}
    out = {}
    for obj in raw:
        obj_path = obj.path() if callable(getattr(obj, "path", None)) else getattr(obj, "path", None)
        if obj_path:
            out[str(obj_path)] = obj
    return out


def _bin_value(bin_obj) -> tuple[float, float]:
    value = bin_obj.val() if callable(getattr(bin_obj, "val", None)) else float(getattr(bin_obj, "sumW", 0.0))
    if callable(getattr(bin_obj, "errAvg", None)):
        error = bin_obj.errAvg()
    elif callable(getattr(bin_obj, "err", None)):
        error = bin_obj.err()
    else:
        error = 0.0
    return float(value), float(error)


def _bins(obj):
    if callable(getattr(obj, "bins", None)):
        return obj.bins()
    return []


def compare_yoda(
    herwig_yoda: Path,
    reference_yoda: Path,
    max_pull: float = 3.0,
    observable_tokens: tuple[str, ...] = DEFAULT_OBSERVABLE_TOKENS,
) -> ComparisonSummary:
    herwig = _objects_by_path(herwig_yoda)
    reference = _objects_by_path(reference_yoda)
    selected_reference = {
        path: obj for path, obj in reference.items() if any(token in path for token in observable_tokens)
    }
    selected_herwig = {
        path: obj for path, obj in herwig.items() if any(token in path for token in observable_tokens)
    }
    common = sorted(set(selected_herwig) & set(selected_reference))
    missing_in_herwig = tuple(sorted(set(selected_reference) - set(selected_herwig)))
    missing_in_reference = tuple(sorted(set(selected_herwig) - set(selected_reference)))

    comparisons: list[BinComparison] = []
    for path in common:
        h_bins = list(_bins(selected_herwig[path]))
        r_bins = list(_bins(selected_reference[path]))
        for idx, (h_bin, r_bin) in enumerate(zip(h_bins, r_bins), start=1):
            h_val, h_err = _bin_value(h_bin)
            r_val, r_err = _bin_value(r_bin)
            sigma = math.hypot(h_err, r_err)
            pull = 0.0 if sigma == 0.0 and h_val == r_val else (h_val - r_val) / sigma if sigma > 0.0 else math.inf
            comparisons.append(BinComparison(path, idx, h_val, r_val, h_err, r_err, pull))

    largest = max((abs(item.pull) for item in comparisons), default=0.0)
    ok = largest <= max_pull and not missing_in_herwig
    return ComparisonSummary(ok, largest, tuple(comparisons), missing_in_herwig, missing_in_reference)
