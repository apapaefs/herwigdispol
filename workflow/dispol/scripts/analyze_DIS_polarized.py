#!/usr/bin/env python3.10
"""Build DIS polarized YODA objects from helicity-resolved Herwig YODAs."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy
import yoda

from poldis_top_to_yoda import new_binned_estimate, set_bin_val_err


BASE_DIS_LABELS = (
    "Q2",
    "Q2PreCut",
    "Pt",
    "Mjj",
    "Eta",
    "Zeta",
    "pT1",
    "pT1PreCut",
    "pT2",
    "pT2PreCut",
    "pT2OverpT1",
    "pTAsym",
    "XBj",
    "XBjPreCut",
    "YPreCut",
)
PS_EXTRA_LABELS = (
    "NJets",
    "pT3",
    "pT3OverpT1",
    "SumPtExtra",
    "Phi3",
    "PhiCurrentHemi",
    "Broadening",
)
ANALYSIS_LABELS = {
    "MC_DIS_BREIT": BASE_DIS_LABELS,
    "MC_DIS_PS": BASE_DIS_LABELS + PS_EXTRA_LABELS,
}
ANALYSIS_ALL_LABELS = {
    "MC_DIS_PS": (
        "pT2OverpT1",
        "pTAsym",
        "pT3",
        "pT3OverpT1",
        "Phi3",
        "PhiCurrentHemi",
        "Broadening",
    ),
}
ANALYSIS_DERIVED_MOMENTS = {
    "MC_DIS_PS": (
        ("Cos2PhiCurrentHemiVsQ2", "Cos2PhiCurrentHemiNumQ2", "Cos2PhiCurrentHemiDenQ2"),
    ),
}
PS_SETUPS = ("SPINVAL", "SPINCOMP")
SETUPS = ("ALL", "GAMMA", "Z", "CC", *PS_SETUPS)
FIVE_HELICITY_SETUPS = ("ALL", "Z", "CC", *PS_SETUPS)
SCALE_VARIATION_FACTORS = {
    "nominal": 1.0,
    "ScaleFactorDown": 0.5,
    "ScaleFactorUp": 2.0,
}


InputSpec = Optional[Union[str, Sequence[str]]]


def _normalize_paths(paths: InputSpec) -> List[str]:
    if not paths:
        return []
    if isinstance(paths, (str, Path)):
        return [str(paths)]
    return [str(path) for path in paths]


def _read_yoda(paths: InputSpec) -> List[dict]:
    return [yoda.read(path) for path in _normalize_paths(paths)]


def _analysis_component_matches(component: str, analysis: str) -> bool:
    if component == analysis:
        return True
    if not component.startswith(analysis) or len(component) == len(analysis):
        return False
    return component[len(analysis)] in ":;[("


def _path_matches_analysis_label(path: str, analysis: str, label: str) -> bool:
    parts = [part for part in path.split("/") if part]
    if not parts:
        return False
    if parts[0] == "RAW":
        if len(parts) < 3:
            return False
        return _analysis_component_matches(parts[1], analysis) and parts[2] == label
    if len(parts) < 2:
        return False
    return _analysis_component_matches(parts[0], analysis) and parts[1] == label


def _resolve_input_path(aos: dict, analysis: str, label: str) -> str:
    candidates = (
        f"/{analysis}/{label}",
        f"/RAW/{analysis}/{label}",
    )
    for candidate in candidates:
        if candidate in aos:
            return candidate
    option_candidates = sorted(path for path in aos.keys() if _path_matches_analysis_label(path, analysis, label))
    if option_candidates:
        return option_candidates[0]
    available = sorted(path for path in aos.keys() if analysis in path)
    preview = ", ".join(available[:10])
    raise KeyError(f"Could not find {label} for analysis {analysis}. Tried {candidates}. Available examples: {preview}")


def _resolve_input_sources(aos_list: Sequence[dict], analysis: str, label: str) -> List[Tuple[dict, str]]:
    if not aos_list:
        return []
    resolved: List[Tuple[dict, str]] = []
    for aos in aos_list:
        try:
            resolved.append((aos, _resolve_input_path(aos, analysis, label)))
        except KeyError:
            continue
    if resolved:
        return resolved
    raise KeyError(f"Could not find {label} for analysis {analysis} in any provided YODA input")


def _x_edges(obj):
    if hasattr(obj, "xEdges"):
        return obj.xEdges()
    raise RuntimeError(f"Object of type {type(obj).__name__} does not expose xEdges()")


def _new_estimate(edges, path: str):
    hist = new_binned_estimate(edges)
    hist.setPath(path)
    try:
        hist.setAnnotation("Legend", "Herwig 7 NLO")
    except Exception:
        pass
    return hist


def _analysis_labels(analysis: str) -> tuple[str, ...]:
    try:
        return ANALYSIS_LABELS[analysis]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported analysis '{analysis}'. Expected one of {tuple(sorted(ANALYSIS_LABELS))}."
        ) from exc


def _analysis_all_labels(analysis: str) -> tuple[str, ...]:
    return ANALYSIS_ALL_LABELS.get(analysis, ())


def _analysis_derived_moments(analysis: str) -> tuple[tuple[str, str, str], ...]:
    return ANALYSIS_DERIVED_MOMENTS.get(analysis, ())


def _write_yoda_compat(path: str, objects):
    try:
        yoda.write(objects, path)
        return
    except TypeError:
        pass

    try:
        yoda.write(path, objects)
        return
    except TypeError as exc:
        raise RuntimeError(f"Could not write YODA file using either signature: {exc}")


def write_yoda_gz(objects: Dict[str, object], output_path: str) -> None:
    path = Path(output_path)
    object_list = list(objects.values())
    if path.suffix == ".gz":
        tmp_path = path.with_suffix("")
        _write_yoda_compat(str(tmp_path), object_list)
        with open(tmp_path, "rb") as src, gzip.open(path, "wb") as dst:
            dst.writelines(src)
        tmp_path.unlink()
        return
    _write_yoda_compat(str(path), object_list)


def build_plot_scatter_objects(objects: Dict[str, object]) -> Dict[str, object]:
    out: Dict[str, object] = {}
    for path, obj in objects.items():
        scatter = yoda.Scatter2D()
        scatter.setPath(path)
        try:
            scatter.setAnnotation("Legend", "Herwig 7 NLO")
        except Exception:
            pass
        for bin_obj in obj.bins():
            xlow = float(bin_obj.xMin())
            xhigh = float(bin_obj.xMax())
            xmid = 0.5 * (xlow + xhigh)
            xerr = 0.5 * (xhigh - xlow)
            if hasattr(bin_obj, "val"):
                yval = float(bin_obj.val())
                try:
                    yerr = float(bin_obj.errAvg())
                except Exception:
                    try:
                        yerr = 0.5 * (abs(float(bin_obj.errMinus())) + abs(float(bin_obj.errPlus())))
                    except Exception:
                        yerr = 0.0
            elif hasattr(bin_obj, "height"):
                yval = float(bin_obj.height())
                yerr = float(bin_obj.heightErr()) if hasattr(bin_obj, "heightErr") else 0.0
            else:
                continue
            if not numpy.isfinite(yval):
                continue
            if not numpy.isfinite(yerr) or yerr < 0.0:
                yerr = 0.0
            point = yoda.Point2D()
            point.setX(xmid)
            point.setY(yval)
            try:
                point.setXErrs(xerr, xerr)
            except Exception:
                pass
            if hasattr(point, "setYErrs"):
                point.setYErrs(yerr, yerr)
            elif hasattr(point, "setErrs"):
                point.setErrs(yerr, yerr)
            elif hasattr(point, "setErr"):
                point.setErr(yerr)
            scatter.addPoint(point)
        if scatter.numPoints() > 0:
            out[path] = scatter
    return out


def variation_object_path(path: str, variation_name: str) -> str:
    if variation_name == "nominal":
        return path
    return f"{path}[{variation_name}]"


def split_variation_object_path(path: str) -> tuple[str, str]:
    match = re.match(r"^(?P<base>.+)\[(?P<variation>[^\[\]]+)\]$", path)
    if not match:
        return path, "nominal"
    return match.group("base"), match.group("variation")


def clone_objects_with_variation(objects: Dict[str, object], variation_name: str) -> Dict[str, object]:
    cloned: Dict[str, object] = {}
    for path, obj in objects.items():
        new_path = variation_object_path(path, variation_name)
        if new_path == path:
            cloned[path] = obj
            continue
        try:
            clone = obj.clone()
        except Exception as exc:
            raise RuntimeError(f"Could not clone YODA object at {path}: {exc}") from exc
        clone.setPath(new_path)
        cloned[new_path] = clone
    return cloned


def _scatter_has_points(scatter: object) -> tuple[bool, str]:
    points = list(scatter.points()) if hasattr(scatter, "points") else []
    if not points:
        return False, "contains no plotted points"
    for index, point in enumerate(points, start=1):
        x = float(point.x())
        y = float(point.y())
        if not numpy.isfinite(x) or not numpy.isfinite(y):
            return False, f"has non-finite coordinates in plotted point {index}"
    return True, ""


def _scatter_points_compatible(nominal_obj: object, variation_obj: object) -> tuple[bool, str]:
    nominal_points = list(nominal_obj.points()) if hasattr(nominal_obj, "points") else []
    variation_points = list(variation_obj.points()) if hasattr(variation_obj, "points") else []
    if not variation_points:
        return False, "contains no finite plotted points"
    if len(variation_points) != len(nominal_points):
        return False, f"has {len(variation_points)} plotted points, expected {len(nominal_points)}"

    for index, (nominal_point, variation_point) in enumerate(zip(nominal_points, variation_points), start=1):
        nominal_x = float(nominal_point.x())
        variation_x = float(variation_point.x())
        variation_y = float(variation_point.y())
        if not numpy.isfinite(variation_x) or not numpy.isfinite(variation_y):
            return False, f"has non-finite coordinates in plotted point {index}"
        if not numpy.isclose(variation_x, nominal_x, rtol=1e-12, atol=1e-12):
            return False, (
                f"uses a different x grid in plotted point {index} "
                f"({variation_x} vs nominal {nominal_x})"
            )
    return True, ""


def filter_compatible_variation_plot_objects(
    nominal_objects: Dict[str, object],
    variation_objects: Dict[str, object],
) -> tuple[Dict[str, object], Dict[str, str]]:
    compatible: Dict[str, object] = {}
    skipped: Dict[str, str] = {}

    for path, variation_obj in variation_objects.items():
        nominal_obj = nominal_objects.get(path)
        if nominal_obj is None:
            skipped[path] = "missing nominal plotted object"
            continue
        ok, reason = _scatter_points_compatible(nominal_obj, variation_obj)
        if ok:
            compatible[path] = variation_obj
        else:
            skipped[path] = reason

    return compatible, skipped


def sanitize_plot_scatter_objects(objects: Dict[str, object]) -> tuple[Dict[str, object], Dict[str, str]]:
    sanitized: Dict[str, object] = {}
    dropped: Dict[str, str] = {}
    nominal_objects: Dict[str, object] = {}
    variation_entries: List[tuple[str, str, object]] = []

    for path, obj in objects.items():
        scatter = obj
        if hasattr(scatter, "type") and "Scatter" not in scatter.type() and hasattr(scatter, "mkScatter"):
            scatter = scatter.mkScatter()
            try:
                scatter.setPath(path)
            except Exception:
                pass

        base_path, variation_name = split_variation_object_path(path)
        if variation_name == "nominal":
            ok, reason = _scatter_has_points(scatter)
            if ok:
                nominal_objects[base_path] = scatter
                sanitized[path] = scatter
            else:
                dropped[path] = reason
            continue

        variation_entries.append((path, base_path, scatter))

    for path, base_path, scatter in variation_entries:
        nominal_obj = nominal_objects.get(base_path)
        if nominal_obj is None:
            dropped[path] = "missing non-empty nominal plotted object"
            continue
        ok, reason = _scatter_points_compatible(nominal_obj, scatter)
        if ok:
            sanitized[path] = scatter
        else:
            dropped[path] = reason

    return sanitized, dropped


def sanitize_plot_yoda(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    drop_band_annotations: bool = False,
) -> tuple[int, int, Dict[str, str]]:
    objects = yoda.read(str(input_path))
    sanitized, dropped = sanitize_plot_scatter_objects(objects)
    if drop_band_annotations:
        for obj in sanitized.values():
            for key in ("ErrorBandColor", "ErrorBandOpacity"):
                try:
                    obj.rmAnnotation(key)
                except Exception:
                    pass
    write_yoda_gz(sanitized, str(output_path))
    return len(sanitized), len(dropped), dropped


def _scatter_bin_key(point: object) -> tuple[float, float]:
    return (round(float(point.xMin()), 12), round(float(point.xMax()), 12))


def _copy_scatter_annotations(source: object, target: object) -> None:
    try:
        annotations = list(source.annotations())
    except Exception:
        annotations = []
    for key in annotations:
        if key == "Path":
            continue
        try:
            target.setAnnotation(key, source.annotation(key))
        except Exception:
            continue


def _extract_err_pair(errs: object) -> tuple[float, float]:
    if hasattr(errs, "minus") and hasattr(errs, "plus"):
        return float(errs.minus), float(errs.plus)
    if isinstance(errs, (tuple, list)) and len(errs) >= 2:
        return float(errs[0]), float(errs[1])
    value = float(errs)
    return value, value


def _filtered_scatter(scatter: object, keep_bins: set[tuple[float, float]]) -> object:
    filtered = yoda.Scatter2D()
    try:
        filtered.setPath(scatter.path())
    except Exception:
        pass
    _copy_scatter_annotations(scatter, filtered)

    for point in scatter.points():
        if _scatter_bin_key(point) not in keep_bins:
            continue
        new_point = yoda.Point2D()
        new_point.setX(float(point.x()))
        new_point.setY(float(point.y()))
        xerrs = point.xErrs()
        yerrs = point.yErrs()
        xerr_dn, xerr_up = _extract_err_pair(xerrs)
        yerr_dn, yerr_up = _extract_err_pair(yerrs)
        try:
            new_point.setXErrs(xerr_dn, xerr_up)
        except Exception:
            pass
        _set_point_y_errs(new_point, yerr_dn, yerr_up)
        filtered.addPoint(new_point)
    return filtered


def harmonize_plot_yoda_bin_grids(yoda_paths: Sequence[Union[str, Path]]) -> tuple[int, int]:
    input_paths = [Path(path) for path in yoda_paths]
    if len(input_paths) < 2:
        return 0, 0

    loaded = [(path, yoda.read(str(path))) for path in input_paths]
    common_object_paths = set.intersection(*(set(objects.keys()) for _, objects in loaded))

    updated_objects = 0
    updated_files = set()

    for object_path in sorted(common_object_paths):
        scatters = []
        scatter_maps = []
        for file_path, objects in loaded:
            obj = objects.get(object_path)
            if obj is None:
                scatters = []
                break
            scatter = obj
            if hasattr(scatter, "type") and "Scatter" not in scatter.type() and hasattr(scatter, "mkScatter"):
                scatter = scatter.mkScatter()
                try:
                    scatter.setPath(object_path)
                except Exception:
                    pass
            if not hasattr(scatter, "points"):
                scatters = []
                break
            scatters.append(scatter)
            scatter_maps.append((file_path, objects))

        if len(scatters) != len(loaded):
            continue

        bin_sets = [{_scatter_bin_key(point) for point in scatter.points()} for scatter in scatters]
        common_bins = set.intersection(*bin_sets)
        if not common_bins:
            continue
        if all(bin_set == common_bins for bin_set in bin_sets):
            continue

        for (file_path, objects), scatter in zip(scatter_maps, scatters):
            filtered = _filtered_scatter(scatter, common_bins)
            objects[object_path] = filtered
            updated_files.add(file_path)
            updated_objects += 1

    for file_path, objects in loaded:
        if file_path in updated_files:
            write_yoda_gz(objects, str(file_path))

    return len(updated_files), updated_objects


def _set_point_y_errs(point: object, err_dn: float, err_up: float) -> None:
    try:
        point.setYErrs(err_dn, err_up)
        return
    except Exception:
        pass
    try:
        point.setYErrs((err_dn, err_up))
        return
    except Exception:
        pass
    try:
        point.setErrs(err_dn, err_up)
        return
    except Exception:
        pass
    try:
        point.setErrs((err_dn, err_up))
        return
    except Exception:
        pass
    if hasattr(point, "setErr"):
        point.setErr(max(err_dn, err_up))


def build_scale_envelope_plot_objects(
    objects: Dict[str, object],
    variation_names: Sequence[str] = ("ScaleFactorDown", "ScaleFactorUp"),
) -> tuple[Dict[str, object], Dict[str, str]]:
    sanitized, dropped = sanitize_plot_scatter_objects(objects)
    grouped: Dict[str, Dict[str, object]] = {}
    for path, obj in sanitized.items():
        base_path, variation_name = split_variation_object_path(path)
        grouped.setdefault(base_path, {})[variation_name] = obj

    envelope_objects: Dict[str, object] = {}
    for base_path, group in grouped.items():
        nominal = group.get("nominal")
        if nominal is None:
            dropped[base_path] = "missing nominal plotted object after sanitization"
            continue

        try:
            envelope = nominal.clone()
        except Exception as exc:
            dropped[base_path] = f"could not clone nominal plotted object: {exc}"
            continue
        envelope.setPath(base_path)

        available_variations = [group[name] for name in variation_names if name in group]
        for index, nominal_point in enumerate(nominal.points()):
            nominal_y = float(nominal_point.y())
            values = [nominal_y]
            for variation in available_variations:
                values.append(float(variation.point(index).y()))
            ymin = min(values)
            ymax = max(values)
            err_dn = nominal_y - ymin
            err_up = ymax - nominal_y

            # Older YODA plotting code may collapse exactly symmetric errors to scalars.
            if numpy.isclose(err_dn, err_up, rtol=1e-15, atol=1e-15):
                tweak = max(abs(nominal_y) * 1e-15, 1e-15)
                err_up += tweak

            _set_point_y_errs(envelope.point(index), err_dn, err_up)

        envelope_objects[base_path] = envelope

    return envelope_objects, dropped


def build_scale_envelope_plot_yoda(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    variation_names: Sequence[str] = ("ScaleFactorDown", "ScaleFactorUp"),
) -> tuple[int, int, Dict[str, str]]:
    objects = yoda.read(str(input_path))
    envelope_objects, dropped = build_scale_envelope_plot_objects(objects, variation_names=variation_names)
    write_yoda_gz(envelope_objects, str(output_path))
    return len(envelope_objects), len(dropped), dropped


def _single_bin_value(aos: dict, path: str, index: int) -> tuple[float, float]:
    hist = aos[path]
    bin_obj = hist.bin(index)
    if hasattr(bin_obj, "val"):
        value = float(bin_obj.val())
        if hasattr(bin_obj, "errAvg"):
            error = 0.0
            try:
                error = float(bin_obj.errAvg("stats"))
            except Exception:
                try:
                    error = float(bin_obj.errAvg())
                except Exception:
                    try:
                        error = 0.5 * (abs(float(bin_obj.errMinus())) + abs(float(bin_obj.errPlus())))
                    except Exception:
                        error = 0.0
        else:
            error = 0.0
        return value, error
    if hasattr(bin_obj, "height"):
        value = float(bin_obj.height())
        error = float(bin_obj.heightErr()) if hasattr(bin_obj, "heightErr") else 0.0
        return value, error
    raise RuntimeError(f"Unsupported YODA bin type: {type(bin_obj).__name__}")


def _combine_measurements(values: Sequence[Tuple[float, float]]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    nshards = len(values)
    mean = float(sum(value for value, _ in values) / nshards)
    error = float(numpy.sqrt(sum(error * error for _, error in values)) / nshards)
    return mean, error


def _bin_value(aos_sources: Sequence[Tuple[dict, str]], index: int) -> tuple[float, float]:
    if not aos_sources:
        return 0.0, 0.0
    values = [_single_bin_value(aos, path, index) for aos, path in aos_sources]
    return _combine_measurements(values)


def _bin_value_with_subtraction(
    aos_pos: Sequence[Tuple[dict, str]],
    aos_neg: Sequence[Tuple[dict, str]],
    index: int,
) -> tuple[float, float]:
    pos_value, pos_error = _bin_value(aos_pos, index)
    if not aos_neg:
        return pos_value, pos_error
    neg_value, neg_error = _bin_value(aos_neg, index)
    return pos_value - neg_value, float(numpy.sqrt(pos_error**2 + neg_error**2))


def _ratio_value_error(numerator: float, numerator_error: float,
                       denominator: float, denominator_error: float) -> tuple[float, float]:
    if (not numpy.isfinite(numerator) or not numpy.isfinite(numerator_error) or
            not numpy.isfinite(denominator) or not numpy.isfinite(denominator_error) or
            abs(denominator) <= 0.0):
        return float("nan"), 0.0
    value = numerator / denominator
    error = float(numpy.sqrt((numerator_error / denominator)**2 +
                             (numerator * denominator_error / (denominator**2))**2))
    return value, error


def _combine_setup_bin(
    setup: str,
    index: int,
    zero_sources: Sequence[Tuple[dict, str]],
    pp_sources: Sequence[Tuple[dict, str]],
    pm_sources: Sequence[Tuple[dict, str]],
    mp_sources: Sequence[Tuple[dict, str]],
    mm_sources: Sequence[Tuple[dict, str]],
    zero_sub_sources: Sequence[Tuple[dict, str]],
    pp_sub_sources: Sequence[Tuple[dict, str]],
    pm_sub_sources: Sequence[Tuple[dict, str]],
    mp_sub_sources: Sequence[Tuple[dict, str]],
    mm_sub_sources: Sequence[Tuple[dict, str]],
) -> tuple[float, float, float, float]:
    n00, e00 = _bin_value_with_subtraction(zero_sources, zero_sub_sources, index)
    npp, epp = _bin_value_with_subtraction(pp_sources, pp_sub_sources, index)
    npm, epm = _bin_value_with_subtraction(pm_sources, pm_sub_sources, index)

    if setup == "GAMMA":
        sigma_value = n00 if zero_sources else 0.5 * (npp + npm)
        sigma_error = e00 if zero_sources else 0.5 * numpy.sqrt(epp**2 + epm**2)
        delta_value = 0.5 * (npp - npm)
        delta_error = 0.5 * numpy.sqrt(epp**2 + epm**2)
    else:
        nmp, emp = _bin_value_with_subtraction(mp_sources, mp_sub_sources, index)
        nmm, emm = _bin_value_with_subtraction(mm_sources, mm_sub_sources, index)
        sigma_value = n00 if zero_sources else 0.25 * (npp + npm + nmp + nmm)
        sigma_error = e00 if zero_sources else 0.25 * numpy.sqrt(epp**2 + epm**2 + emp**2 + emm**2)
        delta_value = 0.25 * (npp + nmm - npm - nmp)
        delta_error = 0.25 * numpy.sqrt(epp**2 + epm**2 + emp**2 + emm**2)
    return sigma_value, sigma_error, delta_value, delta_error


def build_dis_polarized_objects(
    setup: str,
    pp_path: InputSpec,
    pm_path: InputSpec,
    zero_path: InputSpec = None,
    mp_path: InputSpec = None,
    mm_path: InputSpec = None,
    pp_subtract_path: InputSpec = None,
    pm_subtract_path: InputSpec = None,
    zero_subtract_path: InputSpec = None,
    mp_subtract_path: InputSpec = None,
    mm_subtract_path: InputSpec = None,
    analysis: str = "MC_DIS_BREIT",
) -> Dict[str, object]:
    setup = setup.upper()
    if setup not in SETUPS:
        raise ValueError(f"Unsupported setup '{setup}'. Expected one of {SETUPS}.")

    aos_pp = _read_yoda(pp_path)
    aos_pm = _read_yoda(pm_path)
    aos_mp = _read_yoda(mp_path)
    aos_mm = _read_yoda(mm_path)
    aos_zero = _read_yoda(zero_path)
    aos_pp_sub = _read_yoda(pp_subtract_path)
    aos_pm_sub = _read_yoda(pm_subtract_path)
    aos_mp_sub = _read_yoda(mp_subtract_path)
    aos_mm_sub = _read_yoda(mm_subtract_path)
    aos_zero_sub = _read_yoda(zero_subtract_path)

    if not aos_pp or not aos_pm:
        raise ValueError("PP and PM YODA inputs are required")
    if setup in FIVE_HELICITY_SETUPS and (not aos_mp or not aos_mm):
        raise ValueError(f"Setup {setup} requires MP and MM YODA inputs")

    output: Dict[str, object] = {}

    labels = _analysis_labels(analysis)
    all_labels = set(_analysis_all_labels(analysis))
    derived_moments = _analysis_derived_moments(analysis)

    for label in labels:
        pp_sources = _resolve_input_sources(aos_pp, analysis, label)
        pm_sources = _resolve_input_sources(aos_pm, analysis, label)
        mp_sources = _resolve_input_sources(aos_mp, analysis, label)
        mm_sources = _resolve_input_sources(aos_mm, analysis, label)
        zero_sources = _resolve_input_sources(aos_zero, analysis, label)
        pp_sub_sources = _resolve_input_sources(aos_pp_sub, analysis, label)
        pm_sub_sources = _resolve_input_sources(aos_pm_sub, analysis, label)
        mp_sub_sources = _resolve_input_sources(aos_mp_sub, analysis, label)
        mm_sub_sources = _resolve_input_sources(aos_mm_sub, analysis, label)
        zero_sub_sources = _resolve_input_sources(aos_zero_sub, analysis, label)
        edges = _x_edges(pp_sources[0][0][pp_sources[0][1]])
        sigma_path = f"/{analysis}/{label}"
        delta_path = f"/{analysis}/D{label}"
        sigma_hist = _new_estimate(edges, sigma_path)
        delta_hist = _new_estimate(edges, delta_path)
        output[sigma_path] = sigma_hist
        output[delta_path] = delta_hist
        all_hist = None
        if label in all_labels:
            all_path = f"/{analysis}/ALL{label}"
            all_hist = _new_estimate(edges, all_path)
            output[all_path] = all_hist

        for index in range(1, sigma_hist.numBins() + 1):
            sigma_value, sigma_error, delta_value, delta_error = _combine_setup_bin(
                setup,
                index,
                zero_sources,
                pp_sources,
                pm_sources,
                mp_sources,
                mm_sources,
                zero_sub_sources,
                pp_sub_sources,
                pm_sub_sources,
                mp_sub_sources,
                mm_sub_sources,
            )
            set_bin_val_err(sigma_hist.bin(index), sigma_value, sigma_error)
            set_bin_val_err(delta_hist.bin(index), delta_value, delta_error)
            if all_hist is not None:
                all_value, all_error = _ratio_value_error(delta_value, delta_error, sigma_value, sigma_error)
                set_bin_val_err(all_hist.bin(index), all_value, all_error)

    for moment_label, numerator_label, denominator_label in derived_moments:
        num_pp_sources = _resolve_input_sources(aos_pp, analysis, numerator_label)
        num_pm_sources = _resolve_input_sources(aos_pm, analysis, numerator_label)
        num_mp_sources = _resolve_input_sources(aos_mp, analysis, numerator_label)
        num_mm_sources = _resolve_input_sources(aos_mm, analysis, numerator_label)
        num_zero_sources = _resolve_input_sources(aos_zero, analysis, numerator_label)
        num_pp_sub_sources = _resolve_input_sources(aos_pp_sub, analysis, numerator_label)
        num_pm_sub_sources = _resolve_input_sources(aos_pm_sub, analysis, numerator_label)
        num_mp_sub_sources = _resolve_input_sources(aos_mp_sub, analysis, numerator_label)
        num_mm_sub_sources = _resolve_input_sources(aos_mm_sub, analysis, numerator_label)
        num_zero_sub_sources = _resolve_input_sources(aos_zero_sub, analysis, numerator_label)

        den_pp_sources = _resolve_input_sources(aos_pp, analysis, denominator_label)
        den_pm_sources = _resolve_input_sources(aos_pm, analysis, denominator_label)
        den_mp_sources = _resolve_input_sources(aos_mp, analysis, denominator_label)
        den_mm_sources = _resolve_input_sources(aos_mm, analysis, denominator_label)
        den_zero_sources = _resolve_input_sources(aos_zero, analysis, denominator_label)
        den_pp_sub_sources = _resolve_input_sources(aos_pp_sub, analysis, denominator_label)
        den_pm_sub_sources = _resolve_input_sources(aos_pm_sub, analysis, denominator_label)
        den_mp_sub_sources = _resolve_input_sources(aos_mp_sub, analysis, denominator_label)
        den_mm_sub_sources = _resolve_input_sources(aos_mm_sub, analysis, denominator_label)
        den_zero_sub_sources = _resolve_input_sources(aos_zero_sub, analysis, denominator_label)

        edges = _x_edges(num_pp_sources[0][0][num_pp_sources[0][1]])
        moment_path = f"/{analysis}/{moment_label}"
        moment_hist = _new_estimate(edges, moment_path)
        output[moment_path] = moment_hist

        for index in range(1, moment_hist.numBins() + 1):
            sigma_num, sigma_num_error, _, _ = _combine_setup_bin(
                setup,
                index,
                num_zero_sources,
                num_pp_sources,
                num_pm_sources,
                num_mp_sources,
                num_mm_sources,
                num_zero_sub_sources,
                num_pp_sub_sources,
                num_pm_sub_sources,
                num_mp_sub_sources,
                num_mm_sub_sources,
            )
            sigma_den, sigma_den_error, _, _ = _combine_setup_bin(
                setup,
                index,
                den_zero_sources,
                den_pp_sources,
                den_pm_sources,
                den_mp_sources,
                den_mm_sources,
                den_zero_sub_sources,
                den_pp_sub_sources,
                den_pm_sub_sources,
                den_mp_sub_sources,
                den_mm_sub_sources,
            )
            moment_value, moment_error = _ratio_value_error(
                sigma_num, sigma_num_error, sigma_den, sigma_den_error
            )
            set_bin_val_err(moment_hist.bin(index), moment_value, moment_error)

    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--setup", choices=SETUPS, default="ALL", help="DIS setup to analyze.")
    parser.add_argument("--00", dest="zero", help="YODA file with 00 helicity run.")
    parser.add_argument("--PP", required=True, help="YODA file with ++ helicity run.")
    parser.add_argument("--PM", required=True, help="YODA file with +- helicity run.")
    parser.add_argument("--MP", help="YODA file with -+ helicity run.")
    parser.add_argument("--MM", help="YODA file with -- helicity run.")
    parser.add_argument("--00-sub", dest="zero_sub", help="YODA file to subtract from 00.")
    parser.add_argument("--PP-sub", help="YODA file to subtract from ++.")
    parser.add_argument("--PM-sub", help="YODA file to subtract from +-.")
    parser.add_argument("--MP-sub", help="YODA file to subtract from -+.")
    parser.add_argument("--MM-sub", help="YODA file to subtract from --.")
    parser.add_argument("--analysis", default="MC_DIS_BREIT", help="Analysis path prefix.")
    parser.add_argument("output", help="Name of the output YODA file.")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    objects = build_dis_polarized_objects(
        setup=args.setup,
        zero_path=args.zero,
        pp_path=args.PP,
        pm_path=args.PM,
        mp_path=args.MP,
        mm_path=args.MM,
        zero_subtract_path=args.zero_sub,
        pp_subtract_path=args.PP_sub,
        pm_subtract_path=args.PM_sub,
        mp_subtract_path=args.MP_sub,
        mm_subtract_path=args.MM_sub,
        analysis=args.analysis,
    )
    write_yoda_gz(objects, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
