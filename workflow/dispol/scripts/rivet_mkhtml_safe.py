#!/usr/bin/env python3
"""Run rivet-mkhtml with a small compatibility patch for older YODA plotting code."""

from __future__ import annotations

import os
import runpy
import sys
import traceback
from collections.abc import Iterable

import numpy as np


def _coerce_plot_values(vals):
    if isinstance(vals, (str, bytes)):
        return [vals]
    if isinstance(vals, Iterable):
        return vals
    return [vals]


def _annotation_list_literal(vals) -> str:
    return repr(list(_coerce_plot_values(vals)))


def _annotation_color_text(vals) -> str:
    coerced = list(_coerce_plot_values(vals))
    return " ".join(str(val) for val in coerced)


def _normalize_color_sequence(colors):
    if colors is None:
        return colors
    if isinstance(colors, (str, bytes)):
        return [_annotation_color_text(colors)]
    return [_annotation_color_text(color) for color in colors]


def _coerce_boollike(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, (str, bytes)):
        text = str(value).strip().lower()
        if text in {"1", "true", "yes", "on"}:
            return True
        if text in {"0", "false", "no", "off", ""}:
            return False
    return bool(value)


def _normalize_hist_yaml_tree(node, force_no_error_bands: bool = False):
    if isinstance(node, dict):
        normalized = {}
        for key, value in node.items():
            if force_no_error_bands and (key == "BandUncertainty" or str(key).startswith("multiweight")):
                continue
            normalized[key] = _normalize_hist_yaml_tree(value, force_no_error_bands=force_no_error_bands)
        if "histograms" in normalized and isinstance(normalized["histograms"], dict):
            for hist_settings in normalized["histograms"].values():
                if not isinstance(hist_settings, dict):
                    continue
                for key, value in list(hist_settings.items()):
                    if key in {"ErrorBands", "ErrorBars"}:
                        hist_settings[key] = int(_coerce_boollike(value))
                    if key.endswith("ErrorBandOpacity"):
                        hist_settings[key] = _annotation_list_literal(value)
                    if key.endswith("ErrorBandColor"):
                        hist_settings[key] = _annotation_color_text(value)
                if force_no_error_bands:
                    hist_settings["ErrorBands"] = 0
                    hist_settings["ErrorBars"] = 1
                    hist_settings.pop("BandUncertainty", None)
                    for key in list(hist_settings.keys()):
                        if str(key).startswith("multiweight"):
                            del hist_settings[key]
        return normalized
    if isinstance(node, list):
        return [_normalize_hist_yaml_tree(item, force_no_error_bands=force_no_error_bands) for item in node]
    return node


def _patch_yoda_write_lists() -> None:
    import yoda.plotting.fetch_data as fetch_data

    original_write_lists = fetch_data.writeLists

    def _coerce_array_like(vals):
        if vals is None:
            return vals
        if isinstance(vals, np.ndarray):
            return vals
        return np.atleast_1d(vals)

    def _sanitize_plot_values(value):
        if not isinstance(value, dict):
            return _coerce_array_like(value)
        return {
            label: _coerce_array_like(vals)
            for label, vals in value.items()
        }

    def write_lists(cmd, val_type, *value_dicts):
        """Delegate to upstream writeLists after coercing scalar values."""
        return original_write_lists(
            cmd,
            val_type,
            *(_sanitize_plot_values(value_dict) for value_dict in value_dicts),
        )

    fetch_data.writeLists = write_lists


def _patch_fetch_data_band_annotations() -> None:
    import yoda.plotting.fetch_data as fetch_data

    original_mkcurves1d = fetch_data.mkCurves1D
    force_no_error_bands = _coerce_boollike(os.environ.get("DISPOL_FORCE_NO_ERROR_BANDS"))

    def _normalize_annotation_value(ao, key):
        try:
            value = ao.annotation(key, None)
        except Exception:
            return
        if value is None:
            return
        if key.endswith("ErrorBandOpacity"):
            ao.setAnnotation(key, _annotation_list_literal(value))
        elif key.endswith("ErrorBandColor"):
            ao.setAnnotation(key, _annotation_color_text(value))

    def _normalize_ao(ao):
        try:
            keys = list(ao.annotations())
        except Exception:
            keys = []
        for key in keys:
            if key.endswith("ErrorBandOpacity") or key.endswith("ErrorBandColor"):
                _normalize_annotation_value(ao, key)

    def mkcurves1d_with_safe_annotations(
        aos,
        ratio_panels,
        error_bars=True,
        variation_dict=None,
        colors=None,
        deviation=False,
        line_styles=["-", "--", "-.", ":"],
        **kwargs,
    ):
        try:
            ao_values = aos.values()
        except Exception:
            ao_values = aos
        for ao in ao_values:
            _normalize_ao(ao)
        if isinstance(variation_dict, dict):
            for varset in variation_dict.values():
                if isinstance(varset, dict):
                    for ao in varset.values():
                        _normalize_ao(ao)
        if force_no_error_bands:
            variation_dict = {}
        colors = _normalize_color_sequence(colors)
        return original_mkcurves1d(
            aos,
            ratio_panels,
            error_bars=error_bars,
            variation_dict=variation_dict,
            colors=colors,
            deviation=deviation,
            line_styles=line_styles,
            **kwargs,
        )

    fetch_data.mkCurves1D = mkcurves1d_with_safe_annotations


def _patch_script_generator_traceback() -> None:
    import yoda.plotting.script_generator as script_generator

    original_process = script_generator.process
    force_no_error_bands = _coerce_boollike(os.environ.get("DISPOL_FORCE_NO_ERROR_BANDS"))

    def normalize_band_settings(yaml_file):
        return _normalize_hist_yaml_tree(yaml_file, force_no_error_bands=force_no_error_bands)

    def process_with_traceback(*args, **kwargs):
        try:
            if args:
                args = (normalize_band_settings(args[0]), *args[1:])
            if "yaml_file" in kwargs:
                kwargs["yaml_file"] = normalize_band_settings(kwargs["yaml_file"])
            return original_process(*args, **kwargs)
        except Exception:
            traceback.print_exc()
            raise

    script_generator.process = process_with_traceback


def _patch_yoda_parseyaml() -> None:
    import yoda.util as yoda_util

    original_parseyaml = yoda_util._parseyaml
    force_no_error_bands = _coerce_boollike(os.environ.get("DISPOL_FORCE_NO_ERROR_BANDS"))

    def parseyaml_with_sanitized_bands(*args, **kwargs):
        parsed = original_parseyaml(*args, **kwargs)
        return _normalize_hist_yaml_tree(parsed, force_no_error_bands=force_no_error_bands)

    yoda_util._parseyaml = parseyaml_with_sanitized_bands


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        raise SystemExit("Usage: rivet_mkhtml_safe.py <path-to-rivet-mkhtml> [args...]")

    tool = args[0]
    sys.argv = [tool, *args[1:]]
    _patch_yoda_parseyaml()
    _patch_yoda_write_lists()
    _patch_fetch_data_band_annotations()
    _patch_script_generator_traceback()
    runpy.run_path(tool, run_name="__main__")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
