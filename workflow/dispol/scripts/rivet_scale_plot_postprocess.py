#!/usr/bin/env python3.10
"""Patch generated Rivet plot scripts for DIS scale-envelope rendering."""

from __future__ import annotations

import re
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Tuple


PATCH_SENTINEL = "# Codex scale-envelope patch"
FIXED_ORDER_SCALE_COLOR = "#009E73"
LOWER_CENTER_LEGEND_OBJECTS = {
    "/MC_DIS_BREIT/Zeta",
    "/MC_DIS_BREIT/DZeta",
    "/MC_DIS_PS/Zeta",
    "/MC_DIS_PS/DZeta",
}


_HELPER_BLOCK = textwrap.dedent(
    """

    # Codex scale-envelope patch
    desired_reference_label = '__REFERENCE_LABEL__'

    def _label_text(label):
        return str(label).lower()

    def _is_raw_overlay_label(label):
        text = _label_text(label)
        return '.raw_powheg_' in text or text.startswith('raw_powheg')

    def _is_poldis_reference_label(label):
        text = _label_text(label)
        return (
            str(label) == desired_reference_label
            or text == 'data'
            or text.startswith('poldis')
            or 'poldis' in text
            or (
                ('reference' in text or 'ref' in text)
                and 'herwig' not in text
                and 'standalone' not in text
            )
        )

    def _is_scale_band_label(label):
        text = _label_text(label)
        return (
            not _is_raw_overlay_label(label)
            and not _is_poldis_reference_label(label)
            and (
                '.herwig_' in text
                or 'herwig' in text
                or 'standalone' in text
                or 'scale_envelope' in text
            )
        )

    raw_overlay_labels = {label for label in dataf['yvals'] if _is_raw_overlay_label(label)}
    poldis_reference_label = next(
        (label for label in dataf.get('add_legend_handle', []) if _is_poldis_reference_label(label)),
        None,
    )
    if poldis_reference_label is None:
        poldis_reference_label = next((label for label in dataf['yvals'] if _is_poldis_reference_label(label)), None)

    if poldis_reference_label is None:
        scale_band_labels = {label for label in dataf['yvals'] if label not in raw_overlay_labels}
        reference_label = next(
            (
                label
                for label in dataf.get('add_legend_handle', [])
                if 'NOSPIN-UNPOL' in str(label).upper()
            ),
            None,
        )
        if reference_label is None:
            reference_label = next((label for label in dataf['yvals'] if 'NOSPIN-UNPOL' in str(label).upper()), None)
    else:
        scale_band_labels = {label for label in dataf['yvals'] if _is_scale_band_label(label)}
        reference_label = poldis_reference_label
        if not scale_band_labels:
            scale_band_labels = {
                label
                for label in dataf['yvals']
                if label != reference_label and label not in raw_overlay_labels
            }
    original_reference_label = reference_label
    scale_band_opacity = 0.5
    fixed_order_scale_color = '__SCALE_BAND_COLOR__'

    def _rename_plot_label(old, new):
        if old is None or old == new or old not in dataf['yvals']:
            return old
        for key in ('yvals', 'yerrs', 'xpoints', 'xedges', 'xerrs'):
            mapping = dataf.get(key)
            if isinstance(mapping, dict) and old in mapping:
                mapping[new] = mapping.pop(old)
        variation_mapping = dataf.get('variation_yvals')
        if isinstance(variation_mapping, dict):
            renamed_variations = {}
            for label in list(variation_mapping):
                if label.startswith(old):
                    renamed_variations[new + label[len(old):]] = variation_mapping.pop(label)
            variation_mapping.update(renamed_variations)
        if isinstance(dataf.get('add_legend_handle'), list):
            dataf['add_legend_handle'] = [new if label == old else label for label in dataf['add_legend_handle']]
        style_map = globals().get('styles')
        if isinstance(style_map, dict) and old in style_map:
            style_map[new] = style_map.pop(old)
        return new

    if poldis_reference_label is not None:
        reference_label = _rename_plot_label(reference_label, desired_reference_label)

    def _copy_scale_patch_style(source_label, target_label):
        if source_label is None or target_label is None:
            return
        if target_label in styles:
            return
        if source_label in styles:
            styles[target_label] = dict(styles[source_label])
            return
        for candidate in ('Data', 'POLDIS NLO', 'POLDIS LO', 'POLDIS NNLO'):
            if candidate in styles:
                styles[target_label] = dict(styles[candidate])
                return
        for candidate, style in styles.items():
            if candidate not in scale_band_labels:
                styles[target_label] = dict(style)
                return

    def _apply_poldis_reference_style(label):
        if label is None or label not in styles:
            return
        styles[label].update({
            'color': 'black',
            'linestyle': '-',
            'lineopacity': 1.0,
            'linewidth': 1,
            'marker': 'o',
            'markersize': 2,
            'capsize': 0.0,
            'zorder': 5,
            'histstyle': 0,
            'drawstyle': 'steps-pre',
            'xerrorbars': 1,
            'yerrorbars': 1,
            'fillcolor': None,
            'fillopacity': 1.0,
            'ratio0_linestyle': '-',
            'ratio0_lineopacity': 1.0,
            'ratio0_linewidth': 1,
            'ratio0_marker': 'o',
            'ratio0_markersize': 2,
            'ratio0_capsize': 0.0,
            'ratio0_zorder': 5,
            'ratio0_histstyle': 0,
            'ratio0_drawstyle': 'steps-pre',
            'ratio0_xerrorbars': 1,
            'ratio0_yerrorbars': 0,
            'ratio0_1sigmabandcolor': 'black',
            'ratio0_1sigmabandstyle': None,
            'ratio0_1sigmabandopacity': 0.2,
        })

    def _apply_scale_band_style(label):
        if label is None or label not in styles:
            return
        styles[label].update({
            'color': fixed_order_scale_color,
            'linestyle': '-',
            'lineopacity': 1.0,
            'linewidth': 1,
            'marker': 'none',
            'markersize': 2,
            'capsize': 0.0,
            'zorder': 6,
            'histstyle': 1,
            'drawstyle': 'steps-pre',
            'xerrorbars': 1,
            'yerrorbars': 1,
            'fillcolor': None,
            'fillopacity': 1.0,
            'ratio0_linestyle': '-',
            'ratio0_lineopacity': 1.0,
            'ratio0_linewidth': 1,
            'ratio0_marker': 'none',
            'ratio0_markersize': 2,
            'ratio0_capsize': 0.0,
            'ratio0_zorder': 6,
            'ratio0_histstyle': 1,
            'ratio0_drawstyle': 'steps-pre',
            'ratio0_xerrorbars': 1,
            'ratio0_yerrorbars': 0,
            'ratio0_1sigmabandcolor': fixed_order_scale_color,
            'ratio0_1sigmabandstyle': None,
            'ratio0_1sigmabandopacity': 0.2,
        })

    def _ensure_scale_patch_styles():
        if reference_label is not None:
            _copy_scale_patch_style(original_reference_label, reference_label)
        for _label in list(dataf.get('yvals', {})):
            _copy_scale_patch_style(reference_label, _label)
        if poldis_reference_label is not None:
            _apply_poldis_reference_style(reference_label)
            for _label in scale_band_labels:
                _apply_scale_band_style(_label)

    def _symmetrize_reference_yerrs(label):
        if label is None or label not in dataf.get('yerrs', {}):
            return
        yerrs = np.asarray(dataf['yerrs'][label], dtype=float)
        if yerrs.shape[0] != 2:
            return
        down = np.asarray(yerrs[0], dtype=float)
        up = np.asarray(yerrs[1], dtype=float)
        if np.any(down > 0.0) and np.allclose(up, 0.0):
            dataf['yerrs'][label] = np.vstack([down, down])
            return
        if np.any(up > 0.0) and np.allclose(down, 0.0):
            dataf['yerrs'][label] = np.vstack([up, up])

    _symmetrize_reference_yerrs(reference_label)

    def _signed_plot_limits():
        lower = []
        upper = []
        saw_negative = False
        for label, yvals_raw in dataf.get('yvals', {}).items():
            if label in raw_overlay_labels:
                continue
            yvals = np.asarray(yvals_raw, dtype=float)
            finite = yvals[np.isfinite(yvals)]
            if finite.size == 0:
                continue
            if np.any(finite < 0.0):
                saw_negative = True
            if label in dataf.get('yerrs', {}):
                yerrs = np.asarray(dataf['yerrs'][label], dtype=float)
                if yerrs.ndim == 2 and yerrs.shape[0] == 2:
                    low = yvals - yerrs[0]
                    high = yvals + yerrs[1]
                else:
                    low = yvals
                    high = yvals
            else:
                low = yvals
                high = yvals
            low = low[np.isfinite(low)]
            high = high[np.isfinite(high)]
            if low.size:
                lower.extend(low.tolist())
            if high.size:
                upper.extend(high.tolist())
        if not lower or not upper or not saw_negative:
            return
        ymin = min(lower)
        ymax = max(upper)
        if not np.isfinite(ymin) or not np.isfinite(ymax):
            return
        if np.isclose(ymin, ymax):
            span = max(abs(ymin), 1.0)
            pad = 0.1 * span
        else:
            pad = 0.08 * (ymax - ymin)
        globals()['ax_yScale'] = 'linear'
        globals()['yLims'] = (ymin - pad, ymax + pad)

    _signed_plot_limits()

    def _safe_ratio(num, den):
        num = np.asarray(num, dtype=float)
        den = np.asarray(den, dtype=float)
        return np.divide(num, den, out=np.full_like(num, np.nan, dtype=float), where=den != 0)

    def _scale_band_color(label):
        return styles.get(label, {}).get('color', fixed_order_scale_color)

    def _foreground_zorder():
        zorders = []
        for style in styles.values():
            try:
                zorders.append(float(style.get('zorder', 0.0)))
            except Exception:
                continue
        return (max(zorders) if zorders else 0.0) + 2.0

    def _enabled_errorbar_scale(style_key, axis_key):
        style = styles.get(style_key, {})
        try:
            value = float(style.get(axis_key, 1.0))
        except Exception:
            value = 1.0
        return value if value > 0.0 else 1.0

    def _main_band_arrays(label):
        if label not in scale_band_labels:
            return None, None
        yvals = np.asarray(dataf['yvals'][label], dtype=float)
        yerrs = np.asarray(dataf['yerrs'][label], dtype=float)
        raw_dn = yvals - yerrs[0]
        raw_up = yvals + yerrs[1]
        band_dn = np.minimum(raw_dn, raw_up)
        band_up = np.maximum(raw_dn, raw_up)
        if styles[label]['drawstyle']:
            band_dn = np.insert(band_dn, 0, band_dn[0])
            band_up = np.insert(band_up, 0, band_up[0])
        return band_dn, band_up

    def _ratio_arrays(label):
        if reference_label is None or label == reference_label:
            return None
        return _safe_ratio(dataf['yvals'][label], dataf['yvals'][reference_label])

    def _ratio_yerrs(label):
        if reference_label is None or label not in dataf.get('yerrs', {}):
            return None
        yerrs = np.asarray(dataf['yerrs'][label], dtype=float)
        if yerrs.shape[0] != 2:
            return None
        ref = np.asarray(dataf['yvals'][reference_label], dtype=float)
        abs_ref = np.abs(ref)
        err_dn = np.divide(yerrs[0], abs_ref, out=np.zeros_like(yerrs[0], dtype=float), where=abs_ref != 0)
        err_up = np.divide(yerrs[1], abs_ref, out=np.zeros_like(yerrs[1], dtype=float), where=abs_ref != 0)
        err_dn = np.where(np.isfinite(err_dn), err_dn, 0.0)
        err_up = np.where(np.isfinite(err_up), err_up, 0.0)
        return np.vstack([err_dn, err_up])

    def _ratio_band_arrays(label):
        if label not in scale_band_labels or reference_label is None:
            return None, None
        yvals = np.asarray(dataf['yvals'][label], dtype=float)
        yerrs = np.asarray(dataf['yerrs'][label], dtype=float)
        ref = np.asarray(dataf['yvals'][reference_label], dtype=float)
        raw_dn = _safe_ratio(yvals - yerrs[0], ref)
        raw_up = _safe_ratio(yvals + yerrs[1], ref)
        band_dn = np.minimum(raw_dn, raw_up)
        band_up = np.maximum(raw_dn, raw_up)
        if styles[label]['ratio0_drawstyle']:
            band_dn = np.insert(band_dn, 0, band_dn[0])
            band_up = np.insert(band_up, 0, band_up[0])
        return band_dn, band_up

    # End Codex scale-envelope patch
    """
).strip("\n")


_MAIN_PANEL_BLOCK = textwrap.dedent(
    """
    # curve from input yoda files in main panel
    _ensure_scale_patch_styles()
    for label, yvals in dataf['yvals'].items():
        if all(np.isnan(v) for v in dataf['yvals'][label]):
            continue
        line_handle = None
        if styles[label]['histstyle']: # draw as histogram
            xpos = dataf['xedges'][label] if styles[label]['drawstyle'] else dataf['xpoints'][label]
            ypos = np.insert(yvals, 0, yvals[0]) if styles[label]['drawstyle'] else yvals
            if styles[label]['fillcolor']: # fill area below curve
                ax.fill_between(xpos, ypos, step='pre',
                                color=styles[label]['fillcolor'],
                                alpha=styles[label]['fillopacity'])
            if label in scale_band_labels:
                band_dn, band_up = _main_band_arrays(label)
                if band_dn is not None:
                    ax.fill_between(
                        xpos,
                        band_dn,
                        band_up,
                        color=_scale_band_color(label),
                        alpha=scale_band_opacity,
                        step='pre' if styles[label]['drawstyle'] else None,
                        zorder=max(styles[label]['zorder'] - 0.25, 0),
                        edgecolor='none',
                    )
            line_handle, = ax.plot(xpos, ypos,
                                   color=styles[label]['color'],
                                   linestyle=styles[label]['linestyle'],
                                   alpha=styles[label]['lineopacity'],
                                   linewidth=styles[label]['linewidth'],
                                   drawstyle=styles[label]['drawstyle'], solid_joinstyle='miter',
                                   zorder=styles[label]['zorder'], label=label)
        handle = line_handle
        if label not in scale_band_labels:
            tmp = ax.errorbar(dataf['xpoints'][label], yvals,
                              xerr=np.array(dataf['xerrs'][label])*styles[label]['xerrorbars'],
                              yerr=np.array(dataf['yerrs'][label])*styles[label]['yerrorbars'],
                              fmt=styles[label]['marker'], capsize=styles[label]['capsize'],
                              alpha=styles[label]['lineopacity'],
                              markersize=styles[label]['markersize'],
                              ecolor=styles[label]['color'],
                              color=styles[label]['color'], zorder=styles[label]['zorder'])
            tmp[-1][0].set_linestyle(styles[label]['linestyle'])
            tmp[-1][0].set_linewidth(styles[label]['linewidth'])
            handle = tmp
        if label in dataf['add_legend_handle'] and handle is not None:
            legend_handles[label] = handle
        for varLabel in dataf['variation_yvals'].keys():
            if varLabel.startswith(label):
                tmp, = ax.plot(dataf['xedges'][label], dataf['variation_yvals'][varLabel],
                               color=styles[label]['color'],
                               linestyle=styles[label]['linestyle'],
                               linewidth=styles[label]['linewidth'],
                               drawstyle='steps-pre', solid_joinstyle='miter',
                               zorder=styles[label]['zorder'], alpha=0.5)
    """
).strip("\n")


_RATIO_PANEL_BLOCK = textwrap.dedent(
    """
    # plots on ratio panel
    if reference_label is not None:
        for label, yvals in dataf['yvals'].items():
            if label == reference_label:
                continue
            ratio_yvals = _ratio_arrays(label)
            if ratio_yvals is None or np.all(np.isnan(ratio_yvals)):
                continue
            if styles[label]['ratio0_histstyle']: # plot as histogram
                xpos = dataf['xedges'][label] if styles[label]['ratio0_drawstyle'] else dataf['xpoints'][label]
                ypos = np.insert(ratio_yvals, 0, ratio_yvals[0]) if styles[label]['ratio0_drawstyle'] else ratio_yvals
                ratio0_ax.plot(xpos, ypos,
                               color=styles[label]['color'],
                               linewidth=styles[label]['ratio0_linewidth'],
                               linestyle=styles[label]['ratio0_linestyle'],
                               alpha=styles[label]['ratio0_lineopacity'],
                               drawstyle=styles[label]['ratio0_drawstyle'], zorder=styles[label]['zorder'],
                               solid_joinstyle='miter')
            else:
                tmp = ratio0_ax.errorbar(dataf['xpoints'][label], ratio_yvals,
                                         xerr=None,
                                         yerr=None,
                                         fmt=styles[label]['ratio0_marker'], capsize=styles[label]['ratio0_capsize'],
                                         alpha=styles[label]['ratio0_lineopacity'],
                                         markersize=styles[label]['ratio0_markersize'],
                                         ecolor=styles[label]['color'],
                                         color=styles[label]['color'])
                tmp[-1][0].set_linestyle(styles[label]['ratio0_linestyle'])
                tmp[-1][0].set_linewidth(styles[label]['ratio0_linewidth'])
            band_dn, band_up = _ratio_band_arrays(label)
            if band_dn is not None:
                ratio0_ax.fill_between(
                    dataf['xedges'][label],
                    band_dn,
                    band_up,
                    color=_scale_band_color(label),
                    alpha=scale_band_opacity,
                    step='pre' if styles[label]['ratio0_drawstyle'] else None,
                    zorder=max(styles[label]['zorder'] - 0.25, 0),
                    edgecolor='none',
                )
    ratio0_ax.axhline(1.0, color='0.5', linewidth=0.8, linestyle='--', zorder=0)
    ratio0_ax.set_ylim(0.5, 1.5)
    """
).strip("\n")


def scale_envelope_helper_block(
    reference_label: str,
    scale_band_color: str = FIXED_ORDER_SCALE_COLOR,
) -> str:
    return (
        _HELPER_BLOCK.replace("__REFERENCE_LABEL__", reference_label)
        .replace("__SCALE_BAND_COLOR__", scale_band_color)
    )


def _ensure_log_axis_formatter(text: str) -> str:
    marker = "def _codex_log_tick_formatter(value, _pos):"
    if marker in text:
        return text

    old_formatter_block = (
        "if ax_yScale == 'log':\n"
        "    ax.yaxis.set_major_formatter(mpl.ticker.LogFormatterMathtext(base=10.0))\n"
        "    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())"
    )
    formatter_block = (
        "def _codex_log_tick_formatter(value, _pos):\n"
        "    if value <= 0 or not np.isfinite(value):\n"
        "        return ''\n"
        "    exponent = np.log10(value)\n"
        "    rounded = int(np.round(exponent))\n"
        "    if not np.isclose(exponent, rounded):\n"
        "        return ''\n"
        "    if rounded == 0:\n"
        "        return r'$1$'\n"
        "    return rf'$10^{{{rounded}}}$'\n"
        "if ax_yScale == 'log':\n"
        "    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(_codex_log_tick_formatter))\n"
        "    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())"
    )
    if old_formatter_block in text:
        return text.replace(old_formatter_block, formatter_block, 1)

    minor_locator = (
        "ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(\n"
        "                           base=10.0, subs=np.arange(0.1, 1, 0.1), numticks=np.inf))"
    )
    locator_and_formatter_block = (
        minor_locator
        + "\n"
        + formatter_block
    )
    if minor_locator not in text:
        return text
    return text.replace(minor_locator, locator_and_formatter_block, 1)


def _apply_legend_position_overrides(text: str) -> str:
    if not any(f"# Analysis object: {path}" in text for path in LOWER_CENTER_LEGEND_OBJECTS):
        return text
    return text.replace(
        "loc='upper left',\n"
        "                        bbox_to_anchor=(0.0, 0.97)))",
        "loc='lower center',\n"
        "                        bbox_to_anchor=(0.5, 0.03)))",
        1,
    )


def patch_scale_envelope_plot_script_text(
    text: str,
    reference_label: str = "POLDIS NLO",
    scale_band_color: str = FIXED_ORDER_SCALE_COLOR,
) -> str:
    data_import_match = re.search(r"exec\(open\(prefix\+'[^']+__data\.py'\)\.read\(\), dataf\)\n", text)
    if not data_import_match:
        raise ValueError("Could not find generated data import block in Rivet plot script")

    patched = text
    if PATCH_SENTINEL in patched:
        patched = re.sub(
            r"\n\n# Codex scale-envelope patch.*?# End Codex scale-envelope patch\n",
            "\n\n",
            patched,
            count=1,
            flags=re.S,
        )
        if PATCH_SENTINEL in patched:
            patched = re.sub(
                r"\n\n# Codex scale-envelope patch.*?\ndef _ratio_band_arrays\(label\):.*?\n    return band_dn, band_up\n",
                "\n\n",
                patched,
                count=1,
                flags=re.S,
            )
    patched = (
        patched[: data_import_match.end()]
        + "\n\n"
        + scale_envelope_helper_block(reference_label, scale_band_color=scale_band_color)
        + patched[data_import_match.end() :]
    )

    patched = re.sub(r"ratio0_ax\.set_ylim\([^\n]+\)\n", "", patched, count=1)
    patched = re.sub(r"ratio0_ax\.set_ylabel\([^\n]+\)", "ratio0_ax.set_ylabel('MC/POLDIS')", patched, count=1)

    main_start = patched.find("# curve from input yoda files in main panel")
    ratio_start = patched.find("# plots on ratio panel")
    legend_start = patched.find("legend_items = list(legend_handles.values())")
    if min(main_start, ratio_start, legend_start) < 0:
        raise ValueError("Could not locate the expected plotting blocks in Rivet plot script")

    patched = patched[:main_start] + _MAIN_PANEL_BLOCK + "\n\n" + patched[ratio_start:]
    ratio_start = patched.find("# plots on ratio panel")
    legend_start = patched.find("legend_items = list(legend_handles.values())")
    patched = patched[:ratio_start] + _RATIO_PANEL_BLOCK + "\n\n" + patched[legend_start:]
    patched = _ensure_log_axis_formatter(patched)
    patched = _apply_legend_position_overrides(patched)
    return patched


def rewrite_scale_envelope_plot_scripts(
    plot_dir: Path,
    python_executable: str | None = None,
    reference_label: str = "POLDIS NLO",
    scale_band_color: str = FIXED_ORDER_SCALE_COLOR,
) -> Tuple[int, int]:
    plot_dir = Path(plot_dir).resolve()
    python_cmd = python_executable or sys.executable
    patched_count = 0
    rerendered_count = 0

    for script_path in sorted(plot_dir.rglob("*.py")):
        if script_path.name.endswith("__data.py"):
            continue
        original = script_path.read_text()
        updated = patch_scale_envelope_plot_script_text(
            original,
            reference_label=reference_label,
            scale_band_color=scale_band_color,
        )
        if updated != original:
            script_path.write_text(updated)
            patched_count += 1
        proc = subprocess.run(
            [python_cmd, str(script_path)],
            cwd=script_path.parent,
            text=True,
            capture_output=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"Failed to rerender patched Rivet plot script {script_path} with rc={proc.returncode}\n"
                f"stdout:\n{proc.stdout}\n\nstderr:\n{proc.stderr}"
            )
        rerendered_count += 1

    return patched_count, rerendered_count


_RATIO_BAND_BLOCK_RE = re.compile(r"^ratio_band_edges = \{.*?^\}\n?", flags=re.M | re.S)
_RATIO_YERROR_RE = re.compile(r"('ratio\d+_yerrorbars'\s*:\s*)0(\s*,)")
_NOSCALE_PATCH_SENTINEL = "# Codex no-scale ratio patch"

_NOSCALE_HELPER_BLOCK = textwrap.dedent(
    """

    # Codex no-scale ratio patch
    main_reference_label = dataf['add_legend_handle'][0] if dataf.get('add_legend_handle') else next(iter(dataf['yvals']))
    try:
        _initial_ratio_y_lims = tuple(ratio0_ax.get_ylim())
    except Exception:
        _initial_ratio_y_lims = (0.5, 1.5)

    def _is_none_ratio_reference_label(label):
        normalized = str(label).upper()
        return normalized == 'NONE' or 'NOSPIN-UNPOL' in normalized

    ratio_reference_label = next(
        (
            label
            for label in dataf.get('add_legend_handle', [])
            if _is_none_ratio_reference_label(label)
        ),
        None,
    )
    if ratio_reference_label is None:
        ratio_reference_label = next(
            (label for label in dataf['yvals'] if _is_none_ratio_reference_label(label)),
            main_reference_label,
        )

    def _no_scale_ratio_y_lims():
        if (
            isinstance(_initial_ratio_y_lims, tuple)
            and len(_initial_ratio_y_lims) == 2
            and all(np.isfinite(value) for value in _initial_ratio_y_lims)
            and _initial_ratio_y_lims[0] < _initial_ratio_y_lims[1]
        ):
            return _initial_ratio_y_lims
        return (0.5, 1.5)

    def _main_xpoints(label):
        xpts = list(dataf['xpoints'][label])
        target = len(dataf['yvals'][label])
        if len(xpts) == target:
            return xpts
        ref = list(dataf['xpoints'].get(main_reference_label, xpts))
        return ref[:target]

    def _main_xedges(label):
        xeds = list(dataf['xedges'][label])
        target = len(dataf['yvals'][label]) + 1
        if len(xeds) == target:
            return xeds
        ref = list(dataf['xedges'].get(main_reference_label, xeds))
        return ref[:target]

    def _main_xerrs(label):
        xerrs = dataf['xerrs'][label]
        target = len(dataf['yvals'][label])
        if len(xerrs[0]) == target and len(xerrs[1]) == target:
            return xerrs
        ref = dataf['xerrs'].get(main_reference_label, xerrs)
        return [list(ref[0])[:target], list(ref[1])[:target]]

    def _safe_ratio(num, den):
        num = np.asarray(num, dtype=float)
        den = np.asarray(den, dtype=float)
        return np.divide(num, den, out=np.full_like(num, np.nan, dtype=float), where=den != 0)

    def _ratio_yvals(label):
        yvals = np.asarray(dataf['yvals'][label], dtype=float)
        if ratio_reference_label is None:
            return np.full_like(yvals, np.nan, dtype=float)
        if label == ratio_reference_label:
            return np.ones_like(yvals, dtype=float)
        ref = np.asarray(dataf['yvals'][ratio_reference_label], dtype=float)
        return _safe_ratio(yvals, ref)

    def _ratio_yerrs(label):
        yvals = np.asarray(dataf['yvals'][label], dtype=float)
        zeros = np.zeros_like(yvals, dtype=float)
        if ratio_reference_label is None or label == ratio_reference_label:
            return [zeros, zeros]

        ref = np.asarray(dataf['yvals'][ratio_reference_label], dtype=float)
        yerrs = dataf['yerrs'][label]
        ref_yerrs = dataf['yerrs'][ratio_reference_label]
        num_dn = np.asarray(yerrs[0], dtype=float)
        num_up = np.asarray(yerrs[1], dtype=float)
        ref_dn = np.asarray(ref_yerrs[0], dtype=float)
        ref_up = np.asarray(ref_yerrs[1], dtype=float)

        abs_num = np.abs(yvals)
        abs_ref = np.abs(ref)
        rel_num_dn = np.divide(num_dn, abs_num, out=np.zeros_like(num_dn), where=abs_num != 0)
        rel_num_up = np.divide(num_up, abs_num, out=np.zeros_like(num_up), where=abs_num != 0)
        rel_ref_dn = np.divide(ref_dn, abs_ref, out=np.zeros_like(ref_dn), where=abs_ref != 0)
        rel_ref_up = np.divide(ref_up, abs_ref, out=np.zeros_like(ref_up), where=abs_ref != 0)

        ratio_abs = np.abs(_ratio_yvals(label))
        err_dn = ratio_abs * np.sqrt(rel_num_dn**2 + rel_ref_up**2)
        err_up = ratio_abs * np.sqrt(rel_num_up**2 + rel_ref_dn**2)
        err_dn = np.where(np.isfinite(err_dn), err_dn, 0.0)
        err_up = np.where(np.isfinite(err_up), err_up, 0.0)
        return [err_dn, err_up]
    """
).strip("\n")


_NOSCALE_RATIO_PANEL_BLOCK = textwrap.dedent(
    """
    # plots on ratio panel
    for label in dataf['yvals']:
        ratio_yvals = _ratio_yvals(label)
        if all(np.isnan(v) for v in ratio_yvals):
            continue
        if styles[label]['ratio0_histstyle']: # plot as histogram
            xpos = _main_xedges(label) if styles[label]['ratio0_drawstyle'] else _main_xpoints(label)
            ypos = np.insert(ratio_yvals, 0, ratio_yvals[0]) if styles[label]['ratio0_drawstyle'] else ratio_yvals
            ratio0_ax.plot(xpos, ypos,
                           color=styles[label]['color'],
                           linewidth=styles[label]['ratio0_linewidth'],
                           linestyle=styles[label]['ratio0_linestyle'],
                           alpha=styles[label]['ratio0_lineopacity'],
                           drawstyle=styles[label]['ratio0_drawstyle'], zorder=styles[label]['zorder'],
                           solid_joinstyle='miter')
        tmp = ratio0_ax.errorbar(_main_xpoints(label), ratio_yvals,
                                 xerr=None,
                                 yerr=np.array(_ratio_yerrs(label))*styles[label]['ratio0_yerrorbars'],
                                 fmt=styles[label]['ratio0_marker'], capsize=styles[label]['ratio0_capsize'],
                                 alpha=styles[label]['ratio0_lineopacity'],
                                 markersize=styles[label]['ratio0_markersize'],
                                 ecolor=styles[label]['color'],
                                 color=styles[label]['color'])
        tmp[-1][0].set_linestyle(styles[label]['ratio0_linestyle'])
        tmp[-1][0].set_linewidth(styles[label]['ratio0_linewidth'])

    ratio0_ax.axhline(1.0, color='0.5', linewidth=0.8, linestyle='--', zorder=0)
    ratio0_ax.set_ylim(*_no_scale_ratio_y_lims())
    """
).strip("\n")


def patch_no_scale_ratio_data_text(text: str) -> str:
    if "ratio_band_edges" not in text:
        return text
    if _RATIO_BAND_BLOCK_RE.search(text):
        return _RATIO_BAND_BLOCK_RE.sub("ratio_band_edges = {}\n", text, count=1)
    return text


def patch_no_scale_ratio_script_text(text: str) -> str:
    patched = text
    data_import_match = re.search(r"exec\(open\(prefix\+'[^']+__data\.py'\)\.read\(\), dataf\)\n", patched)
    if _NOSCALE_PATCH_SENTINEL in patched:
        patched = re.sub(
            r"\n\n# Codex no-scale ratio patch.*?(?=\nlegend_handles = dict\(\) # keep track of handles for the legend\n)",
            "\n\n",
            patched,
            count=1,
            flags=re.S,
        )
        data_import_match = re.search(r"exec\(open\(prefix\+'[^']+__data\.py'\)\.read\(\), dataf\)\n", patched)
    if data_import_match and _NOSCALE_PATCH_SENTINEL not in patched:
        patched = (
            patched[: data_import_match.end()]
            + "\n\n"
            + _NOSCALE_HELPER_BLOCK
            + patched[data_import_match.end() :]
        )
    patched = _RATIO_YERROR_RE.sub(r"\g<1>1\g<2>", patched)
    patched = re.sub(r"ratio0_ax\.set_ylabel\([^\n]+\)", "ratio0_ax.set_ylabel('MC/None')", patched, count=1)
    patched = patched.replace(
        "dataf['xedges'][label] if styles[label]['drawstyle'] else dataf['xpoints'][label]",
        "_main_xedges(label) if styles[label]['drawstyle'] else _main_xpoints(label)",
    )
    patched = patched.replace(
        "ax.errorbar(dataf['xpoints'][label], yvals,",
        "ax.errorbar(_main_xpoints(label), yvals,",
    )
    patched = patched.replace(
        "xerr=np.array(dataf['xerrs'][label])*styles[label]['xerrorbars'],",
        "xerr=np.array(_main_xerrs(label))*styles[label]['xerrorbars'],",
    )
    patched = patched.replace(
        "ax.plot(dataf['xedges'][label], dataf['variation_yvals'][varLabel],",
        "ax.plot(_main_xedges(label), dataf['variation_yvals'][varLabel],",
    )
    ratio_start = patched.find("# plots on ratio panel")
    legend_start = patched.find("legend_items = list(legend_handles.values())")
    if min(ratio_start, legend_start) < 0:
        raise ValueError("Could not locate the ratio plotting block in Rivet plot script")
    patched = patched[:ratio_start] + _NOSCALE_RATIO_PANEL_BLOCK + "\n\n" + patched[legend_start:]
    patched = _apply_legend_position_overrides(patched)
    return patched


def rewrite_no_scale_ratio_plot_scripts(plot_dir: Path, python_executable: str | None = None) -> Tuple[int, int]:
    plot_dir = Path(plot_dir).resolve()
    python_cmd = python_executable or sys.executable
    patched_count = 0
    rerendered_count = 0

    for data_path in sorted(plot_dir.rglob("*__data.py")):
        original = data_path.read_text()
        updated = patch_no_scale_ratio_data_text(original)
        if updated != original:
            data_path.write_text(updated)
            patched_count += 1

    for script_path in sorted(plot_dir.rglob("*.py")):
        if script_path.name.endswith("__data.py"):
            continue
        original = script_path.read_text()
        updated = patch_no_scale_ratio_script_text(original)
        if updated != original:
            script_path.write_text(updated)
            patched_count += 1

        proc = subprocess.run(
            [python_cmd, str(script_path)],
            cwd=script_path.parent,
            text=True,
            capture_output=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"Failed to rerender patched Rivet plot script {script_path} with rc={proc.returncode}\n"
                f"stdout:\n{proc.stdout}\n\nstderr:\n{proc.stderr}"
            )
        rerendered_count += 1

    return patched_count, rerendered_count
