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


_HELPER_BLOCK = textwrap.dedent(
    """

    # Codex scale-envelope patch
    scale_band_labels = {label for label in dataf['yvals'] if '.herwig_' in label}
    reference_label = next((label for label in dataf['yvals'] if label.startswith('POLDIS')), None)
    scale_band_opacity = 0.5

    def _safe_ratio(num, den):
        num = np.asarray(num, dtype=float)
        den = np.asarray(den, dtype=float)
        return np.divide(num, den, out=np.full_like(num, np.nan, dtype=float), where=den != 0)

    def _scale_band_color(label):
        return styles.get(label, {}).get('color', '#EE3311')

    def _main_band_arrays(label):
        if label not in scale_band_labels:
            return None, None
        yvals = np.asarray(dataf['yvals'][label], dtype=float)
        yerrs = np.asarray(dataf['yerrs'][label], dtype=float)
        band_dn = np.maximum(yvals - yerrs[0], np.finfo(float).tiny)
        band_up = yvals + yerrs[1]
        if styles[label]['drawstyle']:
            band_dn = np.insert(band_dn, 0, band_dn[0])
            band_up = np.insert(band_up, 0, band_up[0])
        return band_dn, band_up

    def _ratio_arrays(label):
        if reference_label is None or label == reference_label:
            return None
        return _safe_ratio(dataf['yvals'][label], dataf['yvals'][reference_label])

    def _ratio_band_arrays(label):
        if label not in scale_band_labels or reference_label is None:
            return None, None
        yvals = np.asarray(dataf['yvals'][label], dtype=float)
        yerrs = np.asarray(dataf['yerrs'][label], dtype=float)
        ref = np.asarray(dataf['yvals'][reference_label], dtype=float)
        band_dn = _safe_ratio(np.maximum(yvals - yerrs[0], 0.0), ref)
        band_up = _safe_ratio(yvals + yerrs[1], ref)
        if styles[label]['ratio0_drawstyle']:
            band_dn = np.insert(band_dn, 0, band_dn[0])
            band_up = np.insert(band_up, 0, band_up[0])
        return band_dn, band_up
    """
).strip("\n")


_MAIN_PANEL_BLOCK = textwrap.dedent(
    """
    # curve from input yoda files in main panel
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
                                         xerr=np.array(dataf['xerrs'][label])*styles[label]['ratio0_xerrorbars'],
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


def patch_scale_envelope_plot_script_text(text: str) -> str:
    data_import_match = re.search(r"exec\(open\(prefix\+'[^']+__data\.py'\)\.read\(\), dataf\)\n", text)
    if not data_import_match:
        raise ValueError("Could not find generated data import block in Rivet plot script")

    patched = text
    if PATCH_SENTINEL in patched:
        patched = re.sub(
            r"\n\n# Codex scale-envelope patch.*?(?=\nlegend_handles = dict\(\) # keep track of handles for the legend\n)",
            "\n\n",
            patched,
            count=1,
            flags=re.S,
        )
    patched = (
        patched[: data_import_match.end()]
        + "\n\n"
        + _HELPER_BLOCK
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
    return patched


def rewrite_scale_envelope_plot_scripts(plot_dir: Path, python_executable: str | None = None) -> Tuple[int, int]:
    python_cmd = python_executable or sys.executable
    patched_count = 0
    rerendered_count = 0

    for script_path in sorted(plot_dir.rglob("*.py")):
        if script_path.name.endswith("__data.py"):
            continue
        original = script_path.read_text()
        updated = patch_scale_envelope_plot_script_text(original)
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
    ratio_reference_label = next(
        (
            label
            for label in dataf.get('add_legend_handle', [])
            if 'RIVETPS-NOSPIN-UNPOL' in label
        ),
        None,
    )
    if ratio_reference_label is None:
        ratio_reference_label = next(
            (label for label in dataf['yvals'] if 'RIVETPS-NOSPIN-UNPOL' in label),
            main_reference_label,
        )

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
                                 xerr=np.array(_main_xerrs(label))*styles[label]['ratio0_xerrorbars'],
                                 yerr=np.array(_ratio_yerrs(label))*styles[label]['ratio0_yerrorbars'],
                                 fmt=styles[label]['ratio0_marker'], capsize=styles[label]['ratio0_capsize'],
                                 alpha=styles[label]['ratio0_lineopacity'],
                                 markersize=styles[label]['ratio0_markersize'],
                                 ecolor=styles[label]['color'],
                                 color=styles[label]['color'])
        tmp[-1][0].set_linestyle(styles[label]['ratio0_linestyle'])
        tmp[-1][0].set_linewidth(styles[label]['ratio0_linewidth'])

    ratio0_ax.axhline(1.0, color='0.5', linewidth=0.8, linestyle='--', zorder=0)
    ratio0_ax.set_ylim(0.5, 1.5)
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
    return patched


def rewrite_no_scale_ratio_plot_scripts(plot_dir: Path, python_executable: str | None = None) -> Tuple[int, int]:
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
