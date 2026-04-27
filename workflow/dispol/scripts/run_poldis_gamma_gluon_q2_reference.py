#!/usr/bin/env python3.10
"""
Build a focused GAMMA channel-separated POLDIS reference with Q2-binned NLO tables.

This helper is intentionally separate from the default validation runner. It
materializes a campaign-local POLDIS work area, patches only the copied
``user_dijet_rivetplots_p2b_components.f`` template, and runs a polarized:

  * GAMMA-only setup         (IBOSON = 0)
  * selectable parton leg    (ICH = 1 for quark, 2 for gluon)
  * Breit-frame analysis     (IFRAME = 1)

The patched user routine prints machine-readable Q2-binned component tables for
the selected channel's NLO contribution before and after the usual dijet
acceptance. That makes the ``nlo_real_tree`` rows a direct channel-separated
fixed-order reference for either the GAMMA BGF-like or Compton-like
positive-real contribution.
"""

from __future__ import annotations

import argparse
import json
import re
import shlex
import shutil
import subprocess
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence

from extract_dis_out_results import Measurement
from run_poldis_gamma_component_reference import (
    DEFAULT_BASE_DIR,
    DEFAULT_EVENTS,
    DEFAULT_POLDIS_DIR,
    PDF_PROFILES,
    CutWindow,
    WINDOWS,
    parse_component_summaries,
    parse_poldis_totals,
    resolve_pdf_choice,
)


DEFAULT_TAG = "plain60-gamma-channel-q2-reference"
RUNTIME_FILES = ("poldis.f", "gbook.f", "jetalg.f")
Q2_COMPONENT_ORDER = (
    "total",
    "bornproj_seed",
    "nlo_real_tree",
    "nlo_subtraction",
    "nlo_finite_virtual",
    "nlo_finite_collin",
    "nlo_p2b_counter",
)

_FLOAT_RE = r"[+-]?\d+(?:\.\d*)?(?:[DdEe][+-]?\d+)?"
Q2_COMPONENT_RE = re.compile(
    r"^\s*POLDIS_Q2_COMPONENT\s+"
    r"region=(?P<region>[A-Za-z0-9_.-]+)\s+"
    r"component=(?P<component>[A-Za-z0-9_.-]+)\s+"
    r"q2_low=\s*(?P<q2_low>" + _FLOAT_RE + r")\s+"
    r"q2_high=\s*(?P<q2_high>" + _FLOAT_RE + r")\s+"
    r"value=\s*(?P<value>" + _FLOAT_RE + r")\s+"
    r"error=\s*(?P<error>" + _FLOAT_RE + r")\s*$"
)
PROGRESS_RE = re.compile(r"^\s*(?P<events>\d+),\s*ISEED=")
PARTONIC_CHANNELS = ("gluon", "quark")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default=DEFAULT_TAG, help="Campaign tag used for the reference work area.")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Directory containing the DISPOL helpers.",
    )
    parser.add_argument(
        "--poldis-dir",
        type=Path,
        default=DEFAULT_POLDIS_DIR,
        help="Directory containing the POLDIS source tree and compile helper.",
    )
    parser.add_argument(
        "--window",
        choices=tuple(WINDOWS),
        default="interior",
        help="DIS window to evaluate.",
    )
    parser.add_argument(
        "--partonic-channel",
        choices=PARTONIC_CHANNELS,
        default="gluon",
        help="Focused POLDIS channel: ICH=2 for gluon, ICH=1 for quark (default: gluon).",
    )
    parser.add_argument(
        "--pdf-profile",
        choices=tuple(PDF_PROFILES),
        default="nnpdf_paired",
        help="Named PDF profile used for the run.",
    )
    parser.add_argument(
        "--unpolarized-pdf",
        help="Optional override for the unused unpolarized branch in the local template.",
    )
    parser.add_argument(
        "--polarized-diff-pdf",
        help="Optional override for the polarized-difference POLDIS PDF set name.",
    )
    parser.add_argument("--events", type=int, default=DEFAULT_EVENTS, help="Number of POLDIS events to run.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned compile/run commands only.")
    parser.add_argument(
        "--collect-only",
        action="store_true",
        help="Skip compile/run and rebuild the parsed reference from existing logs.",
    )
    return parser


def channel_ich(partonic_channel: str) -> int:
    return {"quark": 1, "gluon": 2}[partonic_channel]


def channel_label(partonic_channel: str) -> str:
    return f"{partonic_channel}-channel"


def channel_real_piece_label(partonic_channel: str) -> str:
    return "BGF-like" if partonic_channel == "gluon" else "Compton-like"


def work_dir(base_dir: Path, tag: str, window: CutWindow, partonic_channel: str) -> Path:
    return base_dir / "campaigns" / tag / f"poldis-gamma-{partonic_channel}-q2-{window.label}"


def run_dir(base_dir: Path, tag: str, window: CutWindow, partonic_channel: str) -> Path:
    return work_dir(base_dir, tag, window, partonic_channel) / "polarized"


def reference_txt_path(base_dir: Path, tag: str, window: CutWindow, partonic_channel: str) -> Path:
    return work_dir(base_dir, tag, window, partonic_channel) / "reference.txt"


def reference_json_path(base_dir: Path, tag: str, window: CutWindow, partonic_channel: str) -> Path:
    return work_dir(base_dir, tag, window, partonic_channel) / "reference.json"


def expected_top_name() -> str:
    return "dijets_pol_GAM_p2bc.top"


def measurement_payload(measurement: Measurement) -> dict:
    return {"value_pb": measurement.value_pb, "error_pb": measurement.error_pb}


def fmt_measurement(payload: Mapping[str, float]) -> str:
    return f"{payload['value_pb']:.10f} +- {payload['error_pb']:.10f} pb"


def parse_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def replace_once(text: str, old: str, new: str, *, label: str) -> str:
    if old not in text:
        raise RuntimeError(f"Failed to patch user source: missing anchor for {label}")
    return text.replace(old, new, 1)


def patch_user_source(
    source_text: str,
    *,
    window: CutWindow,
    partonic_channel: str,
    events: int,
    unpolarized_pdf: str,
    polarized_diff_pdf: str,
) -> str:
    replacements = {
        "IPOL=": "      IPOL=1",
        "IBOSON=": "      IBOSON=0",
        "LEPCH=": "      LEPCH=-1D0",
        "IFRAME=": "      IFRAME=1",
        "ICH=": f"      ICH={channel_ich(partonic_channel)}",
        "INNLO=": "      INNLO=1",
        "NEV=": f"      NEV={events}",
        "YMIN=": f"      YMIN={window.y_min}",
        "YMAX=": f"      YMAX={window.y_max}",
        "Q2MIN=": f"      Q2MIN={window.q2_min:.0f}",
        "Q2MAX=": f"      Q2MAX={window.q2_max:.0f}",
    }
    seen = {key: False for key in replacements}
    seen_unpolarized_pdf = False
    seen_polarized_pdf = False
    patched = []
    pdf_branch: str | None = None

    for line in source_text.splitlines():
        stripped = line.strip()
        if stripped.startswith("IF (IPOL.EQ.0) THEN"):
            pdf_branch = "unpolarized"
            patched.append(line)
            continue
        if stripped.startswith("ELSEIF (IPOL.EQ.1) THEN"):
            pdf_branch = "polarized_diff"
            patched.append(line)
            continue
        if stripped.startswith("ENDIF"):
            pdf_branch = None
            patched.append(line)
            continue
        if stripped.startswith("call InitPDFsetByName("):
            if pdf_branch == "unpolarized":
                patched.append(f'         call InitPDFsetByName("{unpolarized_pdf}")')
                seen_unpolarized_pdf = True
                continue
            if pdf_branch == "polarized_diff":
                patched.append(f'         call InitPDFsetByName("{polarized_diff_pdf}")')
                seen_polarized_pdf = True
                continue
        replaced = False
        for prefix, replacement in replacements.items():
            if stripped.startswith(prefix):
                patched.append(replacement)
                seen[prefix] = True
                replaced = True
                break
        if not replaced:
            patched.append(line)

    missing = [prefix for prefix, ok in seen.items() if not ok]
    if missing:
        raise RuntimeError(f"Failed to patch POLDIS source lines: {missing}")
    if not seen_unpolarized_pdf or not seen_polarized_pdf:
        raise RuntimeError("Failed to patch the PDF-set lines in the copied user template")

    rendered = "\n".join(patched)
    if source_text.endswith("\n"):
        rendered += "\n"

    rendered = replace_once(
        rendered,
        "      DOUBLE PRECISION COMPSUM(3,NCOMP),COMPQR(3,NCOMP)\n"
        "      COMMON /COMPDEMSUM/ COMPSUM,COMPQR\n",
        "      DOUBLE PRECISION COMPSUM(3,NCOMP),COMPQR(3,NCOMP)\n"
        "      INTEGER NQ2COMPBIN\n"
        "      PARAMETER (NQ2COMPBIN=6)\n"
        f"C     Q2-binned NLO component tables for the focused GAMMA {partonic_channel} audit.\n"
        "      DOUBLE PRECISION Q2TOTALPRECUTSUM(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2TOTALPRECUTQR(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2TOTALSELSUM(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2TOTALSELQR(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPPRECUTSUM(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPPRECUTQR(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPSELSUM(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPSELQR(NCOMP,NQ2COMPBIN)\n"
        "      COMMON /COMPDEMSUM/ COMPSUM,COMPQR\n"
        "      COMMON /Q2COMPDEMSUM/ Q2TOTALPRECUTSUM,Q2TOTALPRECUTQR,\n"
        "     $ Q2TOTALSELSUM,Q2TOTALSELQR,Q2COMPPRECUTSUM,Q2COMPPRECUTQR,\n"
        "     $ Q2COMPSELSUM,Q2COMPSELQR\n",
        label="DEFINEPLOTS declarations",
    )

    rendered = replace_once(
        rendered,
        "      DO K=1,NCOMP\n"
        "        DO IORD=1,3\n"
        "          COMPSUM(IORD,K)=0D0\n"
        "          COMPQR(IORD,K)=0D0\n"
        "        ENDDO\n"
        "      ENDDO\n",
        "      DO K=1,NCOMP\n"
        "        DO IORD=1,3\n"
        "          COMPSUM(IORD,K)=0D0\n"
        "          COMPQR(IORD,K)=0D0\n"
        "        ENDDO\n"
        "      ENDDO\n"
        "      DO I=1,NQ2COMPBIN\n"
        "        Q2TOTALPRECUTSUM(I)=0D0\n"
        "        Q2TOTALPRECUTQR(I)=0D0\n"
        "        Q2TOTALSELSUM(I)=0D0\n"
        "        Q2TOTALSELQR(I)=0D0\n"
        "        DO K=1,NCOMP\n"
        "          Q2COMPPRECUTSUM(K,I)=0D0\n"
        "          Q2COMPPRECUTQR(K,I)=0D0\n"
        "          Q2COMPSELSUM(K,I)=0D0\n"
        "          Q2COMPSELQR(K,I)=0D0\n"
        "        ENDDO\n"
        "      ENDDO\n",
        label="DEFINEPLOTS Q2 initialization",
    )

    rendered = replace_once(
        rendered,
        "      INTEGER N,NA,NT,I,J,K,IORD,NCOMP,ICOMP\n"
        "      INTEGER nplots,IPOL,VEC(4),IFRAME\n"
        "      INTEGER COMPONENT_SLOT\n",
        "      INTEGER N,NA,NT,I,J,K,IORD,NCOMP,ICOMP,Q2IBIN\n"
        "      INTEGER nplots,IPOL,VEC(4),IFRAME\n"
        "      INTEGER COMPONENT_SLOT,Q2_COMPONENT_BIN\n",
        label="USER integer declarations",
    )

    rendered = replace_once(
        rendered,
        "      DOUBLE PRECISION COMPSUM(3,NCOMP),COMPQR(3,NCOMP),COMP(3,NCOMP)\n"
        "      COMMON /COMPDEMSUM/ COMPSUM,COMPQR\n"
        "      DATA CS/100*0/\n"
        "      DATA COMP/33*0D0/\n",
        "      DOUBLE PRECISION COMPSUM(3,NCOMP),COMPQR(3,NCOMP),COMP(3,NCOMP)\n"
        "      INTEGER NQ2COMPBIN\n"
        "      PARAMETER (NQ2COMPBIN=6)\n"
        "      DOUBLE PRECISION Q2TOTALPRECUTSUM(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2TOTALPRECUTQR(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2TOTALSELSUM(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2TOTALSELQR(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPPRECUTSUM(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPPRECUTQR(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPSELSUM(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPSELQR(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2TOTALPRECUT(NQ2COMPBIN),Q2TOTALSEL(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPPRECUT(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPSEL(NCOMP,NQ2COMPBIN)\n"
        "      COMMON /COMPDEMSUM/ COMPSUM,COMPQR\n"
        "      COMMON /Q2COMPDEMSUM/ Q2TOTALPRECUTSUM,Q2TOTALPRECUTQR,\n"
        "     $ Q2TOTALSELSUM,Q2TOTALSELQR,Q2COMPPRECUTSUM,Q2COMPPRECUTQR,\n"
        "     $ Q2COMPSELSUM,Q2COMPSELQR\n"
        "      DATA CS/100*0/\n"
        "      DATA COMP/33*0D0/\n"
        "      DATA Q2TOTALPRECUT/6*0D0/\n"
        "      DATA Q2TOTALSEL/6*0D0/\n"
        "      DATA Q2COMPPRECUT/66*0D0/\n"
        "      DATA Q2COMPSEL/66*0D0/\n",
        label="USER Q2 bookkeeping declarations",
    )

    rendered = replace_once(
        rendered,
        "        DO K=1,NCOMP\n"
        "          DO IORD=1,3\n"
        "            COMPSUM(IORD,K)=COMPSUM(IORD,K)+COMP(IORD,K)\n"
        "            COMPQR(IORD,K)=COMPQR(IORD,K)+COMP(IORD,K)**2\n"
        "            COMP(IORD,K)=0D0\n"
        "          ENDDO\n"
        "        ENDDO\n",
        "        DO K=1,NCOMP\n"
        "          DO IORD=1,3\n"
        "            COMPSUM(IORD,K)=COMPSUM(IORD,K)+COMP(IORD,K)\n"
        "            COMPQR(IORD,K)=COMPQR(IORD,K)+COMP(IORD,K)**2\n"
        "            COMP(IORD,K)=0D0\n"
        "          ENDDO\n"
        "        ENDDO\n"
        "        DO I=1,NQ2COMPBIN\n"
        "          Q2TOTALPRECUTSUM(I)=Q2TOTALPRECUTSUM(I)+Q2TOTALPRECUT(I)\n"
        "          Q2TOTALPRECUTQR(I)=Q2TOTALPRECUTQR(I)+Q2TOTALPRECUT(I)**2\n"
        "          Q2TOTALSELSUM(I)=Q2TOTALSELSUM(I)+Q2TOTALSEL(I)\n"
        "          Q2TOTALSELQR(I)=Q2TOTALSELQR(I)+Q2TOTALSEL(I)**2\n"
        "          Q2TOTALPRECUT(I)=0D0\n"
        "          Q2TOTALSEL(I)=0D0\n"
        "          DO K=1,NCOMP\n"
        "            Q2COMPPRECUTSUM(K,I)=Q2COMPPRECUTSUM(K,I)\n"
        "     $       +Q2COMPPRECUT(K,I)\n"
        "            Q2COMPPRECUTQR(K,I)=Q2COMPPRECUTQR(K,I)\n"
        "     $       +Q2COMPPRECUT(K,I)**2\n"
        "            Q2COMPSELSUM(K,I)=Q2COMPSELSUM(K,I)+Q2COMPSEL(K,I)\n"
        "            Q2COMPSELQR(K,I)=Q2COMPSELQR(K,I)+Q2COMPSEL(K,I)**2\n"
        "            Q2COMPPRECUT(K,I)=0D0\n"
        "            Q2COMPSEL(K,I)=0D0\n"
        "          ENDDO\n"
        "        ENDDO\n",
        label="USER event-end Q2 accumulation",
    )

    rendered = replace_once(
        rendered,
        "      ICOMP=COMPONENT_SLOT(NA,NT,N)\n"
        "      IF (ICOMP.GT.0) THEN\n"
        "        COMP(1,ICOMP)=COMP(1,ICOMP)+W(0)\n"
        "        COMP(2,ICOMP)=COMP(2,ICOMP)+W(1)\n"
        "        COMP(3,ICOMP)=COMP(3,ICOMP)+W(2)\n"
        "      ENDIF\n",
        "      ICOMP=COMPONENT_SLOT(NA,NT,N)\n"
        "      IF (ICOMP.GT.0) THEN\n"
        "        COMP(1,ICOMP)=COMP(1,ICOMP)+W(0)\n"
        "        COMP(2,ICOMP)=COMP(2,ICOMP)+W(1)\n"
        "        COMP(3,ICOMP)=COMP(3,ICOMP)+W(2)\n"
        "      ENDIF\n"
        f"C     Keep a Q2-binned view of the {partonic_channel}-only NLO callback pieces.\n"
        "      Q2IBIN=Q2_COMPONENT_BIN(Q2)\n"
        "      IF (Q2IBIN.GT.0) THEN\n"
        "        Q2TOTALPRECUT(Q2IBIN)=Q2TOTALPRECUT(Q2IBIN)+W(1)\n"
        "        IF (ICOMP.GT.0) Q2COMPPRECUT(ICOMP,Q2IBIN)=\n"
        "     $    Q2COMPPRECUT(ICOMP,Q2IBIN)+W(1)\n"
        "      ENDIF\n",
        label="USER precut Q2 fill",
    )

    rendered = replace_once(
        rendered,
        "      IF (zkt(n1).gt.5d0.and.yeta(n1).gt.-3.5d0\n"
        "     $ .and.yeta(n1).lt.3.5d0.and.zkt(n2).gt.4d0\n"
        "     $ .and.yeta(n2).gt.-3.5d0.and.yeta(n2).lt.3.5d0) then\n",
        "      IF (zkt(n1).gt.5d0.and.yeta(n1).gt.-3.5d0\n"
        "     $ .and.yeta(n1).lt.3.5d0.and.zkt(n2).gt.4d0\n"
        "     $ .and.yeta(n2).gt.-3.5d0.and.yeta(n2).lt.3.5d0) then\n"
        "        IF (Q2IBIN.GT.0) THEN\n"
        "          Q2TOTALSEL(Q2IBIN)=Q2TOTALSEL(Q2IBIN)+W(1)\n"
        "          IF (ICOMP.GT.0) Q2COMPSEL(ICOMP,Q2IBIN)=\n"
        "     $      Q2COMPSEL(ICOMP,Q2IBIN)+W(1)\n"
        "        ENDIF\n",
        label="USER selected Q2 fill",
    )

    rendered = replace_once(
        rendered,
        "      INTEGER I,J,K,NEV,NCOMP\n"
        "      PARAMETER (NCOMP=11)\n",
        "      INTEGER I,J,K,NEV,NCOMP,NQ2COMPBIN\n"
        "      PARAMETER (NCOMP=11)\n"
        "      PARAMETER (NQ2COMPBIN=6)\n",
        label="PRINTOUT integer declarations",
    )

    rendered = replace_once(
        rendered,
        "      DOUBLE PRECISION COMPSUM(3,NCOMP),COMPQR(3,NCOMP)\n"
        "      COMMON /COMPDEMSUM/ COMPSUM,COMPQR\n"
        "      CHARACTER*20 COMPLABEL(NCOMP)\n",
        "      DOUBLE PRECISION COMPSUM(3,NCOMP),COMPQR(3,NCOMP)\n"
        "      DOUBLE PRECISION Q2TOTALPRECUTSUM(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2TOTALPRECUTQR(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2TOTALSELSUM(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2TOTALSELQR(NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPPRECUTSUM(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPPRECUTQR(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPSELSUM(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2COMPSELQR(NCOMP,NQ2COMPBIN)\n"
        "      DOUBLE PRECISION Q2EDGES(NQ2COMPBIN+1)\n"
        "      COMMON /COMPDEMSUM/ COMPSUM,COMPQR\n"
        "      COMMON /Q2COMPDEMSUM/ Q2TOTALPRECUTSUM,Q2TOTALPRECUTQR,\n"
        "     $ Q2TOTALSELSUM,Q2TOTALSELQR,Q2COMPPRECUTSUM,Q2COMPPRECUTQR,\n"
        "     $ Q2COMPSELSUM,Q2COMPSELQR\n"
        "      CHARACTER*20 COMPLABEL(NCOMP),COMPMACH(NCOMP)\n",
        label="PRINTOUT Q2 declarations",
    )

    rendered = replace_once(
        rendered,
        "      DOUBLE PRECISION VARNNLO,VARNLO,VARLO\n"
        "      DOUBLE PRECISION SIGNNLO,SIGNLO,SIGLO\n"
        "      DOUBLE PRECISION ERRNNLO,ERRNLO,ERRLO\n"
        "      DOUBLE PRECISION VARCLO,VARCNLO,VARCNNLO\n"
        "      DOUBLE PRECISION SIGCLO,SIGCNLO,SIGCNNLO\n"
        "      DOUBLE PRECISION ERRCLO,ERRCNLO,ERRCNNLO\n",
        "      DOUBLE PRECISION VARNNLO,VARNLO,VARLO\n"
        "      DOUBLE PRECISION SIGNNLO,SIGNLO,SIGLO\n"
        "      DOUBLE PRECISION ERRNNLO,ERRNLO,ERRLO\n"
        "      DOUBLE PRECISION VARCLO,VARCNLO,VARCNNLO\n"
        "      DOUBLE PRECISION SIGCLO,SIGCNLO,SIGCNNLO\n"
        "      DOUBLE PRECISION ERRCLO,ERRCNLO,ERRCNNLO\n"
        "      DOUBLE PRECISION SIGBIN,VARBIN,ERRBIN\n",
        label="PRINTOUT Q2 locals",
    )

    rendered = replace_once(
        rendered,
        "      DATA COMPLABEL /\n"
        "     $ 'bornproj_seed       ',\n"
        "     $ 'nlo_real_tree       ',\n"
        "     $ 'nlo_subtraction     ',\n"
        "     $ 'nlo_finite_virtual  ',\n"
        "     $ 'nlo_finite_collin. ',\n"
        "     $ 'nlo_p2b_counter     ',\n"
        "     $ 'nnlo_real_tree      ',\n"
        "     $ 'nnlo_subtraction    ',\n"
        "     $ 'nnlo_finite_virtual ',\n"
        "     $ 'nnlo_finite_collin. ',\n"
        "     $ 'nnlo_p2b_counter    ' /\n",
        "      DATA COMPLABEL /\n"
        "     $ 'bornproj_seed       ',\n"
        "     $ 'nlo_real_tree       ',\n"
        "     $ 'nlo_subtraction     ',\n"
        "     $ 'nlo_finite_virtual  ',\n"
        "     $ 'nlo_finite_collin. ',\n"
        "     $ 'nlo_p2b_counter     ',\n"
        "     $ 'nnlo_real_tree      ',\n"
        "     $ 'nnlo_subtraction    ',\n"
        "     $ 'nnlo_finite_virtual ',\n"
        "     $ 'nnlo_finite_collin. ',\n"
        "     $ 'nnlo_p2b_counter    ' /\n"
        "      DATA COMPMACH /\n"
        "     $ 'bornproj_seed',\n"
        "     $ 'nlo_real_tree',\n"
        "     $ 'nlo_subtraction',\n"
        "     $ 'nlo_finite_virtual',\n"
        "     $ 'nlo_finite_collin',\n"
        "     $ 'nlo_p2b_counter',\n"
        "     $ 'nnlo_real_tree',\n"
        "     $ 'nnlo_subtraction',\n"
        "     $ 'nnlo_finite_virtual',\n"
        "     $ 'nnlo_finite_collin',\n"
        "     $ 'nnlo_p2b_counter' /\n"
        "      DATA Q2EDGES /25D0,50D0,100D0,200D0,500D0,1000D0,2500D0/\n",
        label="PRINTOUT machine labels",
    )

    rendered = replace_once(
        rendered,
        "      DO K=1,NCOMP\n"
        "        SIGCLO   = COMPSUM(1,K)\n"
        "        SIGCNLO  = COMPSUM(2,K)\n"
        "        SIGCNNLO = COMPSUM(3,K)\n"
        "        IF (ABS(SIGCLO)+ABS(SIGCNLO)+ABS(SIGCNNLO).LE.1D-30) GOTO 200\n"
        "\n"
        "        VARCLO   = COMPQR(1,K) - SIGCLO*SIGCLO/DN\n"
        "        VARCNLO  = COMPQR(2,K) - SIGCNLO*SIGCNLO/DN\n"
        "        VARCNNLO = COMPQR(3,K) - SIGCNNLO*SIGCNNLO/DN\n"
        "        IF (VARCLO   .LT. 0D0) VARCLO   = 0D0\n"
        "        IF (VARCNLO  .LT. 0D0) VARCNLO  = 0D0\n"
        "        IF (VARCNNLO .LT. 0D0) VARCNNLO = 0D0\n"
        "        ERRCLO   = SQRT(VARCLO)\n"
        "        ERRCNLO  = SQRT(VARCNLO)\n"
        "        ERRCNNLO = SQRT(VARCNNLO)\n"
        "\n"
        "        WRITE(*,*) 'POLDIS_COMPONENT',\n"
        "     $   COMPLABEL(K),\n"
        "     $   'LO=', SIGCLO,\n"
        "     $   'LOerr=', ERRCLO,\n"
        "     $   'NLO=', SIGCNLO,\n"
        "     $   'NLOerr=', ERRCNLO,\n"
        "     $   'dNLO=', SIGCNLO-SIGCLO,\n"
        "     $   'NNLO=', SIGCNNLO,\n"
        "     $   'NNLOerr=', ERRCNNLO,\n"
        "     $   'dNNLO=', SIGCNNLO-SIGCNLO\n"
        "  200   CONTINUE\n"
        "      ENDDO\n"
        "\n"
        "      END\n",
        "      DO K=1,NCOMP\n"
        "        SIGCLO   = COMPSUM(1,K)\n"
        "        SIGCNLO  = COMPSUM(2,K)\n"
        "        SIGCNNLO = COMPSUM(3,K)\n"
        "        IF (ABS(SIGCLO)+ABS(SIGCNLO)+ABS(SIGCNNLO).LE.1D-30) GOTO 200\n"
        "\n"
        "        VARCLO   = COMPQR(1,K) - SIGCLO*SIGCLO/DN\n"
        "        VARCNLO  = COMPQR(2,K) - SIGCNLO*SIGCNLO/DN\n"
        "        VARCNNLO = COMPQR(3,K) - SIGCNNLO*SIGCNNLO/DN\n"
        "        IF (VARCLO   .LT. 0D0) VARCLO   = 0D0\n"
        "        IF (VARCNLO  .LT. 0D0) VARCNLO  = 0D0\n"
        "        IF (VARCNNLO .LT. 0D0) VARCNNLO = 0D0\n"
        "        ERRCLO   = SQRT(VARCLO)\n"
        "        ERRCNLO  = SQRT(VARCNLO)\n"
        "        ERRCNNLO = SQRT(VARCNNLO)\n"
        "\n"
        "        WRITE(*,*) 'POLDIS_COMPONENT',\n"
        "     $   COMPLABEL(K),\n"
        "     $   'LO=', SIGCLO,\n"
        "     $   'LOerr=', ERRCLO,\n"
        "     $   'NLO=', SIGCNLO,\n"
        "     $   'NLOerr=', ERRCNLO,\n"
        "     $   'dNLO=', SIGCNLO-SIGCLO,\n"
        "     $   'NNLO=', SIGCNNLO,\n"
        "     $   'NNLOerr=', ERRCNNLO,\n"
        "     $   'dNNLO=', SIGCNNLO-SIGCNLO\n"
        "  200   CONTINUE\n"
        "      ENDDO\n"
        "\n"
        f"      WRITE(*,'(/A)') 'POLDIS {channel_label(partonic_channel)} Q2 component tables'\n"
        "      WRITE(*,'(A)') '======================================='\n"
        "      DO I=1,NQ2COMPBIN\n"
        "        VARBIN=Q2TOTALPRECUTQR(I)-Q2TOTALPRECUTSUM(I)\n"
        "     $   *Q2TOTALPRECUTSUM(I)/DN\n"
        "        IF (VARBIN.LT.0D0) VARBIN=0D0\n"
        "        ERRBIN=SQRT(VARBIN)\n"
        "        WRITE(*,'(A,1X,A,1X,A,1X,A,F10.3,1X,A,F10.3,1X,A,\n"
        "     $   ES24.16,1X,A,ES24.16)')\n"
        "     $   'POLDIS_Q2_COMPONENT',\n"
        "     $   'region=precut',\n"
        "     $   'component=total',\n"
        "     $   'q2_low=',Q2EDGES(I),\n"
        "     $   'q2_high=',Q2EDGES(I+1),\n"
        "     $   'value=',Q2TOTALPRECUTSUM(I),\n"
        "     $   'error=',ERRBIN\n"
        "      ENDDO\n"
        "      DO I=1,NQ2COMPBIN\n"
        "        VARBIN=Q2TOTALSELQR(I)-Q2TOTALSELSUM(I)*Q2TOTALSELSUM(I)/DN\n"
        "        IF (VARBIN.LT.0D0) VARBIN=0D0\n"
        "        ERRBIN=SQRT(VARBIN)\n"
        "        WRITE(*,'(A,1X,A,1X,A,1X,A,F10.3,1X,A,F10.3,1X,A,\n"
        "     $   ES24.16,1X,A,ES24.16)')\n"
        "     $   'POLDIS_Q2_COMPONENT',\n"
        "     $   'region=selected',\n"
        "     $   'component=total',\n"
        "     $   'q2_low=',Q2EDGES(I),\n"
        "     $   'q2_high=',Q2EDGES(I+1),\n"
        "     $   'value=',Q2TOTALSELSUM(I),\n"
        "     $   'error=',ERRBIN\n"
        "      ENDDO\n"
        "      DO K=1,6\n"
        "        DO I=1,NQ2COMPBIN\n"
        "          VARBIN=Q2COMPPRECUTQR(K,I)-Q2COMPPRECUTSUM(K,I)\n"
        "     $     *Q2COMPPRECUTSUM(K,I)/DN\n"
        "          IF (VARBIN.LT.0D0) VARBIN=0D0\n"
        "          ERRBIN=SQRT(VARBIN)\n"
        "          WRITE(*,'(A,1X,A,1X,A,1X,A,F10.3,1X,A,F10.3,1X,A,\n"
        "     $     ES24.16,1X,A,ES24.16)')\n"
        "     $     'POLDIS_Q2_COMPONENT',\n"
        "     $     'region=precut',\n"
        "     $     'component='//COMPMACH(K),\n"
        "     $     'q2_low=',Q2EDGES(I),\n"
        "     $     'q2_high=',Q2EDGES(I+1),\n"
        "     $     'value=',Q2COMPPRECUTSUM(K,I),\n"
        "     $     'error=',ERRBIN\n"
        "        ENDDO\n"
        "      ENDDO\n"
        "      DO K=1,6\n"
        "        DO I=1,NQ2COMPBIN\n"
        "          VARBIN=Q2COMPSELQR(K,I)-Q2COMPSELSUM(K,I)\n"
        "     $     *Q2COMPSELSUM(K,I)/DN\n"
        "          IF (VARBIN.LT.0D0) VARBIN=0D0\n"
        "          ERRBIN=SQRT(VARBIN)\n"
        "          WRITE(*,'(A,1X,A,1X,A,1X,A,F10.3,1X,A,F10.3,1X,A,\n"
        "     $     ES24.16,1X,A,ES24.16)')\n"
        "     $     'POLDIS_Q2_COMPONENT',\n"
        "     $     'region=selected',\n"
        "     $     'component='//COMPMACH(K),\n"
        "     $     'q2_low=',Q2EDGES(I),\n"
        "     $     'q2_high=',Q2EDGES(I+1),\n"
        "     $     'value=',Q2COMPSELSUM(K,I),\n"
        "     $     'error=',ERRBIN\n"
        "        ENDDO\n"
        "      ENDDO\n"
        "\n"
        "      END\n",
        label="PRINTOUT Q2 table emission",
    )

    rendered = replace_once(
        rendered,
        "      END\nC#######################################################################\n",
        "      END\n"
        "C#######################################################################\n"
        "      INTEGER FUNCTION Q2_COMPONENT_BIN(Q2VAL)\n"
        "      IMPLICIT NONE\n"
        "      INTEGER I,NQ2COMPBIN\n"
        "      PARAMETER (NQ2COMPBIN=6)\n"
        "      DOUBLE PRECISION Q2VAL,Q2EDGES(NQ2COMPBIN+1)\n"
        "      DATA Q2EDGES /25D0,50D0,100D0,200D0,500D0,1000D0,2500D0/\n"
        "\n"
        "      Q2_COMPONENT_BIN=0\n"
        "      DO I=1,NQ2COMPBIN\n"
        "        IF (Q2VAL.GE.Q2EDGES(I) .AND. Q2VAL.LT.Q2EDGES(I+1)) THEN\n"
        "          Q2_COMPONENT_BIN=I\n"
        "          RETURN\n"
        "        ENDIF\n"
        "      ENDDO\n"
        "      IF (Q2VAL.EQ.Q2EDGES(NQ2COMPBIN+1)) Q2_COMPONENT_BIN=NQ2COMPBIN\n"
        "      END\n"
        "C#######################################################################\n",
        label="Q2_COMPONENT_BIN insertion",
    )

    return rendered


def materialize_sources(
    poldis_dir: Path,
    target_dir: Path,
    *,
    window: CutWindow,
    partonic_channel: str,
    events: int,
    unpolarized_pdf: str,
    polarized_diff_pdf: str,
) -> None:
    target_dir.mkdir(parents=True, exist_ok=True)
    for filename in RUNTIME_FILES:
        source_path = poldis_dir / filename
        if not source_path.exists():
            raise FileNotFoundError(f"Missing POLDIS runtime source {source_path}")
        shutil.copy2(source_path, target_dir / filename)

    template_path = poldis_dir / "user_dijet_rivetplots_p2b_components.f"
    if not template_path.exists():
        raise FileNotFoundError(f"Missing POLDIS user template {template_path}")

    rendered = patch_user_source(
        template_path.read_text(),
        window=window,
        partonic_channel=partonic_channel,
        events=events,
        unpolarized_pdf=unpolarized_pdf,
        polarized_diff_pdf=polarized_diff_pdf,
    )
    rendered_path = target_dir / "user_dijet_rivetplots_p2b_components.f"
    if not rendered_path.exists() or rendered_path.read_text() != rendered:
        rendered_path.write_text(rendered)


def compile_command(poldis_dir: Path) -> list[str]:
    compile_script = poldis_dir / "compile_dijet_rivetplots_p2b_components"
    if not compile_script.exists():
        raise FileNotFoundError(f"Missing compile helper {compile_script}")
    command = compile_script.read_text().strip()
    if not command:
        raise RuntimeError(f"Compile helper {compile_script} is empty")
    return shlex.split(command)


def run_logged(cmd: Sequence[str], cwd: Path, log_path: Path, dry_run: bool) -> None:
    if dry_run:
        print(shlex.join(list(cmd)))
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as handle:
        handle.write(f"$ {' '.join(cmd)}\n")
        handle.flush()
        proc = subprocess.run(cmd, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed in {cwd}: {shlex.join(list(cmd))}\nSee {log_path}")


def run_with_progress(run_dir: Path, log_path: Path, dry_run: bool) -> None:
    cmd = ["./poldis.x"]
    if dry_run:
        print(shlex.join(cmd))
        return

    with log_path.open("w") as handle:
        handle.write(f"$ {' '.join(cmd)}\n")
        handle.flush()
        proc = subprocess.Popen(
            cmd,
            cwd=run_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            handle.write(line)
            handle.flush()
            if PROGRESS_RE.match(line.rstrip()):
                print(line.rstrip(), flush=True)
        returncode = proc.wait()
    if returncode != 0:
        raise RuntimeError(f"Command failed in {run_dir}: {shlex.join(cmd)}\nSee {log_path}")


def ensure_compiled(poldis_dir: Path, target_dir: Path, dry_run: bool) -> None:
    cmd = compile_command(poldis_dir)
    binary = target_dir / "poldis.x"
    if not dry_run and binary.exists():
        source_mtime = max(
            (target_dir / filename).stat().st_mtime
            for filename in (*RUNTIME_FILES, "user_dijet_rivetplots_p2b_components.f")
        )
        if binary.stat().st_mtime >= source_mtime:
            return
    run_logged(cmd, cwd=target_dir, log_path=target_dir / "compile.log", dry_run=dry_run)


def parse_q2_component_tables(text: str) -> Dict[str, Dict[str, list[dict[str, float]]]]:
    tables: Dict[str, Dict[str, list[dict[str, float]]]] = {"precut": {}, "selected": {}}
    for line in text.splitlines():
        match = Q2_COMPONENT_RE.match(line.strip())
        if not match:
            continue
        region = match.group("region")
        component = match.group("component")
        row = {
            "q2_low": parse_float(match.group("q2_low")),
            "q2_high": parse_float(match.group("q2_high")),
            "value_pb": parse_float(match.group("value")),
            "error_pb": parse_float(match.group("error")),
        }
        tables.setdefault(region, {}).setdefault(component, []).append(row)

    for region in tables.values():
        for rows in region.values():
            rows.sort(key=lambda row: row["q2_low"])
    return tables


def load_payload(
    base_dir: Path,
    tag: str,
    window: CutWindow,
    partonic_channel: str,
    events: int,
    pdf_profile: str,
    pdfs: Mapping[str, str],
) -> dict:
    target_dir = run_dir(base_dir, tag, window, partonic_channel)
    log_path = target_dir / "run.log"
    if not log_path.exists():
        raise FileNotFoundError(f"Missing POLDIS run log {log_path}")
    log_text = log_path.read_text()
    totals = parse_poldis_totals(log_text, context=str(log_path))
    components = parse_component_summaries(log_text)
    q2_tables = parse_q2_component_tables(log_text)
    top_path = target_dir / expected_top_name()

    return {
        "tag": tag,
        "mode": "poldis-gamma-q2-reference",
        "partonic_channel": partonic_channel,
        "window": asdict(window),
        "events": events,
        "pdf_profile": pdf_profile,
        "pdfs": {
            "unpolarized": pdfs["unpolarized"],
            "polarized_diff": pdfs["polarized_diff"],
        },
        "work_dir": str(work_dir(base_dir, tag, window, partonic_channel)),
        "run": {
            "log": str(log_path),
            "top": str(top_path) if top_path.exists() else None,
        },
        "totals": {order: measurement_payload(totals[order]) for order in ("LO", "NLO", "NNLO")},
        "components": components,
        "q2_component_tables": q2_tables,
    }


def render_region_table(payload: dict, region: str) -> list[str]:
    tables = payload["q2_component_tables"].get(region, {})
    lines = [f"{region.capitalize()} Q2 tables", "-" * (len(region) + 10)]
    for component in Q2_COMPONENT_ORDER:
        rows = tables.get(component)
        if not rows:
            continue
        lines.append(component)
        lines.append("  q2_low   q2_high        value_pb          error_pb")
        for row in rows:
            lines.append(
                f"  {row['q2_low']:7.1f}  {row['q2_high']:7.1f}  "
                f"{row['value_pb']: .10e}  {row['error_pb']: .10e}"
            )
        lines.append("")
    if lines[-1] == "":
        lines.pop()
    return lines


def render_reference_text(payload: dict) -> str:
    window = payload["window"]
    partonic_channel = payload["partonic_channel"]
    lines = [
        f"POLDIS GAMMA {channel_label(partonic_channel)} Q2 reference",
        "=" * len(f"POLDIS GAMMA {channel_label(partonic_channel)} Q2 reference"),
        f"Tag: {payload['tag']}",
        f"Partonic channel: {partonic_channel} (ICH={channel_ich(partonic_channel)})",
        f"Window: {window['label']}",
        f"Q^2 in [{window['q2_min']:.0f}, {window['q2_max']:.0f}] GeV^2, y in [{window['y_min']:.1f}, {window['y_max']:.1f}]",
        f"Events: {payload['events']}",
        f"PDF profile: {payload['pdf_profile']}",
        f"Polarized diff PDF: {payload['pdfs']['polarized_diff']}",
        "",
        "Setup",
        "-----",
        f"This is a campaign-local POLDIS run with IPOL=1, IBOSON=0, ICH={channel_ich(partonic_channel)}, IFRAME=1.",
        f"The nlo_real_tree rows are the direct fixed-order GAMMA {channel_label(partonic_channel)} reference for the {channel_real_piece_label(partonic_channel)} positive-real piece.",
        "",
        "Totals",
        "------",
        f"LO:   {fmt_measurement(payload['totals']['LO'])}",
        f"NLO:  {fmt_measurement(payload['totals']['NLO'])}",
        f"NNLO: {fmt_measurement(payload['totals']['NNLO'])}",
        f"log:  {payload['run']['log']}",
        f"top:  {payload['run']['top'] or 'missing'}",
        "",
    ]
    lines.extend(render_region_table(payload, "precut"))
    lines.append("")
    lines.extend(render_region_table(payload, "selected"))
    return "\n".join(lines).rstrip() + "\n"


def write_results(base_dir: Path, tag: str, window: CutWindow, partonic_channel: str, payload: dict) -> None:
    output_dir = work_dir(base_dir, tag, window, partonic_channel)
    output_dir.mkdir(parents=True, exist_ok=True)
    text = render_reference_text(payload)
    reference_txt_path(base_dir, tag, window, partonic_channel).write_text(text)
    reference_json_path(base_dir, tag, window, partonic_channel).write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(text, end="")
    print(f"Wrote text reference: {reference_txt_path(base_dir, tag, window, partonic_channel)}")
    print(f"Wrote JSON reference: {reference_json_path(base_dir, tag, window, partonic_channel)}")


def run(args: argparse.Namespace) -> None:
    base_dir = args.base_dir.resolve()
    poldis_dir = args.poldis_dir.resolve()
    window = WINDOWS[args.window]
    partonic_channel = args.partonic_channel
    pdfs = resolve_pdf_choice(args.pdf_profile, args.unpolarized_pdf, args.polarized_diff_pdf)
    target_dir = run_dir(base_dir, args.tag, window, partonic_channel)

    if not args.collect_only:
        materialize_sources(
            poldis_dir,
            target_dir,
            window=window,
            partonic_channel=partonic_channel,
            events=args.events,
            unpolarized_pdf=pdfs["unpolarized"],
            polarized_diff_pdf=pdfs["polarized_diff"],
        )
        ensure_compiled(poldis_dir, target_dir, dry_run=args.dry_run)
        run_with_progress(target_dir, target_dir / "run.log", dry_run=args.dry_run)
        if args.dry_run:
            return

    payload = load_payload(
        base_dir,
        args.tag,
        window,
        partonic_channel,
        events=args.events,
        pdf_profile=args.pdf_profile,
        pdfs=pdfs,
    )
    write_results(base_dir, args.tag, window, partonic_channel, payload)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.base_dir.resolve().exists():
        raise SystemExit(f"Base directory does not exist: {args.base_dir}")
    if not args.poldis_dir.resolve().exists():
        raise SystemExit(f"POLDIS directory does not exist: {args.poldis_dir}")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
