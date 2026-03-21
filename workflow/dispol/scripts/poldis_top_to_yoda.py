#!/usr/bin/env python3
"""
poldis_top_to_yoda.py

Convert POLDIS Topdrawer (.top) output into a Rivet reference YODA file (.yoda.gz)
that overlays with MC BinnedEstimate1D outputs in rivet-mkhtml.

Key points:
- Produces /REF/<analysis>/<name> objects as *BinnedEstimate1D* (not Scatter2D).
- Reads two topdrawer files:
    * dijets_unpol.top (IPOL=0)  -> unpolarized refs (Q2, Pt, ...)
    * dijets_pol.top   (IPOL=1)  -> polarized refs (DQ2, DPt, ...)
- Default order written into the canonical /REF paths is NLO, configurable to NNLO/LO.
- Bin edges inferred from:
    SET LIMITS X low high
    SET SCALE X LOG   (optional)
    number of points in the chosen dataset
- POLDIS Topdrawer dijet plots are booked "per-bin" in user_dijet_rivetplots.f,
  so the printed y values represent bin integrals. We convert them to differential
  bin heights by dividing by the output bin widths before writing the YODA objects.
- Excludes eta1L/eta2L/F2 by default (not in your Rivet analysis), but keeps Eta (=eta*) etc.
"""

import argparse
import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import yoda


# ----------------------------
# Topdrawer parsing
# ----------------------------

_FLOAT_RE = re.compile(r"^[\s]*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)")

def _is_numeric_triplet(line: str) -> bool:
    toks = line.strip().split()
    if len(toks) < 3:
        return False
    return all(_FLOAT_RE.match(toks[i]) for i in range(3))

def _parse_triplet(line: str) -> Tuple[float, float, float]:
    x, y, dy = line.strip().split()[:3]
    return float(x), float(y), float(dy)

def _extract_title(line: str) -> str:
    m = re.search(r"TITLE\s+TOP\s+'([^']*)'", line)
    if m:
        return m.group(1).strip()
    parts = line.split("TITLE TOP", 1)
    if len(parts) == 2:
        return parts[1].strip().strip("'").strip()
    return ""

@dataclass
class Frame:
    title: str = ""
    xlow: Optional[float] = None
    xhigh: Optional[float] = None
    logx: bool = False
    datasets: List[List[Tuple[float, float, float]]] = field(default_factory=list)

def parse_topdrawer(path: str) -> List[Frame]:
    frames: List[Frame] = []
    cur: Optional[Frame] = None
    cur_ds: List[Tuple[float, float, float]] = []
    reading = False

    def flush_ds():
        nonlocal cur_ds, reading
        if cur is not None and cur_ds:
            cur.datasets.append(cur_ds)
        cur_ds = []
        reading = False

    def flush_frame():
        nonlocal cur
        if cur is None:
            return
        if cur_ds:
            cur.datasets.append(cur_ds.copy())
        frames.append(cur)
        cur = None

    with open(path, "rt", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            s = line.strip()

            if s.startswith("NEW FRAME"):
                if cur is not None:
                    flush_ds()
                    flush_frame()
                cur = Frame()
                continue

            if cur is None:
                continue

            if s.startswith("TITLE TOP"):
                cur.title = _extract_title(line)
                continue

            if s.startswith("SET LIMITS X"):
                toks = s.split()
                if len(toks) >= 5:
                    cur.xlow = float(toks[3])
                    cur.xhigh = float(toks[4])
                continue

            if s.startswith("SET SCALE X LOG"):
                cur.logx = True
                continue

            if s.startswith("SET ORDER X Y DY"):
                flush_ds()
                reading = True
                continue

            if s.startswith("HIST"):
                # end of one dataset; next one starts with numbers directly
                flush_ds()
                reading = True
                continue

            if reading and _is_numeric_triplet(line):
                cur_ds.append(_parse_triplet(line))
                continue

    if cur is not None:
        if cur_ds:
            cur.datasets.append(cur_ds.copy())
        frames.append(cur)

    return [fr for fr in frames if fr.title and fr.datasets and fr.xlow is not None and fr.xhigh is not None]


# ----------------------------
# Title mapping -> Rivet names
# ----------------------------

def normalize_title(title: str) -> str:
    t = title.strip().lower().replace(" ", "")
    # strip common order tags
    for suf in ("nnlo", "nlo", "lo"):
        if t.endswith(suf):
            t = t[: -len(suf)]
    return t

def title_to_rivet_name(title: str) -> Optional[str]:
    """
    Map POLDIS plot titles to your Rivet analysis object names.
    """
    t = normalize_title(title)

    # main Rivet observables
    if t.startswith("q2pre"):
        return "Q2PreCut"
    if t.startswith("xpre"):
        return "XBjPreCut"
    if t.startswith("ypre"):
        return "YPreCut"
    if t.startswith("pt1pre"):
        return "pT1PreCut"
    if t.startswith("pt2pre"):
        return "pT2PreCut"
    if t.startswith("q2"):
        return "Q2"
    if t == "x":
        return "XBj"
    if t.startswith("mm") or t.startswith("mjj"):
        return "Mjj"
    if "<pt>" in t or t.startswith("ptavg") or t.startswith("avgpt"):
        return "Pt"
    if t.startswith("eta*") or t.startswith("etastar"):
        return "Eta"
    if t.startswith("zeta") or t.startswith("log(e2)") or t.startswith("log10(") or t.startswith("log10"):
        return "Zeta"
    if t.startswith("pt1"):
        return "pT1"
    if t.startswith("pt2"):
        return "pT2"
    if t.startswith("r21"):
        return "pT2OverpT1"
    if t.startswith("asym"):
        return "pTAsym"

    # keep eta1L/eta2L for consistency in Topdrawer, but they are NOT in Rivet analysis:
    # ignore them here (same for F2)
    return None


# ----------------------------
# Build BinnedEstimate1D refs
# ----------------------------

def make_edges(xlow: float, xhigh: float, nbins: int, logx: bool) -> List[float]:
    if nbins <= 0:
        raise ValueError("nbins must be > 0")
    if not logx:
        w = (xhigh - xlow) / nbins
        return [xlow + i * w for i in range(nbins + 1)]
    if xlow <= 0 or xhigh <= 0:
        raise ValueError("LogX requires xlow,xhigh > 0")
    r = xhigh / xlow
    return [xlow * (r ** (i / nbins)) for i in range(nbins + 1)]


def dataset_edges(frame: Frame, dataset: List[Tuple[float, float, float]]) -> List[float]:
    """
    Reconstruct cell edges from the printed Topdrawer abscissae.

    POLDIS fills these histograms through GFLIN1, which linearly distributes an
    event weight between neighbouring cells. The printed .top values therefore
    behave like nodal values sampled at the booked x positions, not plain box
    bin contents. We keep the booked end points from SET LIMITS X and use
    midpoints between neighbouring printed x values for the internal edges.
    """
    if not dataset:
        return []
    xs = [pt[0] for pt in dataset]
    if len(xs) == 1:
        return [frame.xlow, frame.xhigh]
    edges = [float(frame.xlow)]
    for left, right in zip(xs[:-1], xs[1:]):
        if frame.logx and left > 0.0 and right > 0.0:
            edges.append(math.sqrt(left * right))
        else:
            edges.append(0.5 * (left + right))
    edges.append(float(frame.xhigh))
    return edges


def _eval_linear(xs: List[float], ys: List[float], x: float) -> float:
    if len(xs) == 1:
        return ys[0]
    if x <= xs[0]:
        x1, x2 = xs[0], xs[1]
        y1, y2 = ys[0], ys[1]
    elif x >= xs[-1]:
        x1, x2 = xs[-2], xs[-1]
        y1, y2 = ys[-2], ys[-1]
    else:
        idx = 0
        while idx + 1 < len(xs) and xs[idx + 1] < x:
            idx += 1
        x1, x2 = xs[idx], xs[idx + 1]
        y1, y2 = ys[idx], ys[idx + 1]
    if x2 == x1:
        return y1
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)


def _integrate_linear_piecewise(xs: List[float], ys: List[float], a: float, b: float) -> float:
    if b <= a:
        return 0.0
    if len(xs) == 1:
        return ys[0] * (b - a)

    breaks = [a]
    for x in xs:
        if a < x < b:
            breaks.append(x)
    breaks.append(b)
    area = 0.0
    for left, right in zip(breaks[:-1], breaks[1:]):
        yl = _eval_linear(xs, ys, left)
        yr = _eval_linear(xs, ys, right)
        area += 0.5 * (yl + yr) * (right - left)
    return area


def nodal_to_bin_averages(
    frame: Frame, dataset: List[Tuple[float, float, float]]
) -> Tuple[List[float], List[float], List[float]]:
    """
    Convert the Topdrawer nodal values written by the GFLIN1-backed POLDIS
    histograms into true cell averages on the booked bin edges.

    The returned y values are the average over each output cell; the returned
    uncertainties are propagated assuming uncorrelated nodal errors.
    """
    if not dataset:
        return [], [], []
    xs = [pt[0] for pt in dataset]
    ys = [pt[1] for pt in dataset]
    dys = [pt[2] for pt in dataset]
    edges = dataset_edges(frame, dataset)
    nbins = len(dataset)
    avg_ys: List[float] = []
    avg_dys: List[float] = []

    for ibin in range(nbins):
        left = edges[ibin]
        right = edges[ibin + 1]
        width = right - left
        if width <= 0.0:
            avg_ys.append(0.0)
            avg_dys.append(0.0)
            continue

        weights: List[float] = []
        for inode in range(nbins):
            basis = [0.0] * nbins
            basis[inode] = 1.0
            weights.append(_integrate_linear_piecewise(xs, basis, left, right) / width)

        value = sum(w * y for w, y in zip(weights, ys))
        error = math.sqrt(sum((w * dy) ** 2 for w, dy in zip(weights, dys)))
        avg_ys.append(value)
        avg_dys.append(error)

    return edges, avg_ys, avg_dys

def new_binned_estimate(edges: List[float]):
    # FORCE binned estimate output. If your python-yoda can read MC as BinnedEstimate1D,
    # it must also be able to construct it.
    if hasattr(yoda, "BinnedEstimate1D"):
        return yoda.BinnedEstimate1D(edges)
    # fallback for older naming
    if hasattr(yoda, "Estimate1D"):
        return yoda.Estimate1D(edges)
    raise RuntimeError("This YODA build does not provide BinnedEstimate1D/Estimate1D constructors.")

def set_bin_val_err(binwrap, y: float, dy: float):
    # value setter compatibility
    if hasattr(binwrap, "setVal"):
        binwrap.setVal(float(y))
    elif hasattr(binwrap, "setValue"):
        binwrap.setValue(float(y))
    else:
        raise RuntimeError("BinWrapper has no setVal/setValue")

    dy = 0.0 if (not math.isfinite(dy) or dy < 0.0) else float(dy)

    # error setter compatibility (your build has setErr)
    if hasattr(binwrap, "setErrs"):
        binwrap.setErrs(dy, dy)
    elif hasattr(binwrap, "setErr"):
        binwrap.setErr(dy)
    else:
        # don’t crash if errors unsupported
        pass

def fill_estimate(est, ys: List[float], dys: List[float]):
    bins = est.bins()
    n = min(len(bins), len(ys), len(dys))
    for i in range(n):
        set_bin_val_err(bins[i], ys[i], dys[i])


def to_density(edges: List[float], ys: List[float], dys: List[float]) -> Tuple[List[float], List[float]]:
    """
    Convert per-bin integrals to differential bin heights.
    """
    densities: List[float] = []
    density_errors: List[float] = []
    for index, value in enumerate(ys):
        width = float(edges[index + 1] - edges[index])
        if width <= 0.0:
            densities.append(0.0)
            density_errors.append(0.0)
            continue
        densities.append(value / width)
        density_errors.append(dys[index] / width if index < len(dys) else 0.0)
    return densities, density_errors

def extract_dataset(frame: Frame, order: str) -> List[Tuple[float, float, float]]:
    """
    GBOOK overlay convention: dataset0=NNLO, dataset1=NLO, dataset2=LO
    """
    idx = {"NNLO": 0, "NLO": 1, "LO": 2}[order.upper()]
    if idx < len(frame.datasets):
        return frame.datasets[idx]
    # fallback: last available
    return frame.datasets[-1] if frame.datasets else []

def frames_to_refobjs(frames: List[Frame], analysis: str, prefix: str, order: str) -> Dict[str, object]:
    out: Dict[str, object] = {}

    for fr in frames:
        name = title_to_rivet_name(fr.title)
        if name is None:
            continue

        ds = extract_dataset(fr, order)
        if not ds:
            continue

        nbins = len(ds)
        legacy_edges = make_edges(fr.xlow, fr.xhigh, nbins, fr.logx)
        edges, ys, dys = nodal_to_bin_averages(fr, ds)
        if not edges:
            edges = legacy_edges
            ys = [p[1] for p in ds]
            dys = [p[2] for p in ds]

        ys, dys = to_density(edges, ys, dys)

        est = new_binned_estimate(edges)
        mc_path  = f"/{analysis}/{prefix}{name}"
        ref_path = f"/REF{mc_path}"

        # path setter compatibility
        if hasattr(est, "setPath"):
            est.setPath(ref_path)
        else:
            # most builds have setPath; if not, we can’t proceed safely
            raise RuntimeError("Estimate object has no setPath method")

        fill_estimate(est, ys, dys)

        # optional legend annotation
        try:
            est.setAnnotation("Legend", f"POLDIS {order}")
        except Exception:
            pass

        out[ref_path] = est

    return out


# ----------------------------
# Writing .yoda.gz (old & new APIs)
# ----------------------------

def _yoda_write_compat(filename: str, objects_list):
    """
    Try both YODA write signatures:
      1) yoda.write(filename, objects)
      2) yoda.write(objects, filename)
    """
    try:
        yoda.write(filename, objects_list)
        return
    except TypeError:
        pass

    try:
        yoda.write(objects_list, filename)
        return
    except TypeError as e:
        # If this fails, show the real issue
        raise RuntimeError(f"Could not write YODA file using either signature: {e}")

def write_yoda_gz(objs: Dict[str, object], outpath: str):
    import gzip, os

    objects_list = list(objs.values())

    if outpath.endswith(".gz"):
        tmp = outpath[:-3]
        _yoda_write_compat(tmp, objects_list)

        with open(tmp, "rb") as fin, gzip.open(outpath, "wb") as fout:
            fout.writelines(fin)
        os.remove(tmp)
    else:
        _yoda_write_compat(outpath, objects_list)

def convert_topdrawer_to_ref_objects(
    unpol: str,
    pol: str,
    analysis: str = "MC_DIS_BREIT",
    order: str = "NLO",
) -> Dict[str, object]:
    unpol_frames = parse_topdrawer(unpol)
    pol_frames = parse_topdrawer(pol)

    refobjs: Dict[str, object] = {}
    refobjs.update(frames_to_refobjs(unpol_frames, analysis, prefix="", order=order))
    refobjs.update(frames_to_refobjs(pol_frames, analysis, prefix="D", order=order))

    if not refobjs:
        raise SystemExit("No reference objects created. Check title mapping and input files.")

    # sanity: ensure we are truly producing binned estimates
    bad = [(k, type(v).__name__) for k, v in refobjs.items() if "Estimate" not in type(v).__name__]
    if bad:
        # If this triggers, you will *not* get reliable overlays.
        raise RuntimeError(f"Non-estimate objects found (should not happen): {bad[:5]}")

    return refobjs


def convert_topdrawer_to_yoda(
    unpol: str,
    pol: str,
    out: str,
    analysis: str = "MC_DIS_BREIT",
    order: str = "NLO",
) -> Dict[str, object]:
    refobjs = convert_topdrawer_to_ref_objects(unpol=unpol, pol=pol, analysis=analysis, order=order)
    write_yoda_gz(refobjs, out)
    return refobjs


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unpol", required=True, help="Topdrawer file from IPOL=0 (e.g. dijets_unpol.top)")
    ap.add_argument("--pol", required=True, help="Topdrawer file from IPOL=1 (e.g. dijets_pol.top)")
    ap.add_argument("--analysis", default="MC_DIS_BREIT", help="Analysis name in paths (default: MC_DIS_BREIT)")
    ap.add_argument("--out", required=True, help="Output reference .yoda.gz")
    ap.add_argument("--order", default="NLO", choices=["NNLO", "NLO", "LO"],
                    help="Which order to store in canonical /REF paths (default: NLO)")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    refobjs = convert_topdrawer_to_yoda(
        unpol=args.unpol,
        pol=args.pol,
        out=args.out,
        analysis=args.analysis,
        order=args.order,
    )

    print(f"Wrote {args.out} with {len(refobjs)} objects.")
    for k in sorted(refobjs.keys()):
        print(" ", k, type(refobjs[k]).__name__)


if __name__ == "__main__":
    main()
