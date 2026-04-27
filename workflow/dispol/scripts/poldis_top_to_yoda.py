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


def dataset_to_binned_integrals(
    frame: Frame, dataset: List[Tuple[float, float, float]]
) -> Tuple[List[float], List[float], List[float]]:
    """
    user_dijet_rivetplots fills these observables with GFILL1, and GTOP1 writes
    the resulting bin contents at the booked bin centres. Treat the printed
    triplets as ordinary per-bin integrals on the booked grid, not as nodal
    values to be linearly interpolated, otherwise threshold bins get smeared
    into neighbouring bins that should remain empty.
    """
    if not dataset:
        return [], [], []
    if frame.xlow is None or frame.xhigh is None:
        raise RuntimeError(f"Frame {frame.title!r} is missing SET LIMITS X information.")

    edges = make_edges(float(frame.xlow), float(frame.xhigh), len(dataset), frame.logx)
    ys = [pt[1] for pt in dataset]
    dys = [pt[2] for pt in dataset]
    return edges, ys, dys

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


def _bin_value_error(bin_obj) -> Tuple[float, float]:
    if hasattr(bin_obj, "val"):
        value = float(bin_obj.val())
        try:
            error = float(bin_obj.errAvg())
        except Exception:
            try:
                error = 0.5 * (abs(float(bin_obj.errMinus())) + abs(float(bin_obj.errPlus())))
            except Exception:
                error = 0.0
    elif hasattr(bin_obj, "height"):
        value = float(bin_obj.height())
        error = float(bin_obj.heightErr()) if hasattr(bin_obj, "heightErr") else 0.0
    else:
        raise RuntimeError(f"Unsupported YODA bin type {type(bin_obj).__name__}")

    if not math.isfinite(error) or error < 0.0:
        error = 0.0
    return value, error


def _estimate_edges(obj: object) -> List[float]:
    if hasattr(obj, "xEdges"):
        return [float(edge) for edge in obj.xEdges()]
    bins = list(obj.bins()) if hasattr(obj, "bins") else []
    if not bins:
        raise RuntimeError(f"Object of type {type(obj).__name__} does not expose xEdges() or bins().")
    edges = [float(bins[0].xMin())]
    edges.extend(float(bin_obj.xMax()) for bin_obj in bins)
    return edges


def combine_ref_object_sets(
    ref_object_sets: List[Dict[str, object]],
    generated_events: Optional[List[int]] = None,
) -> Dict[str, object]:
    if not ref_object_sets:
        raise ValueError("Need at least one reference-object set to combine.")
    if generated_events is not None and len(generated_events) != len(ref_object_sets):
        raise ValueError("generated_events must match the number of reference-object sets.")

    common_paths = set(ref_object_sets[0])
    for ref_objects in ref_object_sets[1:]:
        common_paths &= set(ref_objects)
    if not common_paths:
        raise RuntimeError("No common reference-object paths found across shard outputs.")

    weights = [float(events) for events in generated_events] if generated_events is not None else [1.0] * len(ref_object_sets)
    if any(weight <= 0.0 for weight in weights):
        raise ValueError("generated_events must be positive for every shard.")
    total_weight = sum(weights)

    combined: Dict[str, object] = {}
    for path in sorted(common_paths):
        objects = [ref_objects[path] for ref_objects in ref_object_sets]
        edges = _estimate_edges(objects[0])
        for other in objects[1:]:
            other_edges = _estimate_edges(other)
            if len(other_edges) != len(edges) or any(
                not math.isclose(left, right, rel_tol=1.0e-12, abs_tol=1.0e-12)
                for left, right in zip(edges, other_edges)
            ):
                raise RuntimeError(f"Mismatched bin grids for reference object {path}")

        combined_values: List[float] = []
        combined_errors: List[float] = []
        per_object_bins = [list(obj.bins()) for obj in objects]
        bin_count = len(edges) - 1
        for bins in per_object_bins:
            if len(bins) != bin_count:
                raise RuntimeError(f"Mismatched bin counts for reference object {path}")

        for bin_index in range(bin_count):
            bin_measurements = [_bin_value_error(bins[bin_index]) for bins in per_object_bins]
            value = sum(weight * measurement[0] for weight, measurement in zip(weights, bin_measurements)) / total_weight
            error = math.sqrt(
                sum((weight * measurement[1]) ** 2 for weight, measurement in zip(weights, bin_measurements))
            ) / total_weight
            combined_values.append(value)
            combined_errors.append(error)

        estimate = new_binned_estimate(edges)
        if hasattr(estimate, "setPath"):
            estimate.setPath(path)
        else:
            raise RuntimeError("Estimate object has no setPath method")
        fill_estimate(estimate, combined_values, combined_errors)
        try:
            estimate.setAnnotation("Legend", "POLDIS NLO")
        except Exception:
            pass
        combined[path] = estimate

    return combined


def _clone_ref_object(template_obj: object, values: List[float], errors: List[float]) -> object:
    edges = _estimate_edges(template_obj)
    estimate = new_binned_estimate(edges)
    path = getattr(template_obj, "path", None)
    if callable(path):
        estimate.setPath(str(path()))
    elif hasattr(template_obj, "path"):
        estimate.setPath(str(template_obj.path))
    else:
        raise RuntimeError("Template reference object has no path information")
    fill_estimate(estimate, values, errors)
    try:
        legend = template_obj.annotation("Legend")
    except Exception:
        legend = None
    if legend is not None:
        try:
            estimate.setAnnotation("Legend", str(legend))
        except Exception:
            pass
    return estimate


def build_ref_object_error_mode(
    nominal_ref_objects: Dict[str, object],
    *,
    error_mode: str = "stat",
    scale_down_ref_objects: Optional[Dict[str, object]] = None,
    scale_up_ref_objects: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    if error_mode not in {"stat", "scale", "stat+scale"}:
        raise ValueError(f"Unsupported reference error mode {error_mode!r}")

    use_scale_errors = error_mode != "stat"
    if use_scale_errors and (scale_down_ref_objects is None or scale_up_ref_objects is None):
        raise ValueError("Scale-based reference error modes require both scale_down_ref_objects and scale_up_ref_objects")

    if use_scale_errors:
        nominal_paths = set(nominal_ref_objects)
        down_paths = set(scale_down_ref_objects or {})
        up_paths = set(scale_up_ref_objects or {})
        if nominal_paths != down_paths or nominal_paths != up_paths:
            raise RuntimeError("Nominal/down/up reference objects do not share identical object paths")

    updated: Dict[str, object] = {}
    for path in sorted(nominal_ref_objects):
        nominal_obj = nominal_ref_objects[path]
        nominal_bins = list(nominal_obj.bins())
        bin_count = len(nominal_bins)
        values: List[float] = []
        errors: List[float] = []

        if use_scale_errors:
            down_obj = (scale_down_ref_objects or {})[path]
            up_obj = (scale_up_ref_objects or {})[path]
            down_edges = _estimate_edges(down_obj)
            up_edges = _estimate_edges(up_obj)
            nominal_edges = _estimate_edges(nominal_obj)
            for other_edges, label in ((down_edges, "down"), (up_edges, "up")):
                if len(other_edges) != len(nominal_edges) or any(
                    not math.isclose(left, right, rel_tol=1.0e-12, abs_tol=1.0e-12)
                    for left, right in zip(nominal_edges, other_edges)
                ):
                    raise RuntimeError(f"Mismatched {label} bin grid for reference object {path}")
            down_bins = list(down_obj.bins())
            up_bins = list(up_obj.bins())
            if len(down_bins) != bin_count or len(up_bins) != bin_count:
                raise RuntimeError(f"Mismatched bin counts for reference object {path}")
        else:
            down_bins = []
            up_bins = []

        for bin_index, nominal_bin in enumerate(nominal_bins):
            nominal_value, nominal_error = _bin_value_error(nominal_bin)
            scale_error = 0.0
            if use_scale_errors:
                down_value, _ = _bin_value_error(down_bins[bin_index])
                up_value, _ = _bin_value_error(up_bins[bin_index])
                scale_error = max(abs(down_value - nominal_value), abs(up_value - nominal_value))

            if error_mode == "stat":
                error = nominal_error
            elif error_mode == "scale":
                error = scale_error
            else:
                error = math.hypot(nominal_error, scale_error)

            values.append(nominal_value)
            errors.append(error)

        updated[path] = _clone_ref_object(nominal_obj, values, errors)

    return updated


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

        edges, ys, dys = dataset_to_binned_integrals(fr, ds)
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
