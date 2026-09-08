#!/usr/bin/env python3
"""Fixed-seed regression or paired-hard-event LHE closure; never a p-value-only gate."""
import argparse
import itertools
import json
import math
from pathlib import Path
from statistics import NormalDist


def records(path):
    with Path(path).open() as stream:
        for line in stream:
            yield json.loads(line)


def same_hard(first, second):
    if len(first) != len(second):
        return False
    return all(len(a) == len(b) and a[:3] == b[:3] and all(
        math.isclose(x, y, rel_tol=1.e-10, abs_tol=1.e-9)
        for x, y in zip(a[3:], b[3:])) for a, b in zip(first, second))


def hard_key(record):
    # Coarse lookup only; same_hard still enforces tight component tolerances.
    return tuple(tuple(leg[:3])+tuple(round(x,5) for x in leg[3:])
                 for leg in record["hard"])


def pairs(first, second, by_identity):
    if not by_identity:
        yield from itertools.zip_longest(records(first),records(second))
        return
    candidates={}
    for item in records(second):
        key=hard_key(item)
        if key in candidates:
            raise ValueError("duplicate LHE hard input: resampling prevents paired closure")
        candidates[key]=item
    seen=set()
    for item in records(first):
        key=hard_key(item)
        if key in seen:
            raise ValueError("duplicate reference hard input")
        seen.add(key)
        candidate=candidates.pop(key,None)
        if candidate is None:
            # A failed finite-file shower contributes zero to an event rate.
            candidate=dict(item, features=[0.]*len(item["features"]),
                           final_state=[], hard_links=0, shower_links=0,
                           vetoed_hard_input=True)
        yield item,candidate
    if candidates:
        raise ValueError("LHE replay contains unknown or changed hard inputs")


def compare(first, second, regression=False, relative_tolerance=0.05, by_identity=False):
    count = mismatches = unequal_final_states = vetoed = 0
    final_states_available = True
    sums = []
    links = {"reference_hard": 0, "candidate_hard": 0,
             "reference_shower": 0, "candidate_shower": 0}
    for a, b in pairs(first, second, by_identity):
        if a is None or b is None:
            raise ValueError("different accepted event counts: cannot pair hard events")
        count += 1
        vetoed += bool(b.get("vetoed_hard_input",False))
        if regression and (not a["final_state"] or not b["final_state"]):
            raise ValueError("compact records cannot certify fixed-seed regression")
        final_states_available &= bool(a["final_state"] and b["final_state"])
        identical = same_hard(a["hard"], b["hard"]) and a["weight"] == b["weight"]
        if not identical:
            mismatches += 1
        unequal_final_states += a["final_state"] != b["final_state"]
        for prefix, record in (("reference", a), ("candidate", b)):
            for kind in ("hard", "shower"):
                links[prefix+"_"+kind] += record[kind+"_links"]
        if not sums:
            sums = [[0.,0.,0.,0.,0.,0.] for _ in a["features"]]
        if len(a["features"]) != 35 or len(b["features"]) != 35:
            raise ValueError("different feature definitions")
        for item, x, y in zip(sums, a["features"], b["features"]):
            if not math.isfinite(x) or not math.isfinite(y):
                raise ValueError("nonfinite validation feature")
            for index, value in enumerate((x,y,x*x,y*y,x*y,(x-y)**2)):
                item[index] += value
    if count < 2:
        raise ValueError("at least two paired hard events are required")
    names = ["baseline"]
    names += [f"n{j}_pt{t}" for t in (2,3,4,5,6,8,10) for j in (3,4)]
    names += [f"jet{j}_tail{t}" for j in (3,4) for t in (2,3,4,5,6,8,10,15,20,30)]
    # Predeclared well-populated rate gates; all rare rates are still reported.
    gates = {"baseline", "n3_pt2", "n3_pt3"}
    z = NormalDist().inv_cdf(1.-0.05/(2*len(gates)))
    results = []
    for name, (sx,sy,sxx,syy,sxy,sd2) in zip(names,sums):
        delta = (sx-sy)/count
        variance = max(0., (sd2-count*delta*delta)/(count*(count-1)))
        covariance = (sxy-sx*sy/count)/(count*(count-1))
        mean = sy/count
        supported = min(sx,sy) >= 100 and 0 < mean < 1
        upper = (abs(delta)+z*math.sqrt(variance))/mean if supported else None
        results.append({"observable":name,"reference":sx/count,"candidate":mean,
            "difference":delta,"paired_standard_error":math.sqrt(variance),
            "shared_hard_covariance_of_means":covariance,
            "simultaneous_95_relative_upper_bound":upper,
            "gate":name in gates,
            "within_tolerance":supported and upper <= relative_tolerance})
    passed = mismatches == 0 and vetoed/count <= 0.001 and (unequal_final_states == 0 if regression else
        all(item["within_tolerance"] for item in results if item["gate"]))
    return {"kind":"fixed-seed regression" if regression else "paired hard-event LHE closure",
        "events":count,"hard_identity_mismatches":mismatches,
        "vetoed_hard_inputs":vetoed,"veto_fraction":vetoed/count,
        "unequal_final_states":unequal_final_states if final_states_available else None,
        "final_state_comparison_available":final_states_available,"spin_links":links,
        "relative_tolerance":relative_tolerance,"simultaneous_confidence":0.95,
        "passed":passed,"observables":results,
        "scope":"rate closure only; rare rates and untested angular shapes are not certified"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference",type=Path)
    parser.add_argument("candidate",type=Path)
    parser.add_argument("--regression",action="store_true")
    parser.add_argument("--by-identity",action="store_true",
                        help="finite LHE replay: reject duplicates; count vetoed inputs as zero-rate events")
    parser.add_argument("--relative-tolerance",type=float,default=0.05)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.regression and args.by_identity:
        parser.error("fixed-seed regression requires exact event ordering")
    result=compare(args.reference,args.candidate,args.regression,args.relative_tolerance,args.by_identity)
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps({key:value for key,value in result.items() if key != "observables"}))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
