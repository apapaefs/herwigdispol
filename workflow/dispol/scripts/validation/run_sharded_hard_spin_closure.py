#!/usr/bin/env python3
"""Complete a bounded paired LHE closure, preserving finalized diagnostic inputs.

Recovery does not declare an interrupted native run successful. It reuses only
its finalized, audited hard inputs, records that provenance, and generates the
remaining inputs in fresh small chunks. No production campaign is launched.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from compare_hard_spin_validation import compare


def chunks(remaining, size):
    while remaining:
        count = min(remaining, size)
        yield count
        remaining -= count


def unique_file(directory, pattern):
    paths = list(directory.glob(pattern))
    if len(paths) != 1:
        raise ValueError(f"Expected one {pattern} in {directory}")
    return paths[0]


def audited_native(directory):
    summary_path = unique_file(directory, "*.hard-spin-summary.json")
    summary = json.loads(summary_path.read_text())
    audit = unique_file(directory, "*.hard-spin.jsonl")
    lhe = unique_file(directory, "*.hard.lhe")
    if not summary["polarized_pdf_calls"] or any(
            calls["shower"] for calls in summary["polarized_pdf_calls"].values()):
        raise ValueError(f"Native input has unverified polarized shower calls: {directory}")
    with audit.open() as stream:
        first = json.loads(next(stream))
        if first["hard_process_spin"] or first["hard_links"]:
            raise ValueError(f"Not an audited off-mode native input: {directory}")
        count = 1 + sum(1 for _ in stream)
    with lhe.open() as stream:
        exported = 0
        last = ""
        for line in stream:
            exported += line.strip() == "<event>"
            if line.strip():
                last = line.strip()
    if count != summary["events"] or exported != count or last != "</LesHouchesEvents>":
        raise ValueError(f"Native hard export was not finalized consistently: {directory}")
    event_errors = sum(int(value) for logfile in directory.glob("audit-S*.log")
        for value in re.findall(r"eventerror \((\d+) times\)", logfile.read_text(errors="replace")))
    if event_errors/count > 0.001:
        raise ValueError(f"Native event-error fraction exceeds 0.1%: {directory}")
    summary["recorded_event_errors"] = event_errors
    return audit, lhe, summary


def concatenate(paths, destination):
    with destination.open("xb") as target:
        for path in paths:
            with path.open("rb") as source:
                shutil.copyfileobj(source, target)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--pheno", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--recover", action="append", default=[], metavar="HELICITY=DIRECTORY")
    parser.add_argument("--reuse-directory", type=Path,
                        help="reuse valid components of a prior bounded attempt; never overwrite them")
    parser.add_argument("--events", type=int, default=1_000_000, help="total fixed inputs per helicity")
    parser.add_argument("--chunk-events", type=int, default=100_000)
    parser.add_argument("--jobs", type=int, default=24)
    parser.add_argument("--native-seed-base", type=int, default=9507000)
    parser.add_argument("--lhe-seed-base", type=int, default=9567000)
    args = parser.parse_args()
    if not 0 < args.events <= 1_000_000 or not 0 < args.chunk_events <= 100_000 or not 0 < args.jobs <= 32:
        parser.error("validation is capped at 1M/helicity, 100k/chunk and 32 jobs")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    recovery = {}
    for entry in args.recover:
        label, directory = entry.split("=", 1)
        if label not in ("PP", "PM", "MP", "MM") or label in recovery:
            parser.error("each recovered input needs one distinct physical helicity")
        recovery[label] = Path(directory).resolve()
        invocation = json.loads((recovery[label]/"invocation.json").read_text())
        if invocation.get("helicity") != label or invocation.get("mode") != "off":
            parser.error("recovered native helicity or mode differs from the request")
    runner = Path(__file__).with_name("run_hard_spin_validation.py")
    components = {label: [] for label in ("PP", "PM", "MP", "MM")}
    tasks = []
    for hindex, label in enumerate(components):
        remaining = args.events
        if label in recovery:
            audit, lhe, summary = audited_native(recovery[label])
            remaining -= summary["events"]
            if remaining < 0:
                parser.error("recovered inputs exceed the declared total")
            tasks.append((label, hindex, 0, summary["events"], recovery[label]))
        for index, count in enumerate(chunks(remaining, args.chunk_events), 1):
            tasks.append((label, hindex, index, count, None))
    seeds = []
    for label, hindex, index, count, recovered in tasks:
        seeds.append(args.lhe_seed_base+hindex*100+index)
        if recovered is None:
            seeds.append(args.native_seed_base+hindex*100+index)
    if len(seeds) != len(set(seeds)):
        parser.error("native and LHE seed allocations overlap")

    def command(directory, mode, label, count, seed, lhe=None):
        argv = [sys.executable, str(runner), "--prefix", str(args.prefix),
                "--pheno", str(args.pheno), "--output", str(directory),
                "--events", str(count), "--seed", str(seed), "--helicity", label,
                "--mode", mode, "--compact", "--max-errors", "1000"]
        if lhe is not None:
            argv += ["--lhe", str(lhe), "--allow-input-exhaustion"]
        with (args.output/(directory.name+".driver.log")).open("x") as log:
            subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, check=True)

    def execute(task):
        label, hindex, index, count, recovered = task
        stem = f"{label}-{index:02d}"
        existing_native = args.reuse_directory/(stem+"-native") if args.reuse_directory else None
        native = recovered or (existing_native if existing_native and existing_native.is_dir()
                                else args.output/(stem+"-native"))
        if not native.is_dir():
            command(native, "off", label, count, args.native_seed_base+hindex*100+index)
        audit, lhe, native_summary = audited_native(native)
        if native_summary["events"] != count:
            raise ValueError(f"Incomplete native chunk: {native}")
        def replay_record(directory):
            candidate = unique_file(directory, "*.hard-spin.jsonl")
            summary = json.loads(unique_file(directory, "*.hard-spin-summary.json").read_text())
            expected_sigma = native_summary["sigma_pb"]*summary["events"]/count
            if abs(summary["sigma_pb"]-expected_sigma) > 1.e-8*expected_sigma:
                raise ValueError(f"LHE cross-section/veto normalization mismatch: {directory}")
            return candidate, summary
        replay = args.reuse_directory/(stem+"-lhe") if args.reuse_directory else None
        reused_replay = False
        if replay and replay.is_dir():
            try:
                candidate, replay_summary = replay_record(replay)
                reused_replay = True
            except (ValueError, OSError):
                pass
        if not reused_replay:
            replay = args.output/(stem+"-lhe")
            command(replay, "lhe", label, count, args.lhe_seed_base+hindex*100+index, lhe)
            candidate, replay_summary = replay_record(replay)
        return label, index, audit, candidate, {
            "native_directory": str(native), "lhe_directory": str(replay),
            "recovered_finalized_inputs": recovered is not None,
            "reused_replay": reused_replay,
            "native_events": count, "lhe_events": replay_summary["events"],
            "native_summary": native_summary, "lhe_summary": replay_summary,
            "hard_lhe_sha256": sha256(lhe),
            "native_jsonl_sha256": sha256(audit), "lhe_jsonl_sha256": sha256(candidate),
        }

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(execute, task) for task in tasks]
        for future in as_completed(futures):
            label, index, native, replay, metadata = future.result()
            components[label].append((index, native, replay, metadata))
            print(json.dumps({"completed": f"{label}-{index:02d}",
                              "native_events": metadata["native_events"],
                              "lhe_events": metadata["lhe_events"]}), flush=True)
    results = {}
    for label, entries in components.items():
        entries.sort()
        first, second = args.output/(label+"-native.jsonl"), args.output/(label+"-lhe.jsonl")
        concatenate([entry[1] for entry in entries], first)
        concatenate([entry[2] for entry in entries], second)
        result = compare(first, second, by_identity=True)
        if result["events"] != args.events:
            raise ValueError(f"Fixed-input total differs from declared bound: {label}")
        result["components"] = [entry[3] for entry in entries]
        with (args.output/(label+"-closure.json")).open("x") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
        results[label] = {k: v for k, v in result.items() if k not in ("observables", "components")}
    with (args.output/"closure-summary.json").open("x") as stream:
        json.dump(results, stream, indent=2, allow_nan=False)
    print(json.dumps(results), flush=True)
    return 0 if all(result["passed"] for result in results.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
