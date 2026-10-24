#!/usr/bin/env python3
"""Compile a coordinate regression plugin and run it through a real Herwig.

No events are generated and no installation is modified. --library-dir selects
an isolated HwShower; loader output is recorded and checked. --real-pdfs adds
the saved NNPDF failure-state test (requires both named member-0 sets).
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--library-dir", type=Path)
    parser.add_argument("--cxx", required=True,
                        help="C++ compiler matching the selected Herwig ABI")
    parser.add_argument("--real-pdfs", action="store_true")
    args = parser.parse_args()
    prefix, out = args.prefix.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    fixture = Path(__file__).resolve().parent / "fixtures/polarized_qtilde_pdf_coordinates.cc"
    config = subprocess.run([str(prefix / "bin/herwig-config"), "--cppflags"],
                            text=True, capture_output=True)
    # Some installed Herwig/ThePEG config scripts emit valid flags but exit 1.
    if not config.stdout.strip() or not any(flag.startswith("-I") for flag in shlex.split(config.stdout)):
        raise RuntimeError("herwig-config returned no usable include flags: " + config.stderr)
    plugin = out / "PDFCoordinateRegression.so"
    flags = ["-bundle", "-Wl,-undefined,dynamic_lookup"] if sys.platform == "darwin" else ["-shared", "-fPIC"]
    command = shlex.split(args.cxx) + flags + ["-std=c++17", "-O0", "-g"]
    command += shlex.split(config.stdout) + [str(fixture), "-o", str(plugin)]
    build = subprocess.run(command, text=True, capture_output=True)
    (out / "build.log").write_text(build.stdout + build.stderr)
    record = {"compiler_command": command, "config_returncode": config.returncode,
              "fixture": str(fixture),
              "fixture_sha256": sha(fixture), "build_returncode": build.returncode,
              "herwig": str(prefix / "bin/Herwig"),
              "herwig_sha256": sha(prefix / "bin/Herwig")}
    if build.returncode:
        (out / "result.json").write_text(json.dumps(record, indent=2) + "\n")
        print(build.stderr, file=sys.stderr)
        return 2
    card = out / "coordinate-regression.in"
    card.write_text("library ThePEGLHAPDF.so\nlibrary " + str(plugin) + "\n")
    env = os.environ.copy()
    if args.real_pdfs:
        env["HERWIG_COORDINATE_REAL_PDFS"] = "1"
    else:
        env.pop("HERWIG_COORDINATE_REAL_PDFS", None)
    env["DYLD_PRINT_LIBRARIES" if sys.platform == "darwin" else "LD_DEBUG"] = "1" if sys.platform == "darwin" else "libs"
    command = [str(prefix / "bin/Herwig"), "read"]
    library_dir = args.library_dir.resolve() if args.library_dir else prefix / "lib/Herwig"
    if args.library_dir:
        command += ["-L", str(library_dir)]
        if sys.platform.startswith("linux"):
            # PersistentIStream restores saved prepend paths after CLI -L.
            # Save a private repository with the requested path last, so it
            # wins when restored. No installed repository is modified.
            repository = out / "coordinate-regression.rpo"
            setup = out / "coordinate-repository.in"
            setup.write_text("prependpath " + str(library_dir) + "\nsave "
                             + str(repository) + "\n")
            prepare_command = [str(prefix / "bin/Herwig"), "read", str(setup)]
            prepare = subprocess.run(prepare_command, cwd=out, env=env,
                                     text=True, capture_output=True)
            (out / "repository-build.log").write_text(prepare.stdout + prepare.stderr)
            record.update({"repository": str(repository),
                           "repository_prepare_command": prepare_command,
                           "repository_prepare_returncode": prepare.returncode})
            if prepare.returncode or not repository.is_file() or not repository.stat().st_size:
                record["passed"] = False
                (out / "result.json").write_text(json.dumps(record, indent=2) + "\n")
                print("Could not save the private test repository; see repository-build.log", file=sys.stderr)
                return 2
            command += ["--repo", str(repository)]
    command.append(str(card))
    result = subprocess.run(command, cwd=out, env=env, text=True, capture_output=True)
    (out / "stdout.log").write_text(result.stdout)
    (out / "loader.log").write_text(result.stderr)
    match = re.search(r"PDF_COORDINATE_REGRESSION checks=(\d+) failures=(\d+)", result.stdout)
    libraries = {}
    for line in result.stderr.splitlines():
        if sys.platform == "darwin":
            found = re.search(r"(/[^\n]*?/(?:HwShower[^/ ]*|HwMEDIS[^/ ]*|libThePEG[^/ ]*))\s*$", line)
        else:
            # glibc logs attempted search paths too; only accept actual loads.
            found = re.search(r"calling init:\s*(/.*?/(?:HwShower\S*|HwMEDIS\S*|libThePEG\S*))\s*$", line)
        if found and Path(found.group(1)).is_file():
            path = Path(found.group(1)).resolve()
            libraries[str(path)] = sha(path)
    shower_paths = [Path(p) for p in libraries if Path(p).name.startswith("HwShower")]
    # Require proof of library selection, beyond a successful Herwig read.
    loader_ok = len(shower_paths) == 1 and shower_paths[0].parent == library_dir.resolve()
    record.update({"command": command, "returncode": result.returncode,
                   "checks": int(match.group(1)) if match else None,
                   "failures": int(match.group(2)) if match else None,
                   "real_pdfs": args.real_pdfs, "loaded_libraries": libraries,
                   "requested_shower_directory": str(library_dir),
                   "loader_verified": loader_ok})
    record["passed"] = bool(result.returncode == 0 and match and int(match.group(1)) >= 200
                            and int(match.group(2)) == 0 and loader_ok)
    (out / "result.json").write_text(json.dumps(record, indent=2) + "\n")
    print(result.stdout)
    print(json.dumps({k: record[k] for k in ("checks", "failures", "loader_verified", "passed")}, indent=2))
    return 0 if record["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
