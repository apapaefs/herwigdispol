#!/usr/bin/env python3
"""Bounded validation-only hard input export/replay; never launches production."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

SETTINGS = """set /Herwig/Shower/ShowerHandler:MaxPtIsMuF No
set /Herwig/Shower/ShowerHandler:RestrictPhasespace Yes
set /Herwig/Shower/PartnerFinder:PartnerMethod Random
set /Herwig/Shower/PartnerFinder:ScaleChoice Partner
set /Herwig/Shower/KinematicsReconstructor:ReconstructionOption Colour3
set /Herwig/Shower/KinematicsReconstructor:InitialStateReconOption SofterFraction
"""
QUERIES = [
    "/Herwig/Shower/ShowerHandler:SpinCorrelations",
    "/Herwig/Shower/ShowerHandler:MaxPtIsMuF",
    "/Herwig/Shower/ShowerHandler:RestrictPhasespace",
    "/Herwig/Shower/PartnerFinder:PartnerMethod",
    "/Herwig/Shower/PartnerFinder:ScaleChoice",
    "/Herwig/Shower/KinematicsReconstructor:ReconstructionOption",
    "/Herwig/Shower/KinematicsReconstructor:InitialStateReconOption",
    "/Herwig/Shower/ShowerHandler:PDFA",
    "/Herwig/Shower/ShowerHandler:PDFB",
]

def run(command, directory, environment, logfile):
    with (directory / logfile).open("w") as log:
        result = subprocess.run(command, cwd=directory, env=environment,
                                stdout=log, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        raise RuntimeError(f"{command[0]} failed ({result.returncode}); see {directory/logfile}")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True, type=Path)
    parser.add_argument("--pheno", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--events", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=9001001)
    parser.add_argument("--helicity", choices=("PP", "PM", "MP", "MM"), default="PP")
    parser.add_argument("--mode", choices=("on", "off", "default", "lhe"), default="off")
    parser.add_argument("--process", choices=("pp", "dis"), default="pp")
    parser.add_argument("--lhe", type=Path)
    parser.add_argument("--lhe-weight-option", choices=("UnitWeight", "VarWeight"),
                        default="VarWeight", help="VarWeight prevents automatic hard-event skipping")
    parser.add_argument("--allow-input-exhaustion", action="store_true",
                        help="record finite-file shower veto losses; never permit reopening")
    parser.add_argument("--legacy", action="store_true")
    parser.add_argument("--fixtures", action="store_true")
    parser.add_argument("--compact", action="store_true",
                        help="omit final-state strings for larger paired closure runs")
    parser.add_argument("--spin-correlations", choices=("Yes","No"), default="Yes")
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--extra-cxxflags", default="")
    parser.add_argument("--trace-libraries", action="store_true")
    args = parser.parse_args()
    if args.events <= 0 or (args.fixtures and (args.legacy or args.mode == "lhe")):
        parser.error("fixtures require native new-runtime input and positive event count")
    if args.mode == "lhe" and not args.lhe:
        parser.error("--mode lhe requires --lhe")
    if args.process == "dis" and (args.fixtures or args.mode == "lhe"):
        parser.error("DIS is a default-on regression / off-mode rejection test only")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env["PATH"] = str(args.prefix/"bin") + os.pathsep + env.get("PATH","")
    libdirs = [args.prefix/"lib", args.prefix/"lib"/"Herwig", args.prefix/"lib"/"ThePEG"]
    env["LD_LIBRARY_PATH"] = os.pathsep.join(map(str,libdirs)) + os.pathsep + env.get("LD_LIBRARY_PATH","")
    env["DYLD_LIBRARY_PATH"] = os.pathsep.join(map(str,libdirs)) + os.pathsep + env.get("DYLD_LIBRARY_PATH","")
    source = Path(__file__).with_name("HardSpinValidation.cc")
    compile_command = shlex.split(args.compiler) + ["-std=c++17", "-O2", "-shared", "-fPIC",
        f"-I{args.prefix/'include'}", str(source), f"-L{args.prefix/'lib'}", "-lfastjet",
        "-o", "HardSpinValidation.so"]
    if sys.platform == "darwin":
        compile_command += ["-Wl,-undefined,dynamic_lookup"]
    if not args.legacy:
        compile_command += ["-DHARD_SPIN_NEW"]
    if args.fixtures:
        compile_command += ["-DHARD_SPIN_FIXTURES"]
    if args.compact:
        compile_command += ["-DHARD_SPIN_COMPACT"]
    if args.process == "dis":
        compile_command += ["-DHARD_SPIN_DIS"]
    compile_command += shlex.split(args.extra_cxxflags)
    run(compile_command, args.output, env, "compile.log")
    common = (args.pheno/"cards/phenomenology/MC_POLJETSHAPES/MC_POLJETSHAPES-Common.in").read_text()
    common = common.split("read snippets/Rivet.in")[0]
    if args.process == "dis":
        common = Path(__file__).with_name("polarized-dis-regression.in").read_text()
    p1,p2 = {"PP":(1,1),"PM":(1,-1),"MP":(-1,1),"MM":(-1,-1)}[args.helicity]
    common += f"""
set /Herwig/Partons/PPPolarizedExtractor:FirstLongitudinalPolarization {p1}
set /Herwig/Partons/PPPolarizedExtractor:SecondLongitudinalPolarization {p2}
set /Herwig/Shower/ShowerHandler:SpinCorrelations {args.spin_correlations}
set /Herwig/Generators/EventGenerator:RandomNumberGenerator:Seed {args.seed}
set /Herwig/Generators/EventGenerator:MaxErrors 20
"""
    if args.process == "dis":
        handler = "/Herwig/EventHandlers/EventHandler"
    elif args.mode != "lhe":
        common += """insert /Herwig/MatrixElements/SubProcess:MatrixElements[0] /Herwig/MatrixElements/MEQCD2to2
set /Herwig/MatrixElements/MEQCD2to2:Process All
"""
        handler = "/Herwig/EventHandlers/EventHandler"
    else:
        common += f"""
library LesHouches.so
create ThePEG::LesHouchesEventHandler /Herwig/EventHandlers/AuditLHEHandler
create ThePEG::LesHouchesFileReader /Herwig/EventHandlers/AuditLHEReader
create ThePEG::Cuts /Herwig/Cuts/AuditNoCuts
set /Herwig/EventHandlers/AuditLHEHandler:PartonExtractor /Herwig/Partons/PPExtractor
set /Herwig/Partons/PPExtractor:FirstPDF /Herwig/Partons/HardLOPDF
set /Herwig/Partons/PPExtractor:SecondPDF /Herwig/Partons/HardLOPDF
set /Herwig/EventHandlers/AuditLHEHandler:CascadeHandler /Herwig/Shower/ShowerHandler
set /Herwig/EventHandlers/AuditLHEHandler:HadronizationHandler /Herwig/Hadronization/ClusterHadHandler
set /Herwig/EventHandlers/AuditLHEHandler:DecayHandler /Herwig/Decays/DecayHandler
set /Herwig/EventHandlers/AuditLHEHandler:WeightOption {args.lhe_weight_option}
set /Herwig/EventHandlers/AuditLHEReader:FileName {args.lhe.resolve()}
set /Herwig/EventHandlers/AuditLHEReader:AllowedToReOpen No
set /Herwig/EventHandlers/AuditLHEReader:InitPDFs No
set /Herwig/EventHandlers/AuditLHEReader:ReweightPDF No
set /Herwig/EventHandlers/AuditLHEReader:Cuts /Herwig/Cuts/AuditNoCuts
set /Herwig/EventHandlers/AuditLHEReader:PDFA /Herwig/Partons/HardLOPDF
set /Herwig/EventHandlers/AuditLHEReader:PDFB /Herwig/Partons/HardLOPDF
set /Herwig/EventHandlers/AuditLHEReader:MomentumTreatment RescaleEnergy
insert /Herwig/EventHandlers/AuditLHEHandler:LesHouchesReaders 0 /Herwig/EventHandlers/AuditLHEReader
set /Herwig/Generators/EventGenerator:EventHandler /Herwig/EventHandlers/AuditLHEHandler
"""
        handler = "/Herwig/EventHandlers/AuditLHEHandler"
    if not args.legacy and args.mode != "default":
        common += "set /Herwig/Shower/ShowerHandler:HardProcessSpin " + ("No" if args.mode=="off" else "Yes") + "\n"
    common += SETTINGS
    common += f"""
library {args.output/"HardSpinValidation.so"}
create Herwig::AuditPolarizedPDF /Herwig/Partons/AuditPolarizedFirstPDF
create Herwig::AuditPolarizedPDF /Herwig/Partons/AuditPolarizedSecondPDF
set /Herwig/Partons/AuditPolarizedFirstPDF:PDFName NNPDFpol20_nlo_as_01180
set /Herwig/Partons/AuditPolarizedSecondPDF:PDFName NNPDFpol20_nlo_as_01180
set /Herwig/Partons/AuditPolarizedFirstPDF:RemnantHandler /Herwig/Partons/HadronRemnants
set /Herwig/Partons/AuditPolarizedSecondPDF:RemnantHandler /Herwig/Partons/HadronRemnants
set /Herwig/Partons/PPPolarizedExtractor:FirstLongitudinalDifferencePDF /Herwig/Partons/AuditPolarizedFirstPDF
set /Herwig/Partons/PPPolarizedExtractor:SecondLongitudinalDifferencePDF /Herwig/Partons/AuditPolarizedSecondPDF
create Herwig::HardSpinCapture /Herwig/Analysis/HardSpinCapture
create Herwig::HardSpinAudit /Herwig/Analysis/HardSpinAudit
insert {handler}:PreCascadeHandlers 0 /Herwig/Analysis/HardSpinCapture
insert /Herwig/Generators/EventGenerator:AnalysisHandlers 0 /Herwig/Analysis/HardSpinAudit
"""
    for query in QUERIES + ([] if args.legacy else ["/Herwig/Shower/ShowerHandler:HardProcessSpin"]):
        common += f"get {query}\n"
    if args.process == "dis":
        common += "set /Herwig/Partons/EPPolarizedExtractor:SecondLongitudinalDifferencePDF /Herwig/Partons/AuditPolarizedSecondPDF\n"
    common += "saverun audit /Herwig/Generators/EventGenerator\n"
    (args.output/"audit.in").write_text(common)
    (args.output/"invocation.json").write_text(json.dumps(
        {k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}, indent=2)+"\n")
    run([str(args.prefix/"bin/Herwig"),"read","audit.in"],args.output,env,"read.log")
    if args.trace_libraries:
        env["LD_DEBUG"]="libs"
        env["DYLD_PRINT_LIBRARIES"]="1"
    exhausted = False
    try:
        run([str(args.prefix/"bin/Herwig"),"run","audit.run","-N",str(args.events),"-s",str(args.seed)],
            args.output,env,"run.log")
    except RuntimeError:
        output = (args.output/"run.log").read_text(errors="replace")
        exhausted = (args.mode == "lhe" and args.allow_input_exhaustion and
                     "More events requested than available in LesHouchesReader" in output)
        if not exhausted:
            raise
    summaries = list(args.output.glob("audit*.hard-spin-summary.json"))
    if len(summaries) != 1:
        raise RuntimeError("expected one completed audit summary")
    summary = json.loads(summaries[0].read_text())
    if summary["events"] != args.events and not exhausted:
        raise RuntimeError(f"incomplete audit: {summary}")
    if exhausted:
        print("Finite LHE input exhausted; the identity-aware comparison must account for vetoed hard inputs.")
    print(json.dumps(summary))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
