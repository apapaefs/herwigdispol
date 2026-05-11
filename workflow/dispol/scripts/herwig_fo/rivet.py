from __future__ import annotations

import os
from pathlib import Path
import subprocess


def default_analysis_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "analyses" / "rivet" / "dis"


def build_mc_dis_breit(analysis_dir: Path | None = None, output_dir: Path | None = None) -> Path:
    source_dir = analysis_dir or default_analysis_dir()
    source = source_dir / "MC_DIS_BREIT.cc"
    if not source.exists():
        raise FileNotFoundError(f"Could not find MC_DIS_BREIT source at {source}")
    target_dir = output_dir or source_dir
    target_dir.mkdir(parents=True, exist_ok=True)
    plugin = target_dir / "RivetMC_DIS_BREIT.so"
    subprocess.run(["rivet-build", str(plugin), str(source)], check=True)
    return plugin


def run_rivet(
    hepmc_path: Path,
    yoda_path: Path,
    jetinput: str = "MEPARTONS",
    analysis_plugin: Path | None = None,
    rivetweights: bool = False,
) -> None:
    env = os.environ.copy()
    if analysis_plugin is not None:
        env["RIVET_ANALYSIS_PATH"] = str(analysis_plugin.parent)
    analysis = f"MC_DIS_BREIT:JETINPUT={jetinput}"
    if rivetweights:
        analysis += ":RIVETWEIGHTS=YES"
    subprocess.run(["rivet", "-a", analysis, "-o", str(yoda_path), str(hepmc_path)], check=True, env=env)
