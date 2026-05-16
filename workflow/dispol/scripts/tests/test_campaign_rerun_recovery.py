import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path


def scripts_dir() -> Path:
    for root in Path(__file__).resolve().parents:
        for candidate in (
            root / "DISPOL" / "scripts",
            root / "workflow" / "dispol" / "scripts",
            root / "HwPolNotesNew" / "DISPOL",
        ):
            if (candidate / "run_validation_campaign.py").exists():
                return candidate
    raise RuntimeError("Could not locate run_validation_campaign.py")


def load_campaign_module():
    script_dir = scripts_dir()
    module_names = (
        "yoda",
        "analyze_DIS_polarized",
        "poldis_top_to_yoda",
        "extract_dis_out_results",
        "rivet_scale_plot_postprocess",
    )
    sentinel = object()
    previous = {name: sys.modules.get(name, sentinel) for name in module_names}
    sys.path.insert(0, str(script_dir))
    sys.modules["yoda"] = types.ModuleType("yoda")

    try:
        spec = importlib.util.spec_from_file_location(
            "run_validation_campaign_for_rerun_recovery_test",
            script_dir / "run_validation_campaign.py",
        )
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        return module
    finally:
        try:
            sys.path.remove(str(script_dir))
        except ValueError:
            pass
        for name, value in previous.items():
            if value is sentinel:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = value


class CampaignRerunRecoveryTests(unittest.TestCase):

    def test_recovers_successful_shards_from_existing_yoda_outputs(self):
        module = load_campaign_module()
        job = module.JobSpec(
            setup="ALL",
            order="POSNLO",
            helicity="00",
            stem="DIS-POL-POWHEG_00-POSNLO-ALL-RIVETWEIGHTS",
            run_file="DIS-POL-POWHEG_00-POSNLO-ALL-RIVETWEIGHTS.run",
            in_file="DIS-POL-POWHEG_00-POSNLO-ALL-RIVETWEIGHTS.in",
            events=300,
            analysis_variant="RIVETWEIGHTS",
        )
        shards = [
            module.ShardSpec(job=job, shard_index=1, shard_count=3, tag="recover-s001", seed=101, events=100),
            module.ShardSpec(job=job, shard_index=2, shard_count=3, tag="recover-s002", seed=102, events=100),
            module.ShardSpec(job=job, shard_index=3, shard_count=3, tag="recover-s003", seed=103, events=100),
        ]

        with tempfile.TemporaryDirectory() as tmp:
            base_dir = Path(tmp)
            completed_prefix = f"{job.stem}-S101-recover-s001"
            (base_dir / f"{completed_prefix}.yoda.gz").write_text("YODA")
            (base_dir / f"{completed_prefix}.out").write_text("completed")
            (base_dir / f"{job.stem}-S102-recover-s002.out").write_text("started")

            original_collect_artifacts = module.collect_artifacts

            def fail_per_shard_glob(*_args, **_kwargs):
                raise AssertionError("recovery should index the output directory once")

            module.collect_artifacts = fail_per_shard_glob
            try:
                recovered = module.collect_existing_successful_shards_from_filesystem(base_dir, shards)
            finally:
                module.collect_artifacts = original_collect_artifacts

        self.assertEqual([item.spec.tag for item in recovered], ["recover-s001"])
        self.assertEqual(recovered[0].returncode, 0)
        self.assertEqual(recovered[0].out_files[-1].split("/")[-1], f"{completed_prefix}.out")
        self.assertEqual(recovered[0].yoda_files[-1].split("/")[-1], f"{completed_prefix}.yoda.gz")


if __name__ == "__main__":
    unittest.main()
