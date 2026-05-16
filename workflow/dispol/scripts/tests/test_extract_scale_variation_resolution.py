import importlib.util
import sys
import unittest
from pathlib import Path


def scripts_dir() -> Path:
    for root in Path(__file__).resolve().parents:
        for candidate in (
            root / "DISPOL" / "scripts",
            root / "workflow" / "dispol" / "scripts",
            root / "HwPolNotesNew" / "DISPOL",
        ):
            if (candidate / "extract_dis_out_results.py").exists():
                return candidate
    raise RuntimeError("Could not locate extract_dis_out_results.py")


def load_extractor_module():
    script_dir = scripts_dir()
    sys.path.insert(0, str(script_dir))
    try:
        spec = importlib.util.spec_from_file_location(
            "extract_dis_out_results_for_scale_variation_test",
            script_dir / "extract_dis_out_results.py",
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


class ExtractScaleVariationResolutionTests(unittest.TestCase):

    def test_nominal_shard_resolution_rejects_scale_variation_variants(self):
        module = load_extractor_module()

        self.assertTrue(module.variant_matches_tag("S960000-run-s001", "run"))
        self.assertFalse(module.variant_matches_tag("SCALEUP-S960400-run-s001", "run"))
        self.assertFalse(module.variant_matches_tag("SCALEDOWN-S960200-run-s001", "run"))


if __name__ == "__main__":
    unittest.main()
