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
            "run_validation_campaign_for_pdf_profile_test",
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


class PDFProfileInputTests(unittest.TestCase):

    def test_materializer_skips_broken_retired_card_symlinks(self):
        module = load_campaign_module()
        with tempfile.TemporaryDirectory() as tmp:
            base_dir = Path(tmp)
            valid = base_dir / "DIS-POL-POWHEG_00-POSNLO-ALL-RIVETWEIGHTS.in"
            valid.write_text("set /Herwig/Partons/HardNLOPDF:PDFName OLD\n")
            broken = base_dir / "DIS-POL-POWHEG_MM-NEGNLO-ALL-TERMDIAG.in"
            broken.symlink_to(base_dir / "cards" / broken.name)

            root = Path(module.materialize_pdf_profile_inputs(base_dir, "unit", "nnpdf_paired"))

            self.assertTrue((root / valid.name).exists())
            self.assertFalse((root / broken.name).exists())


if __name__ == "__main__":
    unittest.main()
