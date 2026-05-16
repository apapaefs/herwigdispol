import ast
import importlib.util
import json
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
            "campaign_for_rivetplot_poldis_reference_test",
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


class CampaignRivetplotPoldisReferenceTests(unittest.TestCase):

    def test_reused_reference_campaign_controls_default_poldis_error_mode(self):
        text = (scripts_dir() / "run_validation_campaign.py").read_text()

        self.assertIn("def resolve_plot_poldis_error_mode_requested", text)
        self.assertIn('reference_manifest.get("poldis_error_mode")', text)
        self.assertIn("return resolve_poldis_error_mode_requested(None, campaign_manifest)", text)

    def test_rivetplot_does_not_overlay_poldis_scale_variations_as_curves(self):
        text = (scripts_dir() / "run_validation_campaign.py").read_text()

        self.assertNotIn("POLDIS_scale_up", text)
        self.assertNotIn("POLDIS_scale_down", text)
        self.assertNotIn("reference_{variation_name.lower()}_", text)

    def test_rivetplot_cli_passes_explicit_error_mode_only_when_reusing_refs(self):
        tree = ast.parse((scripts_dir() / "run_validation_campaign.py").read_text())
        functions = {
            node.name: node
            for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef)
        }

        for function_name in ("run_rivetplot_command", "run_full_command"):
            with self.subTest(function=function_name):
                source = ast.unparse(functions[function_name])
                self.assertIn("requested_poldis_error_mode = getattr(args, 'poldis_error_mode', None)", source)
                self.assertIn(
                    "plot_poldis_error_mode = requested_poldis_error_mode if external_ref_campaign else poldis_error_mode",
                    source,
                )
                self.assertIn("poldis_error_mode=plot_poldis_error_mode", source)

    def test_external_standalone_poldis_reference_campaign_does_not_need_root_manifest(self):
        campaign = load_campaign_module()
        with tempfile.TemporaryDirectory() as tmp:
            base_dir = Path(tmp)
            campaign_dir = base_dir / "campaigns" / "plain75"
            campaign_dir.mkdir(parents=True)
            campaign_manifest = {"analysis_variant": "RIVETWEIGHTS"}
            reference_dir = base_dir / "campaigns" / "plain74-poldis-nnpdf_paired"
            poldis_dir = reference_dir / "poldis-cc-broad"
            poldis_dir.mkdir(parents=True)
            (poldis_dir / "reference.json").write_text(
                json.dumps(
                    {
                        "setup": "CC",
                        "reference_yoda": "/refs/reference.stat_scale.yoda.gz",
                        "reference_yoda_error_mode": "stat+scale",
                        "reference_yoda_error_modes": {
                            "stat": "/refs/reference.stat.yoda.gz",
                            "stat+scale": "/refs/reference.stat_scale.yoda.gz",
                        },
                        "reference_yoda_variations": {
                            "ScaleFactorDown": "/refs/reference.ScaleFactorDown.yoda.gz",
                            "ScaleFactorUp": "/refs/reference.ScaleFactorUp.yoda.gz",
                            "nominal": "/refs/reference.nominal.yoda.gz",
                        },
                        "reference_yoda_by_order": {
                            "NLO": "/refs/reference.NLO.stat_scale.yoda.gz",
                        },
                        "reference_yoda_error_modes_by_order": {
                            "NLO": {
                                "stat+scale": "/refs/reference.NLO.stat_scale.yoda.gz",
                            },
                        },
                        "reference_yoda_variations_by_order": {
                            "NLO": {
                                "ScaleFactorUp": "/refs/reference.NLO.ScaleFactorUp.yoda.gz",
                            },
                        },
                    }
                )
            )

            resolved_dir, reference_manifest = campaign.resolve_plot_poldis_reference_manifest(
                base_dir,
                campaign_dir,
                campaign_manifest,
                reference_dir.name,
            )

        self.assertEqual(reference_dir, resolved_dir)
        self.assertEqual("stat+scale", reference_manifest["poldis_error_mode"])
        self.assertEqual(
            {"CC": "/refs/reference.stat_scale.yoda.gz"},
            reference_manifest["dynamic_poldis_reference_yoda_by_setup"],
        )
        self.assertEqual(
            {"CC": {"stat+scale": "/refs/reference.NLO.stat_scale.yoda.gz"}},
            reference_manifest["dynamic_poldis_reference_yoda_error_modes_by_order"]["NLO"],
        )


if __name__ == "__main__":
    unittest.main()
