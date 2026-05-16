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


def repo_root() -> Path:
    for root in Path(__file__).resolve().parents:
        if (root / "HerwigSource" / "Herwig-7.3.0").exists():
            return root
    raise RuntimeError("Could not locate Herwig source root")


def load_module(script_name: str, module_name: str):
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
        spec = importlib.util.spec_from_file_location(module_name, script_dir / script_name)
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


class FakePlotObject:
    def __init__(self):
        self.annotations = {}

    def setAnnotation(self, key, value):
        self.annotations[key] = value


class RivetWeightsShowerModeTests(unittest.TestCase):

    def setUp(self):
        self.campaign = load_module("run_validation_campaign.py", "campaign_for_rivetweights_shower_test")

    def test_new_flag_resolves_to_distinct_variant_and_extractor_flag(self):
        parser = self.campaign.build_parser()
        args = parser.parse_args(["campaign", "-t", "unit", "--rivetweights-shower"])

        variant = self.campaign.analysis_variant_from_args(args)

        self.assertEqual(variant, "RIVETWEIGHTS-SHOWER")
        self.assertEqual(self.campaign.analysis_variant_args(variant), ["--rivetweights-shower"])
        self.assertIn(variant, self.campaign.NC_CORRELATED_WEIGHT_VARIANTS)

    def test_shower_card_rewrite_enables_shower_but_keeps_hadronization_off(self):
        source = "\n".join(
            [
                "set EventGenerator:EventHandler:HadronizationHandler  NULL",
                "set /Herwig/EventHandlers/EventHandler:CascadeHandler NULL",
                "set /Herwig/Shower/ShowerHandler:LimitEmissions HardOnly",
                "set /Herwig/EventHandlers/EventHandler:DecayHandler NULL",
                "set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission No",
                "set /Herwig/MatrixElements/PowhegMEDISNCPol:Contribution PositiveNLO",
                "insert /Herwig/Analysis/Rivet:Analyses 0 MC_DIS_BREIT:JETINPUT=MEPARTONS",
                "saverun DIS-POL-POWHEG_00-POSNLO-ALL-RIVETFO EventGenerator",
            ]
        )

        rendered = self.campaign.rewrite_rivetweights_shower_card_text(
            source,
            "DIS-POL-POWHEG_00-POSNLO-ALL-RIVETWEIGHTS-SHOWER",
        )

        self.assertIn("set /Herwig/MatrixElements/PowhegMEDISNCPol:GenerateRivetWeights Yes", rendered)
        self.assertIn("set EventGenerator:EventHandler:HadronizationHandler NULL", rendered)
        self.assertIn("set /Herwig/EventHandlers/EventHandler:DecayHandler NULL", rendered)
        self.assertNotRegex(rendered, r"(?m)^set\s+(?:/Herwig/Shower/)?ShowerHandler:LimitEmissions\s+HardOnly\b")
        self.assertNotRegex(rendered, r"(?m)^set\s+/Herwig/EventHandlers/EventHandler:CascadeHandler\s+NULL\b")
        self.assertIn("MC_DIS_BREIT:JETINPUT=FULL:RIVETWEIGHTS=YES", rendered)
        self.assertIn("saverun DIS-POL-POWHEG_00-POSNLO-ALL-RIVETWEIGHTS-SHOWER EventGenerator", rendered)

    def test_existing_rivetweights_rewrite_keeps_no_shower_card_shape(self):
        source = "\n".join(
            [
                "set EventGenerator:EventHandler:HadronizationHandler  NULL",
                "set /Herwig/Shower/ShowerHandler:LimitEmissions HardOnly",
                "set /Herwig/EventHandlers/EventHandler:DecayHandler NULL",
                "set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission No",
                "set /Herwig/MatrixElements/PowhegMEDISNCPol:Contribution PositiveNLO",
                "insert /Herwig/Analysis/Rivet:Analyses 0 MC_DIS_BREIT:JETINPUT=MEPARTONS",
                "saverun DIS-POL-POWHEG_00-POSNLO-ALL-RIVETFO EventGenerator",
            ]
        )

        rendered = self.campaign.rewrite_rivetweights_card_text(
            source,
            "DIS-POL-POWHEG_00-POSNLO-ALL-RIVETWEIGHTS",
        )

        self.assertRegex(rendered, r"(?m)^set\s+/Herwig/Shower/ShowerHandler:LimitEmissions\s+HardOnly\b")
        self.assertIn("MC_DIS_BREIT:JETINPUT=MEPARTONS:RIVETWEIGHTS=YES", rendered)
        self.assertNotIn("MC_DIS_BREIT:JETINPUT=FULL:RIVETWEIGHTS=YES", rendered)

    def test_materialize_rivetweights_inputs_supports_charged_current(self):
        source = "\n".join(
            [
                "set /Herwig/MatrixElements/PowhegMEDISCC:UseQ2ScaleInPOWHEGEmission No",
                "set /Herwig/MatrixElements/PowhegMEDISCC:Contribution PositiveNLO",
                "insert /Herwig/Analysis/Rivet:Analyses 0 MC_DIS_BREIT:DISMODE=CC:JETINPUT=MEPARTONS",
                "saverun DIS-POL-POWHEG_00-POSNLO-CC-RIVETFO EventGenerator",
            ]
        )

        with tempfile.TemporaryDirectory() as tmp:
            base_dir = Path(tmp)
            for order in ("POSNLO", "NEGNLO"):
                source_stem = self.campaign.build_logical_run_stem(
                    "CC", order, "00", "RIVETFO", "nominal"
                )
                (base_dir / f"{source_stem}.in").write_text(source)

            self.campaign.materialize_rivetweights_inputs(
                base_dir,
                "unit",
                "hybrid",
                ["CC"],
                dry_run=False,
            )

            target_stem = self.campaign.build_logical_run_stem(
                "CC", "POSNLO", "00", "RIVETWEIGHTS", "nominal"
            )
            rendered = (base_dir / f"{target_stem}.in").read_text()

        self.assertIn("set /Herwig/MatrixElements/PowhegMEDISCC:GenerateRivetWeights Yes", rendered)
        self.assertNotIn("PowhegMEDISNCPol", rendered)
        self.assertIn("MC_DIS_BREIT:DISMODE=CC:JETINPUT=MEPARTONS:RIVETWEIGHTS=YES", rendered)
        self.assertIn("saverun DIS-POL-POWHEG_00-POSNLO-CC-RIVETWEIGHTS EventGenerator", rendered)

    def test_charged_current_scale_variation_targets_cc_matrix_element(self):
        source = "\n".join(
            [
                "set PowhegMEDISNCPol:UseFiniteWidthSpacelikeZPropagator Yes",
                "set /Herwig/MatrixElements/PowhegMEDISCC:UseQ2ScaleInPOWHEGEmission No",
                "set /Herwig/MatrixElements/PowhegMEDISCC:Contribution PositiveNLO",
                "insert SubProcess:MatrixElements[0] PowhegMEDISCC",
                "saverun DIS-POL-POWHEG_00-POSNLO-CC-RIVETWEIGHTS EventGenerator",
            ]
        )
        spec = self.campaign.JobSpec(
            setup="CC",
            order="POSNLO",
            helicity="00",
            stem="DIS-POL-POWHEG_00-POSNLO-CC-RIVETWEIGHTS-SCALEUP",
            run_file="DIS-POL-POWHEG_00-POSNLO-CC-RIVETWEIGHTS-SCALEUP.run",
            in_file="DIS-POL-POWHEG_00-POSNLO-CC-RIVETWEIGHTS-SCALEUP.in",
            events=1,
            analysis_variant="RIVETWEIGHTS",
            scale_variation="ScaleFactorUp",
            scale_factor=2.0,
            source_in_file="DIS-POL-POWHEG_00-POSNLO-CC-RIVETWEIGHTS.in",
        )

        rendered = self.campaign.rewrite_scale_variation_card(source, Path.cwd(), spec)

        self.assertIn("set PowhegMEDISCC:ScaleFactor 2.0", rendered)
        self.assertNotIn("PowhegMEDISNCPol:ScaleFactor", rendered)
        self.assertIn(
            "saverun DIS-POL-POWHEG_00-POSNLO-CC-RIVETWEIGHTS-SCALEUP EventGenerator",
            rendered,
        )

    def test_correlated_rivetweights_rewrite_strips_removed_diagnostic_card_settings(self):
        source = "\n".join(
            [
                "set /Herwig/MatrixElements/PowhegMEDISNCPol:DISDiagnostics Off",
                "set /Herwig/Shower/ShowerHandler:LimitEmissions HardOnly",
                "set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission No",
                "set /Herwig/MatrixElements/PowhegMEDISNCPol:Contribution PositiveNLO",
                "insert /Herwig/Analysis/Rivet:Analyses 0 MC_DIS_BREIT:JETINPUT=MEPARTONS",
                "saverun DIS-POL-POWHEG_00-POSNLO-ALL-RIVETFO EventGenerator",
            ]
        )

        rendered = self.campaign.rewrite_rivetweights_shower_card_text(
            source,
            "DIS-POL-POWHEG_00-POSNLO-ALL-RIVETWEIGHTS-SHOWER",
        )

        self.assertNotIn("DISDiagnostics", rendered)

    def test_rivetweights_shower_plot_style_is_blue(self):
        obj = FakePlotObject()

        self.campaign.apply_variant_style_annotation({"plot": obj}, "RIVETWEIGHTS-SHOWER")

        self.assertEqual(obj.annotations["LineColor"], "blue")
        self.assertEqual(obj.annotations["MarkerColor"], "blue")
        self.assertEqual(obj.annotations["ErrorBandColor"], "blue")
        self.assertEqual(self.campaign.scale_band_color_for_variant("RIVETWEIGHTS-SHOWER"), "blue")

    def test_rivetweights_plot_curve_and_scale_band_style_is_red(self):
        obj = FakePlotObject()

        self.campaign.apply_variant_style_annotation({"plot": obj}, "RIVETWEIGHTS")

        self.assertEqual(obj.annotations["LineColor"], "red")
        self.assertEqual(obj.annotations["MarkerColor"], "red")
        self.assertEqual(obj.annotations["ErrorBandColor"], "red")
        self.assertEqual(self.campaign.scale_band_color_for_variant("RIVETWEIGHTS"), "red")


class RivetWeightsShowerExtractorTests(unittest.TestCase):

    def test_extractor_parses_rivetweights_shower_as_its_own_analysis(self):
        extractor = load_module("extract_dis_out_results.py", "extractor_for_rivetweights_shower_test")

        parsed = extractor.parse_requested_name(
            "DIS-POL-POWHEG_00-POSNLO-ALL-RIVETWEIGHTS-SHOWER.out"
        )

        self.assertEqual(parsed, ("00", "NLO", "ALL", "POSNLO", "RIVETWEIGHTS-SHOWER", ""))

    def test_extractor_uses_00_only_charged_current_defaults_for_rivetweights(self):
        extractor = load_module("extract_dis_out_results.py", "extractor_for_cc_rivetweights_defaults_test")

        runs = extractor.default_base_runs_for_analysis(
            "RIVETWEIGHTS",
            {"CC"},
            extractor.DEFAULT_RUNS,
        )

        self.assertEqual(
            [name for name in runs if "-CC.out" in name],
            [
                "DIS-POL-POWHEG_00-POSNLO-CC.out",
                "DIS-POL-POWHEG_00-NEGNLO-CC.out",
            ],
        )
        self.assertNotIn("DIS-POL-POWHEG_MP-POSNLO-CC.out", runs)
        self.assertNotIn("DIS-POL-LO_00-CC.out", runs)


class ChargedCurrentRivetWeightsSourceTests(unittest.TestCase):

    def test_charged_current_matrix_element_overrides_rivet_weight_born_me2(self):
        root = repo_root()
        header = (
            root
            / "HerwigSource"
            / "Herwig-7.3.0"
            / "MatrixElement"
            / "DIS"
            / "MEChargedCurrentDIS.h"
        ).read_text()
        source = (
            root
            / "HerwigSource"
            / "Herwig-7.3.0"
            / "MatrixElement"
            / "DIS"
            / "MEChargedCurrentDIS.cc"
        ).read_text()

        self.assertIn("rivetWeightBornME2(double Pl, double Pq) const override", header)
        self.assertIn("double MEChargedCurrentDIS::rivetWeightBornME2(double Pl, double Pq) const", source)
        self.assertIn("return me2ForPolarizations(Pl, Pq);", source)

    def test_rivet_weight_values_are_default_weight_multipliers(self):
        root = repo_root()
        source = (
            root
            / "HerwigSource"
            / "Herwig-7.3.0"
            / "MatrixElement"
            / "DIS"
            / "DISBase.cc"
        ).read_text()

        self.assertIn("const double targetWeight = born * raw / sampledDenominator;", source)
        self.assertNotIn("eventWeight * born * raw / sampledDenominator", source)

    def test_powheg_emission_scales_follow_scale_factor(self):
        root = repo_root()
        source = (
            root
            / "HerwigSource"
            / "Herwig-7.3.0"
            / "MatrixElement"
            / "DIS"
            / "DISBase.cc"
        ).read_text()

        self.assertIn(
            "return sqr(scaleFact_) * (useQ2ScaleInPOWHEGEmission_ ? q2 : 0.25*q2*sqr(xT));",
            source,
        )
        self.assertIn(
            "return sqr(scaleFact_) * (useQ2ScaleInPOWHEGEmission_ ? q2 : mappedScale);",
            source,
        )
        self.assertIn(
            "Energy2 alphaScale = sqr(scaleFact_) * (useQ2ScaleInPOWHEGEmission_ ? q2_ : scale);",
            source,
        )
        self.assertIn("Energy2 bornPDFScale = sqr(scaleFact_) * q2_;", source)
        self.assertNotIn("pdfScale, q2_,", source)
        self.assertNotIn("Energy2 alphaScale = useQ2ScaleInPOWHEGEmission_ ? q2_ : scale;", source)


if __name__ == "__main__":
    unittest.main()
