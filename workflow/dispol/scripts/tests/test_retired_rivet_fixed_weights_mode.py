import unittest
from pathlib import Path


def layout_candidates(root):
    return (
        (
            root / "DISPOL" / "scripts" / "run_validation_campaign.py",
            root / "HerwigSource" / "Herwig-7.3.0" / "MatrixElement" / "DIS" / "DISBase.cc",
            root / "HerwigSource" / "Herwig-7.3.0" / "MatrixElement" / "DIS" / "DISBase.h",
        ),
        (
            root / "workflow" / "dispol" / "scripts" / "run_validation_campaign.py",
            root / "src" / "herwig" / "MatrixElement" / "DIS" / "DISBase.cc",
            root / "src" / "herwig" / "MatrixElement" / "DIS" / "DISBase.h",
        ),
    )


def find_layout():
    for root in Path(__file__).resolve().parents:
        for campaign, disbase_cc, disbase_h in layout_candidates(root):
            if campaign.exists() and disbase_cc.exists() and disbase_h.exists():
                return campaign, disbase_cc, disbase_h
    raise RuntimeError("Could not locate a supported DISPOL/Herwig DIS layout")


CAMPAIGN, DISBASE_CC, DISBASE_H = find_layout()


class RetiredRivetFixedWeightsModeTests(unittest.TestCase):

    def test_workflow_no_longer_exposes_rivetfixedweights(self):
        campaign = CAMPAIGN.read_text()

        self.assertNotIn("--rivetfixedweights", campaign)
        self.assertNotIn("RIVETFIXEDWEIGHTS", campaign)
        self.assertIn("RIVETWEIGHTS_ANALYSIS_VARIANT", campaign)
        self.assertIn("RIVETWEIGHTS_SHOWER_ANALYSIS_VARIANT", campaign)
        self.assertIn('return ["--rivetweights"]', campaign)

    def test_herwig_no_longer_exposes_fixed_order_validation_comparison_mode(self):
        source = DISBASE_CC.read_text()
        header = DISBASE_H.read_text()

        self.assertNotIn("FixedOrderValidation", source)
        self.assertNotIn("FixedOrderValidation", header)
        self.assertIn("POWHEGEmissionComparisonModeDefault = 0", header)
        self.assertIn("POWHEGEmissionComparisonModeRealOnly = 1", header)
        self.assertIn('"Default"', source)
        self.assertIn('"RealOnly"', source)


if __name__ == "__main__":
    unittest.main()
