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


class CampaignPythonInterpreterTests(unittest.TestCase):

    def test_campaign_helpers_use_current_python_interpreter(self):
        text = (scripts_dir() / "run_validation_campaign.py").read_text()

        self.assertNotIn("python3.10", text)
        self.assertIn("sys.executable", text)


if __name__ == "__main__":
    unittest.main()
