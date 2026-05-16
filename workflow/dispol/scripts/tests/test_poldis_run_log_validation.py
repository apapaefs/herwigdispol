import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


def load_poldis_reference_module():
    for root in Path(__file__).resolve().parents:
        script = root / "DISPOL" / "scripts" / "run_poldis_window_reference.py"
        if script.exists():
            sys.path.insert(0, str(script.parent))
            spec = importlib.util.spec_from_file_location("run_poldis_window_reference", script)
            module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            return module
    raise RuntimeError("Could not locate DISPOL/scripts/run_poldis_window_reference.py")


class POLDISRunLogValidationTests(unittest.TestCase):

    def test_fatal_antikt_log_is_rejected_even_with_zero_exit_status(self):
        module = load_poldis_reference_module()
        with tempfile.TemporaryDirectory() as tmp:
            log_path = Path(tmp) / "run.log"
            log_path.write_text(
                "\n".join(
                    [
                        "$ ./poldis.x",
                        "LHAPDF 6.5.5 loading a/pol/pdf.dat",
                        " Fatal error in subroutine antikt",
                        " imerge=          -1",
                        "Thanks for using LHAPDF 6.5.5.",
                        "",
                    ]
                )
            )

            with self.assertRaisesRegex(RuntimeError, "Fatal error in subroutine antikt"):
                module.validate_poldis_run_log(log_path)

    def test_log_without_all_order_totals_is_rejected(self):
        module = load_poldis_reference_module()
        with tempfile.TemporaryDirectory() as tmp:
            log_path = Path(tmp) / "run.log"
            log_path.write_text(
                "\n".join(
                    [
                        "$ ./poldis.x",
                        "LO = 1.0 +- 0.1",
                        "",
                    ]
                )
            )

            with self.assertRaisesRegex(RuntimeError, "missing NLO, NNLO totals"):
                module.validate_poldis_run_log(log_path)

    def test_log_with_all_order_totals_is_accepted(self):
        module = load_poldis_reference_module()
        with tempfile.TemporaryDirectory() as tmp:
            log_path = Path(tmp) / "run.log"
            log_path.write_text(
                "\n".join(
                    [
                        "$ ./poldis.x",
                        "LO = 1.0 +- 0.1",
                        "NLO = 2.0 +- 0.2",
                        "NNLO = 3.0 +- 0.3",
                        "",
                    ]
                )
            )

            totals = module.validate_poldis_run_log(log_path)

        self.assertEqual(set(totals), {"LO", "NLO", "NNLO"})


if __name__ == "__main__":
    unittest.main()
