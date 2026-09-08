import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1] / "validation" / "compare_hard_spin_validation.py"
spec = importlib.util.spec_from_file_location("hard_spin_compare", SOURCE)
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)
relocate_spec = importlib.util.spec_from_file_location(
    "relocate_validation", SOURCE.with_name("relocate_validation_run.py"))
relocate = importlib.util.module_from_spec(relocate_spec)
relocate_spec.loader.exec_module(relocate)


class HardSpinValidationTests(unittest.TestCase):
    def fixtures(self, directory, changed=False, compact=False):
        paths = [Path(directory)/name for name in ("a.jsonl","b.jsonl")]
        for index, path in enumerate(paths):
            with path.open("w") as stream:
                for event in range(1000):
                    hard = [[21,501,502,0.,0.,10.+event,10.+event,0.]]
                    if changed and index and event == 0:
                        hard[0][0] = 2
                    json.dump({"hard":hard,"weight":1,"features":[event%2]*35,
                        "final_state":[] if compact else [str(event)],
                        "hard_links":0,"shower_links":2}, stream)
                    stream.write("\n")
        return paths

    def test_paired_covariance_cancels_identical_shower_statistics(self):
        with tempfile.TemporaryDirectory() as directory:
            result = compare.compare(*self.fixtures(directory))
            self.assertTrue(result["passed"])
            self.assertEqual(result["observables"][0]["paired_standard_error"],0.)
            self.assertGreater(result["observables"][0]["shared_hard_covariance_of_means"],0.)

    def test_hard_identity_mismatch_fails_even_if_rates_agree(self):
        with tempfile.TemporaryDirectory() as directory:
            result = compare.compare(*self.fixtures(directory,changed=True))
            self.assertFalse(result["passed"])
            self.assertEqual(result["hard_identity_mismatches"],1)

    def test_compact_records_cannot_pass_bitwise_regression(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError,"compact"):
                compare.compare(*self.fixtures(directory,compact=True),regression=True)

    def test_compact_closure_does_not_claim_identical_final_states(self):
        with tempfile.TemporaryDirectory() as directory:
            result = compare.compare(*self.fixtures(directory,compact=True),by_identity=True)
            self.assertIsNone(result["unequal_final_states"])
            self.assertFalse(result["final_state_comparison_available"])

    def test_vetoed_hard_inputs_count_as_zero_not_dropped_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = self.fixtures(directory)
            lines = second.read_text().splitlines(keepends=True)
            second.write_text("".join(lines[:1]+lines[2:]))
            result = compare.compare(first,second,by_identity=True)
            self.assertEqual(result["events"],1000)
            self.assertEqual(result["vetoed_hard_inputs"],1)
            self.assertAlmostEqual(result["observables"][0]["candidate"],0.499)
            self.assertAlmostEqual(result["observables"][0]["difference"],0.001)

    def test_duplicate_or_unknown_lhe_identity_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = self.fixtures(directory)
            lines = second.read_text().splitlines(keepends=True)
            second.write_text("".join(lines)+lines[0])
            with self.assertRaisesRegex(ValueError,"duplicate LHE"):
                compare.compare(first,second,by_identity=True)
            first, second = self.fixtures(directory,changed=True)
            with self.assertRaisesRegex(ValueError,"unknown or changed"):
                compare.compare(first,second,by_identity=True)

    def test_excessive_veto_fraction_fails_closure(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = self.fixtures(directory)
            second.write_text("".join(second.read_text().splitlines(keepends=True)[2:]))
            result = compare.compare(first,second,by_identity=True)
            self.assertFalse(result["passed"])
            self.assertEqual(result["vetoed_hard_inputs"],2)

    def test_saved_fixture_relocation_never_rewrites_object_payload(self):
        data = (b"ThePEG version 1 Database\n0\n3\n2\n/old/lib/Herwig\n/keep/lib\n"
                b"1\n/old/lib/ThePEG\n{ object payload /old/lib/Herwig\n}\n")
        relocated, report = relocate.relocate(data,"/old","/new")
        self.assertEqual(report["changed_loader_paths"],2)
        self.assertIn(b"/new/lib/Herwig",relocated)
        self.assertTrue(relocated.endswith(b"{ object payload /old/lib/Herwig\n}\n"))
        with self.assertRaisesRegex(ValueError,"unrecognized"):
            relocate.relocate(b"not a supported database", "/old", "/new")

    def test_hard_identity_tolerance_does_not_ignore_flavour(self):
        self.assertTrue(compare.same_hard([[21,1,2,1.]], [[21,1,2,1.+1.e-12]]))
        self.assertFalse(compare.same_hard([[21,1,2,1.]], [[2,1,2,1.]]))
        self.assertFalse(compare.same_hard([[21,1,2,1.]], [[21,1,2,1.,0.]]))


if __name__ == "__main__":
    unittest.main()
