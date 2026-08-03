import importlib.util
import math
import unittest
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "make_spinphysical_paper_plots.py"
SPEC = importlib.util.spec_from_file_location("make_spinphysical_paper_plots", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SpinphysicalPaperPlotTests(unittest.TestCase):
    def test_publication_curves_use_markers_without_connecting_lines(self):
        class RecordingAxis:
            def errorbar(self, *_args, **kwargs):
                self.kwargs = kwargs

        axis = RecordingAxis()
        curve = {
            "x": np.asarray([0.5]),
            "values": np.asarray([1.0]),
            "errors": np.asarray([0.1]),
        }

        MODULE._draw_curve(axis, curve, "RIVETPS-SPIN")

        self.assertEqual(axis.kwargs["linestyle"], "none")

    def test_recovered_tag_shard_overrides_stale_spec_index(self):
        spec = {"tag": "spinphysical15-replacement-s212", "shard_index": 1}
        self.assertEqual(MODULE.recover_shard_index(spec), 212)

    def test_ratio_error_retains_numerator_denominator_covariance(self):
        numerator = np.asarray([[2.0, 4.0], [4.0, 8.0], [6.0, 12.0]])
        denominator = np.asarray([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]])
        ratio, error, influence = MODULE.ratio_of_means(numerator, denominator)
        np.testing.assert_allclose(ratio, [2.0, 2.0])
        np.testing.assert_allclose(influence, 0.0, atol=1e-15)
        np.testing.assert_allclose(error, 0.0, atol=1e-15)

    def test_rebin_integrates_histogram_heights(self):
        blocks = np.asarray([[1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0, 8.0]])
        old_edges = np.asarray([0.0, 0.25, 0.5, 0.75, 1.0])
        new_edges = np.asarray([0.0, 0.5, 1.0])
        rebinned = MODULE.integrate_rebin(blocks, old_edges, new_edges)
        np.testing.assert_allclose(rebinned, [[0.75, 1.75], [1.5, 3.5]])

    def test_cumulative_tail_is_width_weighted(self):
        blocks = np.asarray([[1.0, 2.0, 3.0]])
        edges = np.asarray([0.0, 0.2, 0.5, 1.0])
        np.testing.assert_allclose(MODULE.cumulative_tail(blocks, edges), [[2.3, 2.1, 1.5]])

    def test_independent_ratio_error(self):
        ratio, error = MODULE.independent_ratio(
            np.asarray([2.0]),
            np.asarray([0.2]),
            np.asarray([4.0]),
            np.asarray([0.4]),
        )
        self.assertAlmostEqual(float(ratio[0]), 0.5)
        self.assertAlmostEqual(float(error[0]), math.hypot(0.05, 0.05))

    def test_reference_display_label_has_no_empty_braces(self):
        self.assertEqual(MODULE.LABELS[MODULE.REFERENCE_VARIANT], "None")
        self.assertNotIn("{}", MODULE.LABELS[MODULE.REFERENCE_VARIANT])


if __name__ == "__main__":
    unittest.main()
