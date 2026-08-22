import importlib.util
import math
import sys
import types
import unittest
from pathlib import Path

import numpy as np


def scripts_dir() -> Path:
    for root in Path(__file__).resolve().parents:
        for candidate in (
            root / "DISPOL" / "scripts",
            root / "workflow" / "dispol" / "scripts",
        ):
            if (candidate / "analyze_DIS_polarized.py").exists():
                return candidate
    raise RuntimeError("Could not locate the DIS workflow scripts")


def load_rivet_postprocess_module():
    path = scripts_dir() / "rivet_scale_plot_postprocess.py"
    spec = importlib.util.spec_from_file_location("rivet_scale_plot_postprocess_for_test", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class FakeBin:
    def __init__(self, xlow: float, xhigh: float, value: float = 0.0, error: float = 0.0):
        self._xlow = xlow
        self._xhigh = xhigh
        self._value = value
        self._error = error

    def xMin(self):
        return self._xlow

    def xMax(self):
        return self._xhigh

    def val(self):
        return self._value

    def errAvg(self, *_args):
        return self._error


class FakeHist:
    def __init__(self, edges, values=None, errors=None):
        values = list(values if values is not None else [0.0] * (len(edges) - 1))
        errors = list(errors if errors is not None else [0.0] * (len(edges) - 1))
        self._edges = list(edges)
        self._bins = [
            FakeBin(self._edges[i], self._edges[i + 1], values[i], errors[i])
            for i in range(len(values))
        ]
        self.path = None

    def xEdges(self):
        return list(self._edges)

    def bins(self):
        return list(self._bins)

    def numBins(self):
        return len(self._bins)

    def bin(self, index):
        return self._bins[index - 1]

    def setPath(self, path):
        self.path = path

    def setAnnotation(self, *_args):
        pass


class FakePoint:
    def __init__(self, x: float = 0.0, y: float = 0.0, xerr: float = 0.0, yerr: float = 0.0):
        self._x = x
        self._y = y
        self._xerrs = (xerr, xerr)
        self._yerrs = (yerr, yerr)

    def setX(self, value):
        self._x = value

    def setY(self, value):
        self._y = value

    def setXErrs(self, minus, plus):
        self._xerrs = (minus, plus)

    def setYErrs(self, minus, plus):
        self._yerrs = (minus, plus)

    def setErrs(self, minus, plus):
        self.setYErrs(minus, plus)

    def x(self):
        return self._x

    def y(self):
        return self._y

    def xErrs(self):
        return self._xerrs

    def yErrs(self):
        return self._yerrs

    def xMin(self):
        return self._x - self._xerrs[0]

    def xMax(self):
        return self._x + self._xerrs[1]


class FakeScatter:
    def __init__(self, points=None):
        self._points = list(points or [])
        self.path = None
        self._annotations = {}

    def setPath(self, path):
        self.path = path

    def points(self):
        return list(self._points)

    def point(self, index):
        return self._points[index]

    def addPoint(self, point):
        self._points.append(point)

    def numPoints(self):
        return len(self._points)

    def setAnnotation(self, key, value):
        self._annotations[key] = value

    def annotation(self, key):
        return self._annotations[key]

    def annotations(self):
        return list(self._annotations)


class FakeRatioAxis:
    def __init__(self, limits):
        self._limits = limits

    def get_ylim(self):
        return self._limits


def load_analyze_module_with_yoda_stubs():
    script_dir = scripts_dir()
    path = script_dir / "analyze_DIS_polarized.py"

    previous_yoda = sys.modules.get("yoda")
    previous_poldis = sys.modules.get("poldis_top_to_yoda")
    sys.path.insert(0, str(script_dir))

    yoda_stub = types.ModuleType("yoda")
    poldis_stub = types.ModuleType("poldis_top_to_yoda")

    def new_binned_estimate(edges):
        return FakeHist(edges)

    def set_bin_val_err(binwrap, y, dy):
        binwrap._value = y
        binwrap._error = dy

    poldis_stub.new_binned_estimate = new_binned_estimate
    poldis_stub.set_bin_val_err = set_bin_val_err
    yoda_stub.Point2D = FakePoint
    yoda_stub.Scatter2D = FakeScatter
    sys.modules["yoda"] = yoda_stub
    sys.modules["poldis_top_to_yoda"] = poldis_stub

    try:
        spec = importlib.util.spec_from_file_location("analyze_DIS_polarized_for_test", path)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        return module
    finally:
        try:
            sys.path.remove(str(script_dir))
        except ValueError:
            pass
        if previous_yoda is None:
            sys.modules.pop("yoda", None)
        else:
            sys.modules["yoda"] = previous_yoda
        if previous_poldis is None:
            sys.modules.pop("poldis_top_to_yoda", None)
        else:
            sys.modules["poldis_top_to_yoda"] = previous_poldis


class SpinComparisonObservableTests(unittest.TestCase):

    def test_rivetweights_allow_list_keeps_hard_j3_and_fourth_jet_inputs(self):
        text = (scripts_dir() / "run_validation_campaign.py").read_text()
        required_labels = {
            "DeltaPhiHardJ3",
            "DDeltaPhiHardJ3",
            "Cos2DeltaPhiHardJ3NumPt3OverpT1",
            "Cos2DeltaPhiHardJ3DenPt3OverpT1",
            "DCos2DeltaPhiHardJ3NumPt3OverpT1",
            "pT4",
            "pT4OverpT1",
            "DeltaPhiJ13J14",
            "DDeltaPhiJ13J14",
            "Cos2DeltaPhiJ13J14NumPt4OverpT1",
            "Cos2DeltaPhiJ13J14DenPt4OverpT1",
            "DCos2DeltaPhiJ13J14NumPt4OverpT1",
        }

        for label in required_labels:
            self.assertIn(f'"{label}"', text)

    def test_no_scale_ratio_reference_accepts_rivetpsweights_none_labels(self):
        module = load_rivet_postprocess_module()
        full = "sample.RIVETPSWEIGHTS-SPIN.sanitized.yoda.gz"
        born = "sample.RIVETPSWEIGHTS-NOSPIN.sanitized.yoda.gz"
        none = "sample.RIVETPSWEIGHTS-NOSPIN-UNPOL.sanitized.yoda.gz"
        namespace = {
            "np": np,
            "dataf": {
                "add_legend_handle": [full, born, none],
                "yvals": {
                    full: [2.0, 4.0],
                    born: [1.5, 3.0],
                    none: [1.0, 2.0],
                },
                "yerrs": {
                    full: [[0.0, 0.0], [0.0, 0.0]],
                    born: [[0.0, 0.0], [0.0, 0.0]],
                    none: [[0.0, 0.0], [0.0, 0.0]],
                },
            },
        }

        exec(module._NOSCALE_HELPER_BLOCK, namespace)

        self.assertEqual(namespace["ratio_reference_label"], none)
        np.testing.assert_allclose(namespace["_ratio_yvals"](full), [2.0, 2.0])
        np.testing.assert_allclose(namespace["_ratio_yvals"](none), [1.0, 1.0])

    def test_no_scale_ratio_patch_preserves_generated_ratio_limits(self):
        module = load_rivet_postprocess_module()
        namespace = {
            "np": np,
            "ratio0_ax": FakeRatioAxis((0.94, 1.02)),
            "dataf": {
                "add_legend_handle": ["Full", "None"],
                "yvals": {"Full": [0.97], "None": [1.0]},
                "yerrs": {
                    "Full": [[0.0], [0.0]],
                    "None": [[0.0], [0.0]],
                },
            },
        }

        exec(module._NOSCALE_HELPER_BLOCK, namespace)

        self.assertEqual(namespace["_no_scale_ratio_y_lims"](), (0.94, 1.02))

    def test_scale_envelope_patch_detects_plot_safe_standalone_and_poldis_labels(self):
        module = load_rivet_postprocess_module()
        herwig = "standalone_fo_all_correlated_q2gt100.scale_envelope_plot.yoda"
        poldis = "poldis_reference_all_q2gt100.plot_safe.yoda"
        namespace = {
            "np": np,
            "dataf": {
                "add_legend_handle": [herwig, poldis],
                "yvals": {
                    herwig: [2.0, 4.0],
                    poldis: [1.0, 2.0],
                },
                "yerrs": {
                    herwig: [[0.2, 0.4], [0.5, 0.7]],
                    poldis: [[0.0, 0.0], [0.3, 0.4]],
                },
                "xpoints": {
                    herwig: [0.5, 1.5],
                    poldis: [0.5, 1.5],
                },
                "xedges": {
                    herwig: [0.0, 1.0, 2.0],
                    poldis: [0.0, 1.0, 2.0],
                },
                "xerrs": {
                    herwig: [[0.5, 0.5], [0.5, 0.5]],
                    poldis: [[0.5, 0.5], [0.5, 0.5]],
                },
                "variation_yvals": {},
            },
        }

        exec(module.scale_envelope_helper_block("POLDIS NLO"), namespace)
        namespace["styles"] = {
            poldis: {"color": "#EE3311", "histstyle": 1, "marker": "none", "zorder": 5},
            herwig: {"color": "#3366FF", "histstyle": 1, "marker": "none", "zorder": 6},
        }
        namespace["_ensure_scale_patch_styles"]()

        self.assertEqual(namespace["reference_label"], "POLDIS NLO")
        self.assertEqual(namespace["scale_band_labels"], {herwig})
        self.assertIn("POLDIS NLO", namespace["dataf"]["yvals"])
        self.assertNotIn(poldis, namespace["dataf"]["yvals"])
        np.testing.assert_allclose(
            namespace["dataf"]["yerrs"]["POLDIS NLO"],
            [[0.3, 0.4], [0.3, 0.4]],
        )
        self.assertEqual(namespace["styles"]["POLDIS NLO"]["color"], "black")
        self.assertEqual(namespace["styles"]["POLDIS NLO"]["histstyle"], 0)
        self.assertEqual(namespace["styles"]["POLDIS NLO"]["marker"], "o")
        self.assertEqual(namespace["styles"][herwig]["color"], "#009E73")
        self.assertEqual(namespace["styles"][herwig]["ratio0_1sigmabandcolor"], "#009E73")
        self.assertEqual(namespace["styles"][herwig]["histstyle"], 1)

    def test_scale_envelope_patch_preserves_signed_scale_band_errors(self):
        module = load_rivet_postprocess_module()
        herwig = "standalone_fo_all_correlated_q2gt100.scale_envelope_plot.yoda"
        namespace = {
            "np": np,
            "dataf": {
                "add_legend_handle": ["POLDIS NLO", herwig],
                "yvals": {
                    "POLDIS NLO": [1.0, 2.0],
                    herwig: [-1.0, 2.0],
                },
                "yerrs": {
                    "POLDIS NLO": [[0.1, 0.2], [0.1, 0.2]],
                    herwig: [[0.4, 0.5], [0.2, 0.3]],
                },
                "variation_yvals": {},
            },
        }

        exec(module.scale_envelope_helper_block("POLDIS NLO"), namespace)
        namespace["styles"] = {
            "POLDIS NLO": {"drawstyle": None, "ratio0_drawstyle": None},
            herwig: {"drawstyle": None, "ratio0_drawstyle": None},
        }

        band_dn, band_up = namespace["_main_band_arrays"](herwig)
        ratio_dn, ratio_up = namespace["_ratio_band_arrays"](herwig)

        np.testing.assert_allclose(band_dn, [-1.4, 1.5])
        np.testing.assert_allclose(band_up, [-0.8, 2.3])
        np.testing.assert_allclose(ratio_dn, [-1.4, 0.75])
        np.testing.assert_allclose(ratio_up, [-0.8, 1.15])

    def test_dzeta_scale_envelope_legend_moves_to_lower_center(self):
        module = load_rivet_postprocess_module()
        source = "\n".join(
            [
                "#! /usr/bin/env python3",
                "# Analysis object: /MC_DIS_BREIT/DZeta",
                "dataf = dict()",
                "exec(open(prefix+'DZeta__data.py').read(), dataf)",
                "ratio0_ax.set_ylim(0.5, 1.5)",
                "ratio0_ax.set_ylabel('Ratio')",
                "# curve from input yoda files in main panel",
                "old_main_block()",
                "# plots on ratio panel",
                "old_ratio_block()",
                "legend_items = list(legend_handles.values())",
                "ax.add_artist(ax.legend(legend_items,",
                "                        labels['legend'][0],",
                "                        title=labels['legend'][1],",
                "                        alignment='left',",
                "                        loc='upper left',",
                "                        bbox_to_anchor=(0.0, 0.97)))",
            ]
        )

        patched = module.patch_scale_envelope_plot_script_text(source)

        self.assertIn("loc='lower center'", patched)
        self.assertIn("bbox_to_anchor=(0.5, 0.03)", patched)

    def test_zeta_scale_envelope_legend_moves_to_lower_center(self):
        module = load_rivet_postprocess_module()
        source = "\n".join(
            [
                "#! /usr/bin/env python3",
                "# Analysis object: /MC_DIS_BREIT/Zeta",
                "dataf = dict()",
                "exec(open(prefix+'Zeta__data.py').read(), dataf)",
                "ratio0_ax.set_ylim(0.5, 1.5)",
                "ratio0_ax.set_ylabel('Ratio')",
                "# curve from input yoda files in main panel",
                "old_main_block()",
                "# plots on ratio panel",
                "old_ratio_block()",
                "legend_items = list(legend_handles.values())",
                "ax.add_artist(ax.legend(legend_items,",
                "                        labels['legend'][0],",
                "                        title=labels['legend'][1],",
                "                        alignment='left',",
                "                        loc='upper left',",
                "                        bbox_to_anchor=(0.0, 0.97)))",
            ]
        )

        patched = module.patch_scale_envelope_plot_script_text(source)

        self.assertIn("loc='lower center'", patched)
        self.assertIn("bbox_to_anchor=(0.5, 0.03)", patched)

    def test_builds_cumulative_all_pt3_over_pt1_from_tail_integrals(self):
        module = load_analyze_module_with_yoda_stubs()
        objects = {
            "/MC_DIS_PS/pT3OverpT1": FakeHist([0.0, 0.5, 1.0], [2.0, 4.0], [0.2, 0.4]),
            "/MC_DIS_PS/DpT3OverpT1": FakeHist([0.0, 0.5, 1.0], [1.0, 3.0], [0.1, 0.3]),
        }

        derived = module.build_correlated_helicity_derived_objects(objects, "MC_DIS_PS")

        hist = derived["/MC_DIS_PS/ALLpT3OverpT1Cumulative"]
        self.assertTrue(math.isclose(hist.bin(1).val(), 2.0 / 3.0))
        self.assertTrue(math.isclose(hist.bin(2).val(), 3.0 / 4.0))

    def test_plot_safe_derives_missing_precut_all_ratios(self):
        module = load_analyze_module_with_yoda_stubs()
        objects = {
            "/REF/MC_DIS_BREIT/Q2PreCut": FakeScatter([
                FakePoint(150.0, 10.0, 50.0, 1.0),
                FakePoint(300.0, 20.0, 100.0, 2.0),
            ]),
            "/REF/MC_DIS_BREIT/DQ2PreCut": FakeScatter([
                FakePoint(150.0, 2.0, 50.0, 0.2),
                FakePoint(300.0, 6.0, 100.0, 0.6),
            ]),
        }

        derived = module.add_missing_all_ratio_objects(objects, "MC_DIS_BREIT")

        scatter = derived["/REF/MC_DIS_BREIT/ALLQ2PreCut"]
        self.assertTrue(math.isclose(scatter.point(0).y(), 0.2))
        self.assertTrue(math.isclose(scatter.point(1).y(), 0.3))
        self.assertIn("/REF/MC_DIS_BREIT/Q2PreCut", derived)
        self.assertIn("/REF/MC_DIS_BREIT/DQ2PreCut", derived)

    def test_plot_overlay_path_strips_reference_prefix(self):
        module = load_analyze_module_with_yoda_stubs()

        self.assertEqual(
            "/MC_DIS_BREIT/ALLQ2PreCut",
            module._plot_overlay_object_path("/REF/MC_DIS_BREIT/ALLQ2PreCut"),
        )
        self.assertEqual(
            "/MC_DIS_BREIT/ALLQ2PreCut",
            module._plot_overlay_object_path("/MC_DIS_BREIT/ALLQ2PreCut"),
        )

    def test_builds_cumulative_all_pt4_over_pt1_from_tail_integrals(self):
        module = load_analyze_module_with_yoda_stubs()
        objects = {
            "/MC_DIS_PS/pT4OverpT1": FakeHist([0.0, 0.5, 1.0], [3.0, 5.0], [0.3, 0.5]),
            "/MC_DIS_PS/DpT4OverpT1": FakeHist([0.0, 0.5, 1.0], [1.0, 2.0], [0.1, 0.2]),
        }

        derived = module.build_correlated_helicity_derived_objects(objects, "MC_DIS_PS")

        hist = derived["/MC_DIS_PS/ALLpT4OverpT1Cumulative"]
        self.assertTrue(math.isclose(hist.bin(1).val(), 3.0 / 8.0))
        self.assertTrue(math.isclose(hist.bin(2).val(), 2.0 / 5.0))

    def test_builds_all_delta_phi_hard_j3_from_matching_delta_histogram(self):
        module = load_analyze_module_with_yoda_stubs()
        objects = {
            "/MC_DIS_PS/DeltaPhiHardJ3": FakeHist([0.0, 1.0, 2.0], [2.0, 4.0], [0.2, 0.4]),
            "/MC_DIS_PS/DDeltaPhiHardJ3": FakeHist([0.0, 1.0, 2.0], [1.0, 1.0], [0.1, 0.1]),
        }

        derived = module.build_correlated_helicity_derived_objects(objects, "MC_DIS_PS")

        hist = derived["/MC_DIS_PS/ALLDeltaPhiHardJ3"]
        self.assertTrue(math.isclose(hist.bin(1).val(), 0.5))
        self.assertTrue(math.isclose(hist.bin(2).val(), 0.25))

    def test_builds_all_delta_phi_j13_j14_from_matching_delta_histogram(self):
        module = load_analyze_module_with_yoda_stubs()
        objects = {
            "/MC_DIS_PS/DeltaPhiJ13J14": FakeHist([0.0, 1.0, 2.0], [4.0, 6.0], [0.4, 0.6]),
            "/MC_DIS_PS/DDeltaPhiJ13J14": FakeHist([0.0, 1.0, 2.0], [1.0, 3.0], [0.1, 0.3]),
        }

        derived = module.build_correlated_helicity_derived_objects(objects, "MC_DIS_PS")

        hist = derived["/MC_DIS_PS/ALLDeltaPhiJ13J14"]
        self.assertTrue(math.isclose(hist.bin(1).val(), 0.25))
        self.assertTrue(math.isclose(hist.bin(2).val(), 0.5))

    def test_builds_cumulative_delta_phi_hard_j3_cos2_moments_from_tail_integrals(self):
        module = load_analyze_module_with_yoda_stubs()
        objects = {
            "/MC_DIS_PS/Cos2DeltaPhiHardJ3NumPt3OverpT1": FakeHist(
                [0.0, 0.5, 1.0],
                [1.0, -2.0],
                [0.1, 0.2],
            ),
            "/MC_DIS_PS/Cos2DeltaPhiHardJ3DenPt3OverpT1": FakeHist(
                [0.0, 0.5, 1.0],
                [2.0, 4.0],
                [0.2, 0.4],
            ),
            "/MC_DIS_PS/DCos2DeltaPhiHardJ3NumPt3OverpT1": FakeHist(
                [0.0, 0.5, 1.0],
                [0.5, -1.0],
                [0.05, 0.1],
            ),
        }

        derived = module.build_correlated_helicity_derived_objects(objects, "MC_DIS_PS")

        cos2_hist = derived["/MC_DIS_PS/Cos2DeltaPhiHardJ3Cumulative"]
        self.assertTrue(math.isclose(cos2_hist.bin(1).val(), -1.0 / 6.0))
        self.assertTrue(math.isclose(cos2_hist.bin(2).val(), -0.5))

        all_cos2_hist = derived["/MC_DIS_PS/ALLCos2DeltaPhiHardJ3Cumulative"]
        self.assertTrue(math.isclose(all_cos2_hist.bin(1).val(), -1.0 / 12.0))
        self.assertTrue(math.isclose(all_cos2_hist.bin(2).val(), -0.25))

    def test_builds_cumulative_delta_phi_j13_j14_cos2_moments_from_tail_integrals(self):
        module = load_analyze_module_with_yoda_stubs()
        objects = {
            "/MC_DIS_PS/Cos2DeltaPhiJ13J14NumPt4OverpT1": FakeHist(
                [0.0, 0.5, 1.0],
                [2.0, -1.0],
                [0.2, 0.1],
            ),
            "/MC_DIS_PS/Cos2DeltaPhiJ13J14DenPt4OverpT1": FakeHist(
                [0.0, 0.5, 1.0],
                [4.0, 2.0],
                [0.4, 0.2],
            ),
            "/MC_DIS_PS/DCos2DeltaPhiJ13J14NumPt4OverpT1": FakeHist(
                [0.0, 0.5, 1.0],
                [1.0, -0.5],
                [0.1, 0.05],
            ),
        }

        derived = module.build_correlated_helicity_derived_objects(objects, "MC_DIS_PS")

        cos2_hist = derived["/MC_DIS_PS/Cos2DeltaPhiJ13J14Cumulative"]
        self.assertTrue(math.isclose(cos2_hist.bin(1).val(), 1.0 / 6.0))
        self.assertTrue(math.isclose(cos2_hist.bin(2).val(), -0.5))

        all_cos2_hist = derived["/MC_DIS_PS/ALLCos2DeltaPhiJ13J14Cumulative"]
        self.assertTrue(math.isclose(all_cos2_hist.bin(1).val(), 1.0 / 12.0))
        self.assertTrue(math.isclose(all_cos2_hist.bin(2).val(), -0.25))


if __name__ == "__main__":
    unittest.main()
