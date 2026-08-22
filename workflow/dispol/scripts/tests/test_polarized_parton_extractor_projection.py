import math
import re
import unittest
from pathlib import Path


def find_source_layout():
    for root in Path(__file__).resolve().parents:
        source = root / "src" / "thepeg"
        if source.exists():
            return source
    raise RuntimeError("Could not locate the ThePEG polarized extractor source")


THEPEG_ROOT = find_source_layout()
EXTRACTOR_CC = THEPEG_ROOT / "PDF" / "PolarizedPartonExtractor.cc"
EXTRACTOR_H = THEPEG_ROOT / "PDF" / "PolarizedPartonExtractor.h"


def project(longitudinal, transverse):
    components = (longitudinal, transverse.real, transverse.imag)
    if not all(math.isfinite(value) for value in components):
        raise ValueError("nonfinite polarization")
    scale = max(abs(value) for value in components)
    if scale == 0.0:
        return longitudinal, transverse
    scaled_norm = math.sqrt(sum((value / scale) ** 2 for value in components))
    if scale <= 1.0 / scaled_norm:
        return longitudinal, transverse
    factor = (1.0 / scale) / scaled_norm
    return longitudinal * factor, transverse * factor


class PolarizedPartonExtractorProjectionTests(unittest.TestCase):

    def test_inside_ball_is_unchanged(self):
        value = (0.3, complex(0.4, -0.2))
        self.assertEqual(project(*value), value)

    def test_full_vector_is_projected_radially(self):
        longitudinal, transverse = project(2.0, complex(1.0, -2.0))
        norm = math.sqrt(
            longitudinal * longitudinal
            + transverse.real * transverse.real
            + transverse.imag * transverse.imag
        )
        self.assertAlmostEqual(norm, 1.0, places=14)
        self.assertAlmostEqual(transverse.real / longitudinal, 0.5, places=14)
        self.assertAlmostEqual(transverse.imag / longitudinal, -1.0, places=14)

    def test_projection_is_stable_for_large_finite_components(self):
        longitudinal, transverse = project(1.0e308, complex(-1.0e308, 1.0e308))
        self.assertTrue(math.isfinite(longitudinal))
        self.assertTrue(math.isfinite(transverse.real))
        self.assertTrue(math.isfinite(transverse.imag))
        self.assertAlmostEqual(
            longitudinal * longitudinal
            + transverse.real * transverse.real
            + transverse.imag * transverse.imag,
            1.0,
            places=14,
        )

    def test_nonfinite_components_are_rejected(self):
        for value in (
            (math.nan, 0j),
            (0.0, complex(math.inf, 0.0)),
            (0.0, complex(0.0, -math.inf)),
        ):
            with self.subTest(value=value), self.assertRaises(ValueError):
                project(*value)

    def test_cpp_applies_projection_after_pdf_reweighting_on_both_beams(self):
        source = re.sub(r"\s+", "", EXTRACTOR_CC.read_text())
        self.assertEqual(
            source.count("constautophysical=projectPhysicalPolarization(pL,pT);"),
            2,
        )
        self.assertEqual(source.count("pL*=diff/sum;"), 2)
        self.assertLess(
            source.index("pL*=diff/sum;"),
            source.index("constautophysical=projectPhysicalPolarization(pL,pT);"),
        )
        self.assertIn("if(scale>1./scaledMagnitude)", source)
        self.assertIn("!std::isfinite(longitudinal)", source)

    def test_public_helper_contract_is_declared(self):
        header = re.sub(r"\s+", "", EXTRACTOR_H.read_text())
        self.assertIn(
            "staticpair<double,Complex>projectPhysicalPolarization"
            "(doublelongitudinal,Complextransverse);",
            header,
        )


if __name__ == "__main__":
    unittest.main()
