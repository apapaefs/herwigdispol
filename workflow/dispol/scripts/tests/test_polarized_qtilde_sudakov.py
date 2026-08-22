import math
import re
import unittest
from pathlib import Path


CPP_TARGETS = (
    Path("Shower/QTilde/SplittingFunctions/SudakovFormFactor.h"),
    Path("Shower/QTilde/SplittingFunctions/SudakovFormFactor.cc"),
    Path("Shower/QTilde/SplittingFunctions/Sudakov1to2FormFactor.cc"),
    Path("Shower/QTilde/Kinematics/IS_QTildeShowerKinematics1to2.cc"),
)


def find_source_layout():
    standalone = None
    for root in Path(__file__).resolve().parents:
        active = root / "HerwigSource" / "Herwig-7.3.0"
        mirror = root / "herwigdispol" / "src" / "herwig"
        if active.exists():
            return active, mirror if mirror.exists() else None

        curated = root / "src" / "herwig"
        if curated.exists() and standalone is None:
            standalone = curated

    if standalone is not None:
        return standalone, None
    raise RuntimeError("Could not locate the Herwig QTilde source layout")


HERWIG_ROOT, MIRROR_ROOT = find_source_layout()
SUDAKOV_H = HERWIG_ROOT / CPP_TARGETS[0]
SUDAKOV_CC = HERWIG_ROOT / CPP_TARGETS[1]
SUDAKOV_1TO2_CC = HERWIG_ROOT / CPP_TARGETS[2]
IS_KINEMATICS_CC = HERWIG_ROOT / CPP_TARGETS[3]


def compact(source):
    return re.sub(r"\s+", "", source)


def source_between(source, start, end):
    begin = source.index(start)
    finish = source.index(end, begin)
    return source[begin:finish]


def density_contraction(h_matrix, rho_matrix):
    return sum(
        h_matrix[row][column] * rho_matrix[row][column]
        for row in range(len(h_matrix))
        for column in range(len(h_matrix[row]))
    )


def normalized_polarized_factor(numerator, denominator):
    den_real = denominator.real
    den_imag = denominator.imag
    reality_tolerance = 1.0e-10 * den_real
    if (
        not math.isfinite(numerator)
        or not math.isfinite(den_real)
        or not math.isfinite(den_imag)
        or abs(den_imag) > reality_tolerance
        or den_real <= 0.0
        or numerator < -1.0e-12
    ):
        raise ValueError("invalid polarized PDF normalization")
    physical_numerator = max(0.0, numerator)
    factor = physical_numerator / (2.0 * den_real)
    if not math.isfinite(factor) or factor < 0.0:
        raise ValueError("invalid polarized PDF ratio")
    return factor


class PolarizedQTildeSudakovTests(unittest.TestCase):

    def test_unpolarized_density_matrix_closes_to_unit_weight(self):
        h_matrix = ((0.5, 0.0), (0.0, 0.5))
        rho_matrix = (
            (0.7, complex(0.02, -0.04)),
            (complex(0.02, 0.04), 0.3),
        )

        denominator = density_contraction(h_matrix, rho_matrix)

        self.assertAlmostEqual(denominator.real, 0.5, places=14)
        self.assertAlmostEqual(
            normalized_polarized_factor(1.0, denominator),
            1.0,
            places=14,
        )

    def test_nontrivial_hermitian_density_contraction_is_normalized(self):
        h_matrix = (
            (0.4, complex(0.1, 0.05)),
            (complex(0.1, -0.05), 0.6),
        )
        rho_matrix = (
            (0.7, complex(0.02, -0.04)),
            (complex(0.02, 0.04), 0.3),
        )

        denominator = density_contraction(h_matrix, rho_matrix)

        self.assertAlmostEqual(denominator.real, 0.468, places=14)
        self.assertAlmostEqual(denominator.imag, 0.0, places=14)
        self.assertAlmostEqual(
            normalized_polarized_factor(0.936, denominator),
            1.0,
            places=14,
        )

    def test_normalized_factor_divides_by_old_density_trace(self):
        numerator = 0.84
        denominator = complex(0.35, 5.0e-12)

        factor = normalized_polarized_factor(numerator, denominator)

        self.assertAlmostEqual(factor, 1.2, places=14)
        self.assertNotAlmostEqual(factor, numerator, places=14)

    def test_normalized_factor_rejects_invalid_denominators_and_weights(self):
        invalid_inputs = (
            (math.nan, complex(0.5, 0.0)),
            (1.0, complex(math.nan, 0.0)),
            (1.0, complex(math.inf, 0.0)),
            (1.0, complex(0.5, math.nan)),
            (1.0, complex(0.5, math.inf)),
            (1.0, complex(0.5, 6.0e-11)),
            (1.0, complex(0.0, 0.0)),
            (1.0, complex(-0.5, 0.0)),
            (-0.1, complex(0.5, 0.0)),
            (1.0e308, complex(1.0e-308, 0.0)),
        )

        for numerator, denominator in invalid_inputs:
            with self.subTest(numerator=numerator, denominator=denominator):
                with self.assertRaises(ValueError):
                    normalized_polarized_factor(numerator, denominator)

        self.assertEqual(
            normalized_polarized_factor(-5.0e-13, complex(0.5, 0.0)),
            0.0,
        )

    def test_pdf_veto_uses_guarded_normalized_polarized_factor(self):
        source = SUDAKOV_CC.read_text()
        body = source_between(
            source,
            "double SudakovFormFactor::PDFVetoRatio",
            "bool SudakovFormFactor::alphaSVeto",
        )
        code = compact(body)

        self.assertIn(
            "constdoublepolarizedRatio=physicalNumerator/(2.*denReal);",
            code,
        )
        self.assertIn(
            "constdoublephysicalNumerator=numerator<0.?0.:numerator;",
            code,
        )
        self.assertIn("ratio*=polarizedRatio;", code)
        self.assertNotIn(
            "ratio*=initialStatePolarizationWeight(z,HNew,rho);",
            code,
        )
        for guard in (
            "!std::isfinite(numerator)",
            "!std::isfinite(denReal)",
            "!std::isfinite(denImag)",
            "denReal<=0.",
            "abs(denImag)>realityTolerance",
            "numerator<-1e-12",
            "!std::isfinite(polarizedRatio)",
            "polarizedRatio<0.",
            "!std::isfinite(ratio)",
            "ratio<0.",
        ):
            with self.subTest(guard=guard):
                self.assertIn(guard, code)

    def test_pdf_scale_is_centralized_and_used_by_pdf_veto(self):
        header = SUDAKOV_H.read_text()
        source = SUDAKOV_CC.read_text()
        helper = source_between(
            source,
            "Energy2 SudakovFormFactor::effectivePDFScale",
            "double SudakovFormFactor::PDFVetoRatio",
        )
        veto = source_between(
            source,
            "double SudakovFormFactor::PDFVetoRatio",
            "bool SudakovFormFactor::alphaSVeto",
        )

        self.assertRegex(
            header,
            r"Energy2\s+effectivePDFScale\s*\(\s*Energy2\s+\w+\s*,"
            r"\s*double\s+\w+\s*=\s*1(?:\.0?)?\s*\)\s*const\s*;",
        )
        self.assertIn("factorizationScaleFactor()", helper)
        self.assertIn("sqr(freeze_)", helper)
        self.assertIn("max(", helper)
        self.assertIn(
            "constEnergy2theScale=effectivePDFScale(t,factor);",
            compact(veto),
        )

    def test_backward_azimuth_uses_parent_x_and_pdf_scale_not_split_invariant(self):
        source = SUDAKOV_1TO2_CC.read_text()
        body = source_between(
            source,
            "double Sudakov1to2FormFactor::generatePhiBackward",
            "double Sudakov1to2FormFactor::generatePhiDecay",
        )
        code = compact(body)

        self.assertIn(
            "constEnergy2qtilde2=sqr(kinematics->scale());",
            code,
        )
        self.assertIn("Energy2t=(1.-z)*qtilde2/z;", code)
        self.assertIn("backwardPhiWeights(z,t,ids,rho,H);", code)
        self.assertIn(
            "calculateHMatrix(beam,ids[0],x/z,"
            "effectivePDFScale(qtilde2));",
            code,
        )
        self.assertNotIn("calculateHMatrix(beam,ids[0],x,t)", code)

    def test_terminal_density_update_uses_effective_pdf_scale(self):
        source = IS_KINEMATICS_CC.read_text()
        body = source_between(
            source,
            "void IS_QTildeShowerKinematics1to2::\nupdateLast",
            "void IS_QTildeShowerKinematics1to2::\nresetChildren",
        )
        code = compact(body)

        self.assertIn(
            "calculateHMatrix(beam,theLast->dataPtr(),theLast->x(),"
            "SudakovFormFactor()->effectivePDFScale(sqr(scale())))",
            code,
        )
        self.assertNotIn(
            "calculateHMatrix(beam,theLast->dataPtr(),theLast->x(),"
            "sqr(scale()))",
            code,
        )

    def test_active_and_curated_qtilde_sources_are_byte_identical(self):
        if MIRROR_ROOT is None:
            self.skipTest("standalone source layout has no active-tree mirror")

        for relative in CPP_TARGETS:
            with self.subTest(source=str(relative)):
                self.assertEqual(
                    (HERWIG_ROOT / relative).read_bytes(),
                    (MIRROR_ROOT / relative).read_bytes(),
                )


if __name__ == "__main__":
    unittest.main()
