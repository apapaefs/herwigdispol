import unittest
from pathlib import Path


FORBIDDEN_SOURCE_TOKENS = (
    "DISDiagnostics",
    "DumpNLO",
    "DumpPOWHEG",
    "DiagnosePOWHEG",
    "POWHEGRealSpinDiag",
    "NLO_AUDIT",
    "NLO_TERM",
    "POWHEG_RAW",
    "LO_GAMMA_POINT",
    "POWHEG_SPIN",
    "NLOWeightRawComponents",
    "UseRawFinitePolarizedNLODeltas",
    "UseUniformPolarizedNLORepresentation",
    "FixedOrderValidation",
)

FORBIDDEN_CARD_TOKENS = (
    "DISDiagnostics",
    "DumpNLO",
    "DumpPOWHEG",
    "DiagnosePOWHEG",
    "POWHEGRealSpinDiag",
    "FixedOrderValidation",
)


def layout_candidates(root):
    return (
        (
            root / "DISPOL",
            root / "HerwigSource" / "Herwig-7.3.0",
        ),
        (
            root / "workflow" / "dispol",
            root / "src" / "herwig",
        ),
    )


def find_layout():
    for root in Path(__file__).resolve().parents:
        for workflow_root, herwig_root in layout_candidates(root):
            if workflow_root.exists() and herwig_root.exists():
                return workflow_root, herwig_root
    raise RuntimeError("Could not locate a supported DISPOL/Herwig layout")


WORKFLOW_ROOT, HERWIG_ROOT = find_layout()


def display_path(path):
    for root in (WORKFLOW_ROOT.parent, HERWIG_ROOT.parent, path.anchor):
        try:
            return str(path.relative_to(root))
        except ValueError:
            continue
    return str(path)


class ReleaseCleanupTests(unittest.TestCase):

    def test_release_sources_do_not_expose_diagnostic_interfaces_or_logs(self):
        source_roots = (
            HERWIG_ROOT / "MatrixElement" / "DIS",
            HERWIG_ROOT / "Shower" / "QTilde",
        )
        files = [
            path
            for source_root in source_roots
            for path in source_root.glob("*.[ch]*")
        ]
        self.assertGreater(len(files), 0)

        offenders = []
        for path in files:
            text = path.read_text(errors="replace")
            for token in FORBIDDEN_SOURCE_TOKENS:
                if token in text:
                    offenders.append(f"{display_path(path)}: {token}")

        self.assertEqual([], offenders)

    def test_active_cards_do_not_set_removed_diagnostic_interfaces(self):
        card_root = WORKFLOW_ROOT / "cards"
        files = list(card_root.glob("*.in")) + list(card_root.glob("*.template"))
        self.assertGreater(len(files), 0)

        offenders = []
        for path in files:
            text = path.read_text(errors="replace")
            for token in FORBIDDEN_CARD_TOKENS:
                if token in text:
                    offenders.append(f"{path.relative_to(WORKFLOW_ROOT)}: {token}")

        self.assertEqual([], offenders)

    def test_active_workflow_scripts_do_not_generate_removed_diagnostics(self):
        script_root = WORKFLOW_ROOT / "scripts"
        files = [
            path
            for path in script_root.glob("*.py")
            if path.name != Path(__file__).name
        ]
        self.assertGreater(len(files), 0)

        offenders = []
        for path in files:
            text = path.read_text(errors="replace")
            for token in FORBIDDEN_SOURCE_TOKENS:
                if token in text:
                    offenders.append(f"{path.relative_to(WORKFLOW_ROOT)}: {token}")

        self.assertEqual([], offenders)

    def test_fixed_order_no_shower_release_mode_is_retained(self):
        qtilde_cc = HERWIG_ROOT / "Shower" / "QTilde" / "QTildeShowerHandler.cc"
        qtilde_h = HERWIG_ROOT / "Shower" / "QTilde" / "QTildeShowerHandler.h"
        source = qtilde_cc.read_text(errors="replace")
        header = qtilde_h.read_text(errors="replace")

        self.assertIn("POWHEGEmissionMode", source)
        self.assertIn("ShowerReconstructed", source)
        self.assertIn("FixedOrderNoShower", source)
        self.assertIn("insertFixedOrderPOWHEGRealEmission", source)
        self.assertIn("insertFixedOrderPOWHEGRealEmission", header)


if __name__ == "__main__":
    unittest.main()
