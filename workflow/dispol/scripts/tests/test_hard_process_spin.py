"""Source contracts supplement (not replace) compiled and event-level tests."""
import unittest
from pathlib import Path

def source_root():
    for root in Path(__file__).resolve().parents:
        for relative in ("src/herwig", "HerwigSource/Herwig-7.3.0"):
            candidate = root / relative
            if (candidate / "Shower/ShowerHandler.cc").is_file():
                return candidate
    raise RuntimeError("Herwig source not found")

class HardProcessSpinSourceTests(unittest.TestCase):
    def read(self, relative):
        return (source_root() / relative).read_text()

    def test_only_shower_polarization_is_gated(self):
        source = self.read("Shower/ShowerHandler.cc")
        self.assertIn("if (polExtractor && hardProcessSpin())", source)
        reset = source.index("isPolarized_ = make_pair(false,false);")
        gate = source.index("if (polExtractor && hardProcessSpin())")
        self.assertLess(reset, gate)
        self.assertIn("longpdfs_  = make_pair(PDFPtr(),PDFPtr());", source[reset:gate])
        self.assertNotIn("polExtractor->isPolarized(false)", source)

    def test_default_and_versioned_persistence(self):
        source = self.read("Shower/QTilde/QTildeShowerHandler.cc")
        self.assertIn("_hardProcessSpin(true)", source)
        self.assertIn("if (version >= 1) is >> _hardProcessSpin;", source)
        self.assertIn('"HwShower.so", 1)', source)
        self.assertIn("<< _partnerfinder << _hardProcessSpin", source)
        self.assertIn("_hardProcessSpin = true;", source)

    def test_detach_copies_and_clean_retries(self):
        source = self.read("Shower/QTilde/Base/ShowerTree.cc")
        self.assertIn("copy[ix]->spinInfo(SpinPtr())", source)
        self.assertNotIn("original[ix]->spinInfo(SpinPtr())", source)
        clear = source.split("void ShowerTree::clear()")[1].split("void ShowerTree::resetShowerProducts")[0]
        self.assertIn("if (!_inheritHardSpin) orig->spinInfo(SpinPtr())", clear)
        self.assertIn("ShowerTree(hard,inheritHardSpin)", source)
        self.assertIn("ShowerTree(it->second)", source)  # decay defaults unchanged

    def test_supported_scope_and_independent_spin_switch(self):
        source = self.read("Shower/QTilde/QTildeShowerHandler.cc")
        body = source.split("void QTildeShowerHandler::checkFlags()")[1]
        for guard in ("orderInAlphaS() == 2", "orderInAlphaEW() == 0",
                      "!isMPIOn()", "!currentTree()->truncatedShower()",
                      "!_hardme->hasPOWHEGCorrection()", "id <= 5"):
            self.assertIn(guard, body)
        self.assertNotIn("spinOpt_ =", source)

    def test_unpolarized_terminal_matrix_retains_physical_gluons(self):
        source = self.read("Shower/QTilde/SplittingFunctions/SudakovFormFactor.cc")
        body = source.split("RhoDMatrix SudakovFormFactor::calculateHMatrix")[1]
        self.assertLess(body.index("H(1,1) = 0.0"), body.index("if(!isPol_ || !pdf_) return H;"))
        self.assertLess(body.index("if(!isPol_ || !pdf_) return H;"), body.index("longPDF_->xfl"))

if __name__ == "__main__":
    unittest.main()
