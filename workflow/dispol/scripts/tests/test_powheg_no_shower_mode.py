import ast
import unittest
from pathlib import Path


def find_layout():
    for root in Path(__file__).resolve().parents:
        campaign = root / "DISPOL" / "scripts" / "run_validation_campaign.py"
        qtilde_cc = root / "HerwigSource" / "Herwig-7.3.0" / "Shower" / "QTilde" / "QTildeShowerHandler.cc"
        qtilde_h = root / "HerwigSource" / "Herwig-7.3.0" / "Shower" / "QTilde" / "QTildeShowerHandler.h"
        if campaign.exists() and qtilde_cc.exists() and qtilde_h.exists():
            return campaign, qtilde_cc, qtilde_h

        campaign = root / "workflow" / "dispol" / "scripts" / "run_validation_campaign.py"
        qtilde_cc = root / "src" / "herwig" / "Shower" / "QTilde" / "QTildeShowerHandler.cc"
        qtilde_h = root / "src" / "herwig" / "Shower" / "QTilde" / "QTildeShowerHandler.h"
        if campaign.exists() and qtilde_cc.exists() and qtilde_h.exists():
            return campaign, qtilde_cc, qtilde_h

    raise RuntimeError("Could not locate a supported DISPOL/Herwig source layout")


CAMPAIGN, QTILDE_CC, QTILDE_H = find_layout()


def load_campaign_functions(*names):
    tree = ast.parse(CAMPAIGN.read_text())
    selected = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
    module = ast.Module(
        body=[
            ast.Import(names=[ast.alias(name="re")]),
            *selected,
        ],
        type_ignores=[],
    )
    namespace = {}
    exec(compile(ast.fix_missing_locations(module), str(CAMPAIGN), "exec"), namespace)
    return namespace


class FixedOrderPOWHEGNoShowerTests(unittest.TestCase):

    def test_fixed_order_powheg_card_patch_inserts_runtime_mode(self):
        namespace = load_campaign_functions(
            "set_or_insert_card_setting",
            "rewrite_fixed_order_powheg_card_text",
        )
        text = "\n".join(
            [
                "set /Herwig/Shower/ShowerHandler:LimitEmissions HardOnly",
                "set /Herwig/Shower/ShowerHandler:HardEmission POWHEG",
                "set EventGenerator:EventHandler:HadronizationHandler  NULL",
                "set /Herwig/EventHandlers/EventHandler:DecayHandler NULL",
                "",
            ]
        )

        rendered = namespace["rewrite_fixed_order_powheg_card_text"](text, True)

        self.assertIn("set /Herwig/Shower/ShowerHandler:POWHEGEmissionMode FixedOrderNoShower", rendered)
        self.assertEqual(rendered.count("POWHEGEmissionMode"), 1)
        self.assertIn("set /Herwig/Shower/ShowerHandler:HardEmission POWHEG", rendered)

    def test_qtilde_declares_fixed_order_powheg_switch_and_helper(self):
        header = QTILDE_H.read_text()
        source = QTILDE_CC.read_text()

        self.assertIn("POWHEGEmissionMode", source)
        self.assertIn("ShowerReconstructed", source)
        self.assertIn("FixedOrderNoShower", source)
        self.assertIn("POWHEGEmissionMode FixedOrderNoShower requires", source)
        self.assertIn("insertFixedOrderPOWHEGRealEmission", header)
        self.assertIn("insertFixedOrderPOWHEGRealEmission", source)
        self.assertIn("_powhegEmissionMode", header)

    def test_fixed_order_powheg_setup_exits_before_evolution_partner_reconstruction(self):
        source = QTILDE_CC.read_text()
        setup_start = source.index("vector<ShowerProgenitorPtr> QTildeShowerHandler::setupShower")
        setup_end = source.index("void QTildeShowerHandler::setEvolutionPartners", setup_start)
        setup_body = source[setup_start:setup_end]

        helper_index = setup_body.index("insertFixedOrderPOWHEGRealEmission")
        return_index = setup_body.index("return vector<ShowerProgenitorPtr>();", helper_index)
        partners_index = setup_body.index("setEvolutionPartners")

        self.assertLess(helper_index, return_index)
        self.assertLess(return_index, partners_index)


if __name__ == "__main__":
    unittest.main()
