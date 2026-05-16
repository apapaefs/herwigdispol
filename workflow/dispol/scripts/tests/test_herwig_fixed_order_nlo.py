import importlib.util
import io
import math
import sys
import tempfile
import unittest
from pathlib import Path


def scripts_dir() -> Path:
    for root in Path(__file__).resolve().parents:
        candidate = root / "DISPOL" / "scripts"
        if candidate.exists():
            return candidate
    raise RuntimeError("Could not locate DISPOL/scripts")


SCRIPT_DIR = scripts_dir()
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


class HerwigFixedOrderNLOTests(unittest.TestCase):

    def test_default_config_matches_validation_plan(self):
        from herwig_fo import config

        run = config.default_run_config()

        self.assertEqual((18.0, 275.0), (run.beams.lepton_energy, run.beams.hadron_energy))
        self.assertEqual((100.0, 2500.0), (run.cuts.q2_min, run.cuts.q2_max))
        self.assertEqual((0.2, 0.6), (run.cuts.y_min, run.cuts.y_max))
        self.assertEqual("NNPDF40_nlo_pch_as_01180", run.pdf_profile.unpolarized)
        self.assertEqual("NNPDFpol20_nlo_as_01180", run.pdf_profile.polarized_difference)
        self.assertEqual(("GAMMA", "Z", "ALL"), config.NC_SETUPS)
        self.assertEqual((1.0, -1.0), config.HELICITIES["PM"])

    def test_plain66_window_matches_reference_campaign_cuts(self):
        from herwig_fo import config

        cuts = config.cuts_for_window("plain66")

        self.assertEqual((49.0, 2500.0), (cuts.q2_min, cuts.q2_max))
        self.assertEqual((0.2, 0.6), (cuts.y_min, cuts.y_max))

    def test_matchbox_alpha_s_is_normalized_at_mz(self):
        from herwig_fo.alphas import MatchboxAlphaS

        running = MatchboxAlphaS(alpha_s_mz=0.118, mz=91.188)

        self.assertAlmostEqual(0.118, running.alpha_s(91.188**2), places=12)
        self.assertGreater(running.alpha_s(100.0), running.alpha_s(2500.0))

    def test_photon_limit_and_a_pol(self):
        from herwig_fo import ew

        ell = 2.0 / 0.4 - 1.0

        self.assertAlmostEqual(-2.0, ew.a_pol("GAMMA", 11, 2, 500.0, 1.0, -1.0))
        self.assertAlmostEqual(
            1.0 + ell * ell - 2.0 * ell,
            ew.sigma_born_factor("GAMMA", 11, 2, 500.0, 1.0, -1.0, ell),
        )
        coeff = ew.nc_coefficients("GAMMA", 11, 2, 500.0)
        self.assertAlmostEqual(4.0 / 9.0, coeff.D0)
        self.assertAlmostEqual(8.0 / 9.0, coeff.Nlq)

    def test_born_weight_uses_dis_dq2dy_structure_function_normalization(self):
        from herwig_fo import ew
        from herwig_fo.born import DISPoint, PB_PER_GEV2, born_weight
        from herwig_fo.pdfs import PDFPair, ToyPDFProvider

        point = DISPoint(q2=500.0, y=0.4, x_b=0.025, flavor=2, setup="GAMMA")
        pair = PDFPair(ToyPDFProvider("sum"), ToyPDFProvider("difference"))

        y_plus = 1.0 + (1.0 - point.y) ** 2
        xf = pair.unpolarized.xfx(point.flavor, point.x_b, point.q2)
        expected = (
            2.0
            * math.pi
            * ew.DEFAULT_EW_INPUTS.alpha_em
            * ew.DEFAULT_EW_INPUTS.alpha_em
            * PB_PER_GEV2
            / (point.y * point.q2 * point.q2)
            * xf
            * (4.0 / 9.0)
            * y_plus
        )

        self.assertAlmostEqual(expected, born_weight(point, pair), delta=abs(expected) * 1.0e-12)

    def test_pdf_ratio_helpers_keep_flux_unpolarized(self):
        from herwig_fo.pdfs import PDFPair, ToyPDFProvider, ratios_for_flavor

        pair = PDFPair(ToyPDFProvider("sum"), ToyPDFProvider("difference"))
        ratios = ratios_for_flavor(pair, flavor=2, x_b=0.02, x_p=0.7, q2=400.0, hadron_pol=1.0)

        self.assertGreater(ratios.lo_pdf, 0.0)
        self.assertGreater(ratios.q_ratio, 0.0)
        self.assertGreater(ratios.g_ratio, 0.0)
        self.assertLessEqual(abs(ratios.pq_born), 1.0)
        self.assertLessEqual(abs(ratios.pq_mapped), 1.0)
        self.assertLessEqual(abs(ratios.pg_mapped), 1.0)

    def test_nlo_terms_decompose_into_herwig_pieces(self):
        from herwig_fo.alphas import MatchboxAlphaS
        from herwig_fo.born import DISPoint
        from herwig_fo.nlo import nlo_terms
        from herwig_fo.pdfs import PDFPair, ToyPDFProvider

        point = DISPoint(q2=500.0, y=0.35, x_b=0.025, flavor=2, setup="ALL", lepton_pol=1.0, hadron_pol=1.0)
        pair = PDFPair(ToyPDFProvider("sum"), ToyPDFProvider("difference"))
        terms = nlo_terms(point, x_p=0.65, pdfs=pair, alpha_s=MatchboxAlphaS().alpha_s)

        pieces = terms.virtual + terms.collinear_quark + terms.collinear_gluon + terms.real_qcdc + terms.real_bgf
        self.assertTrue(math.isfinite(terms.total))
        self.assertAlmostEqual(pieces, terms.total)
        self.assertNotEqual(0.0, terms.virtual)
        self.assertNotEqual(0.0, terms.real_qcdc)

    def test_photon_local_real_kernels_match_breit_spin_tensors(self):
        from herwig_fo.born import DISPoint
        from herwig_fo.ew import sigma_born_factor
        from herwig_fo.pdfs import PDFRatios
        from herwig_fo.real import (
            _real_variables,
            bgf_azimuthal_coefficients,
            compton_azimuthal_coefficients,
        )

        def ratios(pq: float, pg: float) -> PDFRatios:
            return PDFRatios(
                lo_pdf=1.0,
                q_pdf=1.0,
                g_pdf=1.0,
                dq_pdf=0.0,
                dg_pdf=0.0,
                dlo_pdf=0.0,
                q_ratio=1.0,
                g_ratio=1.0,
                deltaq_over_lo=0.0,
                deltag_over_lo=0.0,
                deltaq_over_lo_eff=0.0,
                deltag_over_lo_eff=0.0,
                pq_born=pq,
                pq_mapped=pq,
                pg_mapped=pg,
                has_difference_pdf=True,
                has_stable_difference_ratio=True,
            )

        def dot(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> float:
            return a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]

        def breit_tensor_ratios(q2: float, y: float, x_p: float, z_p: float, phi: float) -> tuple[float, float]:
            q = math.sqrt(q2)
            x_perp, x1, x2, x3 = _real_variables(x_p, z_p)
            p1 = (0.0, 0.0, -0.5 * x1 * q, -0.5 * x1 * q)
            p2 = (
                0.5 * q * x_perp * math.cos(phi),
                0.5 * q * x_perp * math.sin(phi),
                -0.5 * q * x2,
                0.5 * q * math.hypot(x_perp, x2),
            )
            p3 = (-p2[0], -p2[1], -0.5 * q * x3, 0.5 * q * math.hypot(x_perp, x3))
            ell = 2.0 / y - 1.0
            root = math.sqrt(ell * ell - 1.0)
            lepton_in = (0.5 * q * root, 0.0, -0.5 * q, 0.5 * q * ell)
            lepton_out = (0.5 * q * root, 0.0, 0.5 * q, 0.5 * q * ell)

            qcdc_unpol = (
                dot(p1, lepton_in) ** 2
                + dot(p1, lepton_out) ** 2
                + dot(p2, lepton_out) ** 2
                + dot(p2, lepton_in) ** 2
            )
            qcdc_pol = (
                dot(p1, lepton_in) ** 2
                - dot(p1, lepton_out) ** 2
                + dot(p2, lepton_out) ** 2
                - dot(p2, lepton_in) ** 2
            )
            bgf_unpol = (
                dot(p3, lepton_in) ** 2
                + dot(p3, lepton_out) ** 2
                + dot(p2, lepton_out) ** 2
                + dot(p2, lepton_in) ** 2
            )
            bgf_pol = (
                -dot(p3, lepton_in) ** 2
                + dot(p3, lepton_out) ** 2
                + dot(p2, lepton_out) ** 2
                - dot(p2, lepton_in) ** 2
            )
            return qcdc_pol / qcdc_unpol, bgf_pol / bgf_unpol

        def herwig_ratio(channel: str, q2: float, y: float, x_p: float, z_p: float, phi: float) -> float:
            ell = 2.0 / y - 1.0

            def born_factor(lepton_pol: float, quark_pol: float) -> float:
                return sigma_born_factor("GAMMA", 11, 2, q2, lepton_pol, quark_pol, ell)

            def weighted_kernel(lepton_pol: float, parton_pol: float) -> float:
                point = DISPoint(q2=q2, y=y, x_b=0.01, flavor=2, setup="GAMMA", lepton_pol=lepton_pol)
                if channel == "qcdc":
                    kernel = compton_azimuthal_coefficients(point, x_p, z_p, ratios(parton_pol, 0.0)).value(phi)
                    return born_factor(lepton_pol, parton_pol) * kernel
                kernel = bgf_azimuthal_coefficients(point, x_p, z_p, ratios(0.0, parton_pol)).value(phi)
                return born_factor(lepton_pol, 0.0) * kernel

            pp = weighted_kernel(1.0, 1.0)
            pm = weighted_kernel(1.0, -1.0)
            mp = weighted_kernel(-1.0, 1.0)
            mm = weighted_kernel(-1.0, -1.0)
            sigma = 0.25 * math.fsum((pp, pm, mp, mm))
            delta = 0.25 * math.fsum((pp, mm, -pm, -mp))
            return delta / sigma

        for q2, y, x_p, z_p, phi in (
            (200.0, 0.4, 0.6, 0.3, 1.1),
            (1000.0, 0.25, 0.8, 0.7, 2.2),
            (500.0, 0.55, 0.4, 0.2, 0.0),
        ):
            qcdc_tensor, bgf_tensor = breit_tensor_ratios(q2, y, x_p, z_p, phi)
            self.assertAlmostEqual(qcdc_tensor, herwig_ratio("qcdc", q2, y, x_p, z_p, phi), places=12)
            self.assertAlmostEqual(bgf_tensor, herwig_ratio("bgf", q2, y, x_p, z_p, phi), places=12)

    def test_qcdc_leading_term_uses_mapped_polarized_luminosity(self):
        from herwig_fo.born import DISPoint
        from herwig_fo.ew import sigma_born_factor
        from herwig_fo.pdfs import PDFRatios
        from herwig_fo.real import compton_azimuthal_coefficients

        q2 = 500.0
        y = 0.35
        x_p = 0.62
        z_p = 0.41
        phi = 1.3
        ell = 2.0 / y - 1.0
        lepton_pol = 1.0
        pq_born = 0.15
        pq_mapped = 0.55

        def ratios(born_pol: float, mapped_pol: float) -> PDFRatios:
            return PDFRatios(
                lo_pdf=1.0,
                q_pdf=1.0,
                g_pdf=1.0,
                dq_pdf=0.0,
                dg_pdf=0.0,
                dlo_pdf=0.0,
                q_ratio=1.0,
                g_ratio=1.0,
                deltaq_over_lo=0.0,
                deltag_over_lo=0.0,
                deltaq_over_lo_eff=0.0,
                deltag_over_lo_eff=0.0,
                pq_born=born_pol,
                pq_mapped=mapped_pol,
                pg_mapped=0.0,
                has_difference_pdf=True,
                has_stable_difference_ratio=True,
            )

        point = DISPoint(q2=q2, y=y, x_b=0.01, flavor=2, setup="GAMMA", lepton_pol=lepton_pol)
        mixed = compton_azimuthal_coefficients(point, x_p, z_p, ratios(pq_born, pq_mapped)).value(phi)
        mapped = compton_azimuthal_coefficients(point, x_p, z_p, ratios(pq_mapped, pq_mapped)).value(phi)
        born_norm = sigma_born_factor("GAMMA", 11, 2, q2, lepton_pol, pq_born, ell)
        mapped_norm = sigma_born_factor("GAMMA", 11, 2, q2, lepton_pol, pq_mapped, ell)

        self.assertAlmostEqual(born_norm * mixed, mapped_norm * mapped, places=12)

    def test_herwig_xp_power_map_matches_disbase_jacobian(self):
        from herwig_fo.nlo import herwig_xp_power_map

        x_b = 0.025
        power = 0.1
        random_unit = 0.37
        x_p, jacobian = herwig_xp_power_map(random_unit, x_b, power)
        rhomin = (1.0 - x_b) ** (1.0 - power)
        expected_xp = 1.0 - (random_unit * rhomin) ** (1.0 / (1.0 - power))
        expected_jacobian = rhomin / (1.0 - power) * (1.0 - expected_xp) ** power

        self.assertAlmostEqual(expected_xp, x_p)
        self.assertAlmostEqual(expected_jacobian, jacobian)
        self.assertGreater(x_p, x_b)
        self.assertLess(x_p, 1.0)

    def test_generate_events_uses_full_nlo_seed_and_deterministic_flavour_sum(self):
        from herwig_fo.alphas import MatchboxAlphaS
        from herwig_fo.born import born_weight
        from herwig_fo.config import ACTIVE_QUARK_FLAVORS, default_run_config
        from herwig_fo.nlo import herwig_xp_power_map, nlo_terms
        from herwig_fo.pdfs import load_pdf_pair
        from herwig_fo.real import local_real_weight_ratios
        import herwig_fixed_order_nlo as cli

        args = cli.parse_args(["generate", "--events", "1", "--output", "unused.hepmc3", "--toy-pdfs", "--sampler", "random"])
        events, sumw = cli.generate_events(args)

        run = cli.run_config_from_args(args)
        rng = __import__("random").Random(args.seed)
        q2, y, x_b, q2_jacobian = cli._sample_coordinates(rng, run, args.q2_sampling)
        x_p, x_p_jacobian = herwig_xp_power_map(rng.random(), x_b)
        z_p = min(max(rng.random(), 1.0e-12), 1.0 - 1.0e-12)
        phi = 2.0 * math.pi * rng.random()
        q2_volume, _ = cli._q2_volume_and_jacobian(run.cuts.q2_min, run.cuts, args.q2_sampling)
        phase_space_volume = q2_volume * (run.cuts.y_max - run.cuts.y_min)
        pdfs = load_pdf_pair(run.pdf_profile, use_toy=True)
        alpha_s = MatchboxAlphaS(run.alpha_s_mz, run.mz).alpha_s
        expected_seed = 0.0
        for flavor in ACTIVE_QUARK_FLAVORS:
            point = cli._point_from_coordinates(q2, y, x_b, flavor, args.setup, args.helicity, run)
            terms = nlo_terms(point, x_p, pdfs, alpha_s, x_p_jacobian)
            expected_seed += phase_space_volume * q2_jacobian * born_weight(point, pdfs) * x_p_jacobian * terms.total

        first_point = cli._point_from_coordinates(q2, y, x_b, ACTIVE_QUARK_FLAVORS[0], args.setup, args.helicity, run)
        first_base = (
            phase_space_volume
            * q2_jacobian
            * born_weight(first_point, pdfs)
            * x_p_jacobian
            / args.events
        )
        first_local_real = local_real_weight_ratios(first_point, x_p, z_p, phi, pdfs, alpha_s)

        self.assertEqual(5 * len(ACTIVE_QUARK_FLAVORS), len(events))
        self.assertAlmostEqual(first_base * first_local_real.qcdc, events[1].weight)
        self.assertAlmostEqual(-(first_base * first_local_real.qcdc), events[2].weight)
        self.assertAlmostEqual(first_base * first_local_real.bgf, events[3].weight)
        self.assertAlmostEqual(-(first_base * first_local_real.bgf), events[4].weight)
        self.assertAlmostEqual(expected_seed, sum(event.weight for event in events[0::5]), delta=max(1.0, abs(expected_seed)) * 1.0e-12)
        self.assertAlmostEqual(expected_seed, sumw, delta=max(1.0, abs(sumw)) * 1.0e-12)

    def test_qcdc_breit_momenta_are_transversely_balanced(self):
        from herwig_fo.real import breit_qcdc_momenta

        momenta = breit_qcdc_momenta(q2=400.0, x_p=0.55, z_p=0.3, phi=0.7)
        p1 = momenta.outgoing[0]
        p2 = momenta.outgoing[1]

        self.assertAlmostEqual(0.0, p1.px + p2.px, places=12)
        self.assertAlmostEqual(0.0, p1.py + p2.py, places=12)
        self.assertGreater(p1.pt(), 0.0)
        self.assertGreater(p2.pt(), 0.0)
        self.assertEqual((2, 21), tuple(p.pid for p in momenta.outgoing))

    def test_bgf_breit_momenta_preserve_born_flavor_orientation(self):
        from herwig_fo.real import breit_bgf_momenta

        momenta = breit_bgf_momenta(q2=400.0, x_p=0.55, z_p=0.3, phi=0.7, flavor=-2)

        self.assertEqual((-2, 2), tuple(p.pid for p in momenta.outgoing))

    def test_local_real_weight_ratios_use_full_herwig_kernels(self):
        from herwig_fo.alphas import MatchboxAlphaS
        from herwig_fo.born import DISPoint
        from herwig_fo.pdfs import PDFPair, ToyPDFProvider
        from herwig_fo.real import local_real_weight_ratios

        point = DISPoint(q2=500.0, y=0.35, x_b=0.025, flavor=2, setup="GAMMA")
        pair = PDFPair(ToyPDFProvider("sum"), ToyPDFProvider("difference"))
        ratios = local_real_weight_ratios(point, x_p=0.65, z_p=0.4, phi=0.7, pdfs=pair, alpha_s=MatchboxAlphaS().alpha_s)

        self.assertTrue(math.isfinite(ratios.qcdc))
        self.assertTrue(math.isfinite(ratios.bgf))
        self.assertGreater(ratios.qcdc, 0.0)
        self.assertGreater(ratios.bgf, 0.0)

    def test_bgf_real_weight_uses_quark_antiquark_endpoint_partition(self):
        from herwig_fo.alphas import MatchboxAlphaS
        from herwig_fo.born import DISPoint
        from herwig_fo.pdfs import PDFPair
        from herwig_fo.real import bgf_real_weight_ratio

        class SymmetricPDF:
            name = "symmetric"

            def xfx(self, pid, x, q2):
                if x <= 0.0 or x >= 1.0:
                    return 0.0
                return x * (5.0 if pid == 21 else 2.0)

            def alpha_s(self, q2):
                return None

        pair = PDFPair(SymmetricPDF())
        alpha_s = MatchboxAlphaS().alpha_s
        kwargs = dict(x_p=0.65, z_p=0.4, phi=0.7, pdfs=pair, alpha_s=alpha_s)
        quark = DISPoint(q2=500.0, y=0.35, x_b=0.025, flavor=2, setup="GAMMA")
        antiquark = DISPoint(q2=500.0, y=0.35, x_b=0.025, flavor=-2, setup="GAMMA")

        q_weight = bgf_real_weight_ratio(quark, **kwargs)
        qbar_weight = bgf_real_weight_ratio(antiquark, **kwargs)

        self.assertAlmostEqual((1.0 - kwargs["z_p"]) / kwargs["z_p"], qbar_weight / q_weight)

    def test_real_event_uses_mapped_parton_remnant_and_conserves_momentum(self):
        from herwig_fo.born import DISPoint
        from herwig_fo.real import real_event

        point = DISPoint(q2=400.0, y=0.4, x_b=400.0 / (4.0 * 18.0 * 275.0 * 0.4), flavor=2, setup="GAMMA")
        event = real_event(point, "QCDC", x_p=0.55, z_p=0.3, phi=0.7, weight=1.0, event_number=1)

        incoming = event.incoming[0].momentum + event.incoming[1].momentum
        outgoing = event.outgoing[0].momentum
        for particle in event.outgoing[1:]:
            outgoing = outgoing + particle.momentum

        self.assertAlmostEqual(incoming.px, outgoing.px, delta=1.0e-9)
        self.assertAlmostEqual(incoming.py, outgoing.py, delta=1.0e-9)
        self.assertAlmostEqual(incoming.pz, outgoing.pz, delta=1.0e-9)
        self.assertAlmostEqual(incoming.e, outgoing.e, delta=1.0e-9)
        self.assertLess(event.outgoing[-1].momentum.pz, 275.0 * (1.0 - point.x_b))

    def test_hepmc_writer_emits_signed_weights_and_me_tags(self):
        from herwig_fo.born import DISPoint, born_event
        from herwig_fo.hepmc import write_hepmc3

        point = DISPoint(q2=300.0, y=0.4, x_b=0.0151515151515, flavor=2, setup="GAMMA")
        event = born_event(point, weight=-2.5, event_number=7)

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "events.hepmc3"
            write_hepmc3(path, [event])
            text = path.read_text()

        self.assertIn("HepMC::Asciiv3-START_EVENT_LISTING", text)
        self.assertIn("W Weight", text)
        self.assertIn("GenCrossSection -2.500000000000e+00", text)
        self.assertIn(" 990072", text)
        self.assertIn("herwig_powheg_me_parton 2", text)

    def test_hepmc_writer_emits_named_helicity_weights(self):
        from herwig_fo.born import DISPoint, born_event
        from herwig_fo.hepmc import write_hepmc3

        point = DISPoint(q2=300.0, y=0.4, x_b=0.0151515151515, flavor=2, setup="GAMMA")
        event = born_event(point, weight=2.0, event_number=7).with_weights(
            {
                "HERWIG_DIS_PP": 1.1,
                "HERWIG_DIS_PM": 0.9,
                "HERWIG_DIS_SIGMA": 1.0,
                "HERWIG_DIS_DELTA_LL": 0.1,
            }
        )

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "events.hepmc3"
            write_hepmc3(path, [event])
            text = path.read_text()

        self.assertIn("W Weight HERWIG_DIS_PP HERWIG_DIS_PM HERWIG_DIS_SIGMA HERWIG_DIS_DELTA_LL", text)
        self.assertIn(
            "W 2.000000000000e+00 1.100000000000e+00 9.000000000000e-01 "
            "1.000000000000e+00 1.000000000000e-01",
            text,
        )

    def test_scale_factor_sets_dispoint_mu2(self):
        import herwig_fixed_order_nlo as cli
        from herwig_fo.config import default_run_config

        point = cli._point_from_coordinates(
            q2=500.0,
            y=0.35,
            x_b=0.025,
            flavor=2,
            setup="ALL",
            helicity="PP",
            run=default_run_config(),
            scale_factor=2.0,
        )

        self.assertEqual(2000.0, point.mu2)
        self.assertEqual(2000.0, point.scale2)

    def test_correlated_weight_bundle_matches_independent_helicity_weights(self):
        import herwig_fixed_order_nlo as cli
        from herwig_fo.alphas import MatchboxAlphaS
        from herwig_fo.config import default_run_config
        from herwig_fo.pdfs import load_pdf_pair

        run = default_run_config()
        pdfs = load_pdf_pair(run.pdf_profile, use_toy=True)
        alpha_s = MatchboxAlphaS(run.alpha_s_mz, run.mz).alpha_s
        q2 = 500.0
        y = 0.35
        x_b = cli.xb_from_q2_y(q2, y, run.beams)
        x_p = 0.65
        z_p = 0.4
        phi = 0.7
        base = 3.0

        bundle = cli.correlated_helicity_weights(
            q2=q2,
            y=y,
            x_b=x_b,
            flavor=2,
            setup="ALL",
            contribution="inclusive",
            base=base,
            x_p=x_p,
            z_p=z_p,
            phi=phi,
            pdfs=pdfs,
            alpha_s=alpha_s,
            run=run,
        )

        pp = cli.single_helicity_contribution_weight("PP", q2, y, x_b, 2, "ALL", "inclusive", base, x_p, z_p, phi, pdfs, alpha_s, run=run)
        pm = cli.single_helicity_contribution_weight("PM", q2, y, x_b, 2, "ALL", "inclusive", base, x_p, z_p, phi, pdfs, alpha_s, run=run)
        mp = cli.single_helicity_contribution_weight("MP", q2, y, x_b, 2, "ALL", "inclusive", base, x_p, z_p, phi, pdfs, alpha_s, run=run)
        mm = cli.single_helicity_contribution_weight("MM", q2, y, x_b, 2, "ALL", "inclusive", base, x_p, z_p, phi, pdfs, alpha_s, run=run)

        self.assertAlmostEqual(0.25 * (pp + pm + mp + mm), bundle.sigma)
        self.assertAlmostEqual(0.25 * (pp + mm - pm - mp), bundle.delta_ll)
        self.assertAlmostEqual(bundle.sigma, bundle.nominal)
        self.assertAlmostEqual(bundle.pp, pp)
        self.assertAlmostEqual(bundle.pm, pm)
        self.assertAlmostEqual(bundle.mp, mp)
        self.assertAlmostEqual(bundle.mm, mm)
        self.assertAlmostEqual(1.0, bundle.named_multipliers()["HERWIG_DIS_SIGMA"])

    def test_projected_nlo_correlated_real_weights_keep_carrier_and_use_inclusive_ratio(self):
        import herwig_fixed_order_nlo as cli
        from herwig_fo.alphas import MatchboxAlphaS
        from herwig_fo.config import default_run_config
        from herwig_fo.pdfs import load_pdf_pair

        run = default_run_config()
        pdfs = load_pdf_pair(run.pdf_profile, use_toy=True)
        alpha_s = MatchboxAlphaS(run.alpha_s_mz, run.mz).alpha_s
        q2 = 500.0
        y = 0.35
        x_b = cli.xb_from_q2_y(q2, y, run.beams)
        x_p = 0.65
        z_p = 0.4
        phi = 0.7
        base = 3.0

        local = cli.correlated_helicity_weights(q2, y, x_b, 2, "ALL", "bgf", base, x_p, z_p, phi, pdfs, alpha_s, run=run)
        inclusive = cli.correlated_helicity_weights(
            q2, y, x_b, 2, "ALL", "inclusive", base, x_p, z_p, phi, pdfs, alpha_s, run=run
        )
        projected = cli.projected_nlo_correlated_helicity_weights(
            q2, y, x_b, 2, "ALL", "bgf", base, x_p, z_p, phi, pdfs, alpha_s, run=run
        )

        self.assertAlmostEqual(local.nominal, projected.nominal)
        self.assertAlmostEqual(local.nominal, projected.sigma)
        self.assertAlmostEqual(
            inclusive.delta_ll / inclusive.sigma,
            projected.delta_ll / projected.sigma,
            delta=1.0e-12,
        )
        self.assertAlmostEqual(
            inclusive.pp / inclusive.sigma,
            projected.pp / projected.sigma,
            delta=1.0e-12,
        )

    def test_correlated_generation_preserves_unpolarized_sumw_with_xp_jacobian(self):
        import herwig_fixed_order_nlo as cli

        generate_args = cli.parse_args(
            [
                "generate",
                "--setup",
                "GAMMA",
                "--events",
                "3",
                "--output",
                "unused.hepmc3",
                "--toy-pdfs",
            ]
        )
        correlated_args = cli.parse_args(
            [
                "polarized-full",
                "--setup",
                "GAMMA",
                "--events",
                "3",
                "--work-dir",
                "work",
                "--output",
                "polarized.yoda",
                "--toy-pdfs",
            ]
        )

        _, generate_sumw = cli.generate_events(generate_args)
        _, correlated_sumw = cli.generate_correlated_events(correlated_args)

        self.assertAlmostEqual(generate_sumw, correlated_sumw, delta=max(1.0, abs(generate_sumw)) * 1.0e-12)

    def test_run_rivet_builds_rivetweights_analysis_string(self):
        import herwig_fo.rivet as rivet

        calls = []
        original_run = rivet.subprocess.run

        def fake_run(command, check, env):
            calls.append((command, check, env))

        try:
            rivet.subprocess.run = fake_run
            rivet.run_rivet(
                Path("events.hepmc3"),
                Path("out.yoda"),
                jetinput="MEPARTONS",
                rivetweights=True,
                analysis_plugin=Path("/tmp/RivetMC_DIS_BREIT.so"),
            )
        finally:
            rivet.subprocess.run = original_run

        command = calls[0][0]
        self.assertIn("MC_DIS_BREIT:JETINPUT=MEPARTONS:RIVETWEIGHTS=YES", command)

    def test_cli_argument_defaults_are_fixed_order(self):
        script = SCRIPT_DIR / "herwig_fixed_order_nlo.py"
        spec = importlib.util.spec_from_file_location("herwig_fixed_order_nlo_for_test", script)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        old_module = sys.modules.get("herwig_fixed_order_nlo_for_test")
        sys.modules["herwig_fixed_order_nlo_for_test"] = module
        try:
            spec.loader.exec_module(module)
        finally:
            if old_module is None:
                sys.modules.pop("herwig_fixed_order_nlo_for_test", None)
            else:
                sys.modules["herwig_fixed_order_nlo_for_test"] = old_module

        args = module.parse_args(["generate", "--events", "10", "--output", "out.hepmc3", "--toy-pdfs"])

        self.assertEqual("generate", args.command)
        self.assertEqual("ALL", args.setup)
        self.assertEqual("00", args.helicity)
        self.assertEqual("validation", args.window)
        self.assertEqual("log", args.q2_sampling)
        self.assertEqual("halton", args.sampler)
        self.assertFalse(args.powheg)
        self.assertEqual("MEPARTONS", args.jetinput)
        self.assertEqual(1.0, args.scale_factor)

    def test_polarize_command_defaults_to_breit_helicity_combination(self):
        import herwig_fixed_order_nlo as cli

        args = cli.parse_args(
            [
                "polarize",
                "--PP",
                "pp.yoda",
                "--PM",
                "pm.yoda",
                "--MP",
                "mp.yoda",
                "--MM",
                "mm.yoda",
                "--output",
                "polarized.yoda",
            ]
        )

        self.assertEqual("polarize", args.command)
        self.assertEqual("ALL", args.setup)
        self.assertEqual("MC_DIS_BREIT", args.analysis)
        self.assertEqual("polarized.yoda", str(args.output))

    def test_polarized_full_uses_required_helicity_set(self):
        import herwig_fixed_order_nlo as cli

        self.assertEqual(("00", "PP", "PM", "MP", "MM"), cli.required_polarized_helicities("ALL"))
        self.assertEqual(("00", "PP", "PM", "MP", "MM"), cli.required_polarized_helicities("Z"))
        self.assertEqual(("00", "PP", "PM"), cli.required_polarized_helicities("GAMMA"))

    def test_polarized_full_defaults_to_fixed_order_correlated_weights_and_exposes_scale_variations(self):
        import herwig_fixed_order_nlo as cli

        args = cli.parse_args(
            [
                "polarized-full",
                "--events",
                "10",
                "--work-dir",
                "work",
                "--output",
                "polarized.yoda",
                "--scale-variations",
            ]
        )

        self.assertFalse(args.independent_helicities)
        self.assertTrue(args.scale_variations)
        self.assertEqual(1.0, args.scale_factor)
        self.assertEqual("local-real", args.polarized_real_weight_mode)

    def test_standalone_campaign_defaults_match_fixed_order_plot_workflow(self):
        import herwig_fixed_order_nlo as cli

        args = cli.parse_args(["standalone-campaign", "-t", "standalone_fo_q2gt100_rivet05_localreal"])

        self.assertEqual("standalone-campaign", args.command)
        self.assertEqual(20000, args.events)
        self.assertTrue(args.scale_variations)
        self.assertEqual("local-real", args.polarized_real_weight_mode)
        self.assertIn("rivetweights_noshowerMac03", str(args.poldis_reference))
        self.assertIn("DQ2PreCut", args.match)
        self.assertIn("ALLZeta", args.match)

    def test_standalone_campaign_exposes_parallel_shard_controls(self):
        import herwig_fixed_order_nlo as cli

        args = cli.parse_args(
            [
                "standalone-campaign",
                "-t",
                "standalone_fo_q2gt100_highstat",
                "--events",
                "101",
                "--jobs",
                "4",
                "--shards",
                "10",
                "--keep-shard-hepmc",
                "--yoda-merge-tool",
                "yodamerge",
                "--shard-monitor-interval",
                "2.5",
                "--shard-status-file",
                "status.tsv",
            ]
        )

        self.assertEqual(4, args.jobs)
        self.assertEqual(10, args.shards)
        self.assertTrue(args.keep_shard_hepmc)
        self.assertEqual("yodamerge", args.yoda_merge_tool)
        self.assertEqual(2.5, args.shard_monitor_interval)
        self.assertEqual(Path("status.tsv"), args.shard_status_file)

        specs = cli._standalone_shard_specs(total_events=args.events, shard_count=args.shards, seed=args.seed)
        self.assertEqual(args.shards, len(specs))
        self.assertEqual(args.events, sum(spec.events for spec in specs))
        self.assertEqual([11], [spec.events for spec in specs[:1]])
        self.assertEqual(10, specs[-1].events)
        self.assertEqual(args.seed, specs[0].seed)
        self.assertEqual(args.seed + specs[1].start, specs[1].seed)

    def test_shard_status_file_records_monitor_state(self):
        import herwig_fixed_order_nlo as cli

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            spec = cli.StandaloneShardSpec(index=3, start=30, events=10, seed=12375)
            task = cli.StandaloneShardTask(
                variation_name="ScaleFactorUp",
                scale_factor=2.0,
                spec=spec,
                work_dir=root / "shard_0003",
                hepmc_path=root / "shard_0003" / "events.hepmc3",
                raw_yoda_path=root / "shard_0003" / "raw.yoda",
                stdout_path=root / "shard_0003" / "stdout.log",
                stderr_path=root / "shard_0003" / "stderr.log",
                command=("python3", "script.py"),
            )
            status_path = root / "status.tsv"

            cli._write_shard_status_file(
                status_path,
                [task],
                {cli._standalone_shard_key(task): "running"},
                {},
                {cli._standalone_shard_key(task): 10.0},
                {},
                now=13.25,
            )

            lines = status_path.read_text().splitlines()
            self.assertEqual(
                "variation\tshard\tevents\tseed\tstatus\treturncode\telapsed_s\traw_yoda\tstdout\tstderr",
                lines[0],
            )
            fields = lines[1].split("\t")
            self.assertEqual("ScaleFactorUp", fields[0])
            self.assertEqual("3", fields[1])
            self.assertEqual("10", fields[2])
            self.assertEqual("12375", fields[3])
            self.assertEqual("running", fields[4])
            self.assertEqual("", fields[5])
            self.assertEqual("3.250", fields[6])
            self.assertEqual(str(task.raw_yoda_path), fields[7])
            self.assertEqual(str(task.stdout_path), fields[8])
            self.assertEqual(str(task.stderr_path), fields[9])

    def test_standalone_campaign_invokes_generation_plot_prep_and_mkhtml(self):
        import herwig_fixed_order_nlo as cli

        with tempfile.TemporaryDirectory() as tmp:
            base_dir = Path(tmp) / "DISPOL"
            rivet_dir = base_dir / "analyses" / "rivet" / "dis"
            rivet_dir.mkdir(parents=True)
            (rivet_dir / "MC_DIS_BREIT.plot").write_text("# plot config\n")
            reference = Path(tmp) / "poldis.yoda.gz"
            reference.write_text("# fake reference\n")

            args = cli.parse_args(
                [
                    "standalone-campaign",
                    "-t",
                    "smoke",
                    "--base-dir",
                    str(base_dir),
                    "--events",
                    "7",
                    "--poldis-reference",
                    str(reference),
                    "--rivet-mkhtml",
                    "/usr/bin/rivet-mkhtml",
                    "--no-progress",
                ]
            )

            calls = []
            original_polarized_full = cli.command_polarized_full
            original_plot_prep = cli.command_plot_prep
            original_run = cli.subprocess.run
            original_rewrite_scale = cli.rewrite_scale_envelope_plot_scripts
            original_rewrite_no_scale = cli.rewrite_no_scale_ratio_plot_scripts

            def fake_polarized_full(polarized_args):
                calls.append(("polarized-full", polarized_args))
                return 0

            def fake_plot_prep(plot_args):
                calls.append(("plot-prep", plot_args))
                return 0

            def fake_run(command, check, env):
                calls.append(("mkhtml", command, check, env))

            def fake_rewrite_scale(plot_dir, python_executable=None, reference_label="POLDIS NLO"):
                calls.append(("postprocess-scale", plot_dir, python_executable, reference_label))
                return 3, 2

            def fake_rewrite_no_scale(plot_dir, python_executable=None):
                calls.append(("postprocess-noscale", plot_dir, python_executable))
                return 0, 0

            try:
                cli.command_polarized_full = fake_polarized_full
                cli.command_plot_prep = fake_plot_prep
                cli.subprocess.run = fake_run
                cli.rewrite_scale_envelope_plot_scripts = fake_rewrite_scale
                cli.rewrite_no_scale_ratio_plot_scripts = fake_rewrite_no_scale
                from contextlib import redirect_stdout

                with redirect_stdout(io.StringIO()):
                    rc = cli.command_standalone_campaign(args)
            finally:
                cli.command_polarized_full = original_polarized_full
                cli.command_plot_prep = original_plot_prep
                cli.subprocess.run = original_run
                cli.rewrite_scale_envelope_plot_scripts = original_rewrite_scale
                cli.rewrite_no_scale_ratio_plot_scripts = original_rewrite_no_scale

            self.assertEqual(0, rc)
            self.assertEqual(
                ["polarized-full", "plot-prep", "mkhtml", "postprocess-scale"],
                [call[0] for call in calls],
            )
            polarized_args = calls[0][1]
            self.assertEqual(7, polarized_args.events)
            self.assertTrue(polarized_args.scale_variations)
            self.assertEqual("local-real", polarized_args.polarized_real_weight_mode)
            self.assertEqual(base_dir / "campaigns" / "smoke" / "work", polarized_args.work_dir)
            self.assertEqual(
                base_dir / "campaigns" / "smoke" / "analysis" / "standalone_fo_all_correlated_q2gt100.yoda",
                polarized_args.output,
            )
            plot_args = calls[1][1]
            self.assertTrue(plot_args.herwig_scale_envelope)
            self.assertTrue(plot_args.poldis.samefile(reference))
            mkhtml_command = calls[2][1]
            mkhtml_env = calls[2][3]
            self.assertIn("--verbose", mkhtml_command)
            self.assertIn("Standalone Python HERWIG FO", " ".join(mkhtml_command))
            self.assertIn(str(rivet_dir), mkhtml_env["RIVET_PLOT_PATH"])
            postprocess_call = calls[3]
            self.assertEqual(base_dir / "campaigns" / "smoke" / "plots_mc_vs_poldis_all_nlo_with_polarized", postprocess_call[1])
            self.assertEqual("POLDIS NLO", postprocess_call[3])

    def test_standalone_campaign_parallel_path_runs_shards_before_plotting(self):
        import herwig_fixed_order_nlo as cli

        with tempfile.TemporaryDirectory() as tmp:
            base_dir = Path(tmp) / "DISPOL"
            rivet_dir = base_dir / "analyses" / "rivet" / "dis"
            rivet_dir.mkdir(parents=True)
            (rivet_dir / "MC_DIS_BREIT.plot").write_text("# plot config\n")
            reference = Path(tmp) / "poldis.yoda.gz"
            reference.write_text("# fake reference\n")

            args = cli.parse_args(
                [
                    "standalone-campaign",
                    "-t",
                    "parallel-smoke",
                    "--base-dir",
                    str(base_dir),
                    "--events",
                    "13",
                    "--jobs",
                    "2",
                    "--shards",
                    "3",
                    "--poldis-reference",
                    str(reference),
                    "--rivet-mkhtml",
                    "/usr/bin/rivet-mkhtml",
                    "--no-progress",
                ]
            )

            calls = []
            original_parallel = cli._run_parallel_standalone_generation
            original_polarized_full = cli.command_polarized_full
            original_plot_prep = cli.command_plot_prep
            original_run = cli.subprocess.run
            original_rewrite_scale = cli.rewrite_scale_envelope_plot_scripts
            original_rewrite_no_scale = cli.rewrite_no_scale_ratio_plot_scripts

            def fake_parallel(parallel_args, paths):
                calls.append(("parallel-generation", parallel_args, paths))
                return None

            def fail_polarized_full(_polarized_args):
                raise AssertionError("serial polarized-full should not be used for a parallel standalone campaign")

            def fake_plot_prep(plot_args):
                calls.append(("plot-prep", plot_args))
                return 0

            def fake_run(command, check, env):
                calls.append(("mkhtml", command, check, env))

            def fake_rewrite_scale(plot_dir, python_executable=None, reference_label="POLDIS NLO"):
                calls.append(("postprocess-scale", plot_dir, python_executable, reference_label))
                return 1, 1

            def fake_rewrite_no_scale(plot_dir, python_executable=None):
                calls.append(("postprocess-noscale", plot_dir, python_executable))
                return 0, 0

            try:
                cli._run_parallel_standalone_generation = fake_parallel
                cli.command_polarized_full = fail_polarized_full
                cli.command_plot_prep = fake_plot_prep
                cli.subprocess.run = fake_run
                cli.rewrite_scale_envelope_plot_scripts = fake_rewrite_scale
                cli.rewrite_no_scale_ratio_plot_scripts = fake_rewrite_no_scale
                from contextlib import redirect_stdout

                with redirect_stdout(io.StringIO()):
                    rc = cli.command_standalone_campaign(args)
            finally:
                cli._run_parallel_standalone_generation = original_parallel
                cli.command_polarized_full = original_polarized_full
                cli.command_plot_prep = original_plot_prep
                cli.subprocess.run = original_run
                cli.rewrite_scale_envelope_plot_scripts = original_rewrite_scale
                cli.rewrite_no_scale_ratio_plot_scripts = original_rewrite_no_scale

            self.assertEqual(0, rc)
            self.assertEqual(
                ["parallel-generation", "plot-prep", "mkhtml", "postprocess-scale"],
                [call[0] for call in calls],
            )
            parallel_args = calls[0][1]
            self.assertEqual(13, parallel_args.events)
            self.assertEqual(2, parallel_args.jobs)
            self.assertEqual(3, parallel_args.shards)
            plot_args = calls[1][1]
            self.assertTrue(plot_args.herwig_scale_envelope)

    def test_progress_reporter_emits_fraction_rate_and_eta(self):
        import herwig_fixed_order_nlo as cli

        class FakeClock:
            def __init__(self):
                self.value = 0.0

            def __call__(self):
                return self.value

        clock = FakeClock()
        stream = io.StringIO()
        progress = cli.ProgressReporter(
            "generate ALL/00",
            total=4,
            unit="seeds",
            stream=stream,
            clock=clock,
            min_interval=0.0,
        )

        progress.start()
        clock.value = 2.0
        progress.update(2)
        clock.value = 4.0
        progress.done()
        text = stream.getvalue()

        self.assertIn("[generate ALL/00] start: 4 seeds", text)
        self.assertIn("50.0% (2/4 seeds)", text)
        self.assertIn("1.00 seeds/s", text)
        self.assertIn("ETA 2s", text)
        self.assertIn("[generate ALL/00] done: 4 seeds in 4s", text)

    def test_breit_polarized_equivalents_include_delta_and_all_histograms(self):
        fake_yoda = type("FakeYODA", (), {})()
        fake_poldis = type(
            "FakePoldis",
            (),
            {
                "new_binned_estimate": staticmethod(lambda edges: None),
                "set_bin_val_err": staticmethod(lambda bin_obj, value, error: None),
            },
        )()
        old_yoda = sys.modules.get("yoda")
        old_poldis = sys.modules.get("poldis_top_to_yoda")
        sys.modules["yoda"] = fake_yoda
        sys.modules["poldis_top_to_yoda"] = fake_poldis
        try:
            script = SCRIPT_DIR / "analyze_DIS_polarized.py"
            spec = importlib.util.spec_from_file_location("analyze_dis_polarized_for_breit_all_test", script)
            module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(module)
        finally:
            if old_yoda is None:
                sys.modules.pop("yoda", None)
            else:
                sys.modules["yoda"] = old_yoda
            if old_poldis is None:
                sys.modules.pop("poldis_top_to_yoda", None)
            else:
                sys.modules["poldis_top_to_yoda"] = old_poldis

        labels = set(module.ANALYSIS_LABELS["MC_DIS_BREIT"])
        all_labels = set(module.ANALYSIS_ALL_LABELS["MC_DIS_BREIT"])

        self.assertIn("Q2", labels)
        self.assertIn("Zeta", labels)
        self.assertIn("Q2", all_labels)
        self.assertIn("pT1", all_labels)
        self.assertIn("Zeta", all_labels)

    def test_plot_safe_scatter_symmetrizes_one_sided_errors(self):
        class FakeErrs:
            def __init__(self, minus, plus):
                self.minus = minus
                self.plus = plus

        class FakePoint2D:
            def __init__(self):
                self._x = 0.0
                self._y = 0.0
                self._xerrs = FakeErrs(0.0, 0.0)
                self._yerrs = FakeErrs(0.0, 0.0)

            def setX(self, value):
                self._x = value

            def setY(self, value):
                self._y = value

            def setXErrs(self, minus, plus):
                self._xerrs = FakeErrs(minus, plus)

            def setYErrs(self, minus, plus):
                self._yerrs = FakeErrs(minus, plus)

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
                return self._x - self._xerrs.minus

            def xMax(self):
                return self._x + self._xerrs.plus

        class FakeScatter2D:
            def __init__(self):
                self._path = ""
                self._points = []
                self._annotations = {}

            def setPath(self, path):
                self._path = path

            def path(self):
                return self._path

            def type(self):
                return "Scatter2D"

            def addPoint(self, point):
                self._points.append(point)

            def points(self):
                return list(self._points)

            def point(self, index):
                return self._points[index]

            def numPoints(self):
                return len(self._points)

            def setAnnotation(self, key, value):
                self._annotations[key] = value

            def annotation(self, key):
                return self._annotations[key]

            def annotations(self):
                return list(self._annotations)

        fake_yoda = type("FakeYODA", (), {"Scatter2D": FakeScatter2D, "Point2D": FakePoint2D})()
        fake_poldis = type(
            "FakePoldis",
            (),
            {
                "new_binned_estimate": staticmethod(lambda edges: None),
                "set_bin_val_err": staticmethod(lambda bin_obj, value, error: None),
            },
        )()
        old_yoda = sys.modules.get("yoda")
        old_poldis = sys.modules.get("poldis_top_to_yoda")
        sys.modules["yoda"] = fake_yoda
        sys.modules["poldis_top_to_yoda"] = fake_poldis
        try:
            script = SCRIPT_DIR / "analyze_DIS_polarized.py"
            spec = importlib.util.spec_from_file_location("analyze_dis_polarized_for_plot_safe_test", script)
            module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(module)
        finally:
            if old_yoda is None:
                sys.modules.pop("yoda", None)
            else:
                sys.modules["yoda"] = old_yoda
            if old_poldis is None:
                sys.modules.pop("poldis_top_to_yoda", None)
            else:
                sys.modules["poldis_top_to_yoda"] = old_poldis

        point = FakePoint2D()
        point.setX(150.0)
        point.setY(2.0)
        point.setXErrs(10.0, 10.0)
        point.setYErrs(0.0, 0.4)
        scatter = FakeScatter2D()
        scatter.setPath("/REF/MC_DIS_BREIT/DQ2")
        scatter.addPoint(point)

        safe = module.symmetrize_plot_scatter_objects({"/REF/MC_DIS_BREIT/DQ2": scatter})
        yerrs = safe["/REF/MC_DIS_BREIT/DQ2"].point(0).yErrs()

        self.assertEqual(0.4, yerrs.minus)
        self.assertEqual(0.4, yerrs.plus)
        self.assertEqual("0", safe["/REF/MC_DIS_BREIT/DQ2"].annotation("ErrorBands"))
        self.assertEqual("1", safe["/REF/MC_DIS_BREIT/DQ2"].annotation("ErrorBars"))

    def test_integrate_command_exposes_plain66_window(self):
        import herwig_fixed_order_nlo as cli

        args = cli.parse_args(
            [
                "integrate",
                "--window",
                "plain66",
                "--toy-pdfs",
                "--order",
                "LO",
                "--q2-bins",
                "2",
                "--y-bins",
                "2",
            ]
        )
        run = cli.run_config_from_args(args)

        self.assertEqual(49.0, run.cuts.q2_min)
        self.assertEqual("sigma0", args.observable)
        self.assertGreater(
            cli.integrate_observable(
                run,
                args.setup,
                args.observable,
                args.helicity,
                args.order,
                cli.load_pdf_pair(run.pdf_profile, use_toy=True),
                cli.MatchboxAlphaS(run.alpha_s_mz, run.mz).alpha_s,
                args.q2_sampling,
                args.q2_bins,
                args.y_bins,
                args.xp_bins,
            ),
            0.0,
        )


if __name__ == "__main__":
    unittest.main()
