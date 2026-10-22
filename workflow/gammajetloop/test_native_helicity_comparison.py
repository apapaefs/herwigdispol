#!/usr/bin/env python3
"""Independent controls for phase conversion and coherent spin contractions."""
import cmath
import itertools
import math
import unittest

from compare_native_helicities import COLOUR, HELICITIES, convert_basis, density, require_finite, weighted


class NativeHelicityComparisonTest(unittest.TestCase):
    def test_basis_swap_phases_and_gauge_shifts(self):
        # Circular vectors for four null momenta along the z axis. The expected
        # result uses the prescribed permutation/phases, independently of overlaps.
        pol, native_pol, momenta = {}, {}, {}
        phases = {}
        for leg in range(4):
            sign = 1 if leg % 2 == 0 else -1
            momenta[leg] = [2., 0., 0., 2.*sign]
            for h in (0, 2):
                pol[leg, h] = [0j, 1/math.sqrt(2), (1j if h == 0 else -1j)/math.sqrt(2), 0j]
            for h in (0, 2):
                phases[leg, h] = cmath.exp(1j*(leg+1)*(h+1)/7)
                original = 2-h if leg == 1 else h
                native_pol[leg, h] = [phases[leg, h]*v + (0.07+0.03j)*p
                                       for v, p in zip(pol[leg, original], momenta[leg])]
        amplitudes = {h: complex(i+1, 3-i) for i, h in enumerate(HELICITIES)}
        mg = dict(pol=pol, amp={h: [2*z, 2*z] for h, z in amplitudes.items()}, momenta_used=momenta)
        native = dict(pol=native_pol, momenta=momenta)
        mapped, checks = convert_basis(mg, native)
        for h in HELICITIES:
            original = (h[0], 2-h[1], h[2], h[3])
            expected = amplitudes[original]*math.prod(phases[leg, h[leg]] for leg in range(4))
            self.assertLess(abs(mapped[h]-expected), 2e-14)
        self.assertTrue(all(value < 1e-13 for value in checks.values()))

    def test_complex_incoming_coherence(self):
        # Incoming |psi>=(|->+i|+>)/sqrt(2). Only the equal-helicity
        # incoming states contribute to |gamma-,g->, with relative phase i.
        amp = dict.fromkeys(HELICITIES, 0j)
        amp[0, 0, 0, 0] = 1
        amp[2, 2, 0, 0] = 1j
        rho = {(leg, i, j): 0.5*(1j if i == 2 else 1)*(-1j if j == 2 else 1)
               for leg in range(2) for i, j in itertools.product((0, 2), repeat=2)}
        self.assertAlmostEqual(weighted(amp, rho), COLOUR*abs((1-1j)/2)**2)
        self.assertEqual(density(amp, rho, 0), {(0, 0): 1, (0, 2): 0, (2, 0): 0, (2, 2): 0})

    def test_phase_error_changes_spin_interference_but_not_rate(self):
        # |gamma-,g-> + exp(i phi)|gamma+,g-> has a pure photon density.
        rho = {(leg, i, j): float(i == j == 0)
               for leg in range(2) for i, j in itertools.product((0, 2), repeat=2)}
        amp = dict.fromkeys(HELICITIES, 0j)
        amp[0, 0, 0, 0] = 1
        amp[0, 0, 2, 0] = 1
        baseline = weighted(amp, rho)
        amp[0, 0, 2, 0] = cmath.exp(0.1j)
        self.assertAlmostEqual(weighted(amp, rho), baseline)
        outgoing = density(amp, rho, 0)
        self.assertAlmostEqual(abs(outgoing[0, 2] - cmath.exp(-0.1j)/2), 0.)
        self.assertGreater(abs(outgoing[0, 2]-0.5), 0.04)

    def test_nonfinite_values_are_not_valid_agreement(self):
        for value in (float('nan'), float('inf'), complex(1, float('nan'))):
            with self.assertRaisesRegex(ValueError, 'Non-finite'):
                require_finite([1, value], 'control')


if __name__ == '__main__':
    unittest.main()
