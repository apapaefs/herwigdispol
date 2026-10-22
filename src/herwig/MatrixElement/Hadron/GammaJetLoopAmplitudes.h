// -*- C++ -*-
// Native quark-box amplitudes. GPL version 3; see COPYING.
#ifndef HERWIG_GammaJetLoopAmplitudes_H
#define HERWIG_GammaJetLoopAmplitudes_H

#include <array>
#include <complex>

namespace Herwig { namespace GammaJetLoop {
using Real = long double;
using Complex = std::complex<Real>;
using Vector = std::array<Complex,4>; // (energy, x, y, z), metric +---
using Momenta = std::array<Vector,4>; // all incoming
// Reduced ACS helicities: ++++, -+++, --++, +--+, +-+-.
using HelicityFunctions = std::array<Complex,5>;

HelicityFunctions massless(Real s, Real t, Real u);
HelicityFunctions massive(Real s, Real t, Real u, Real mass);

/** Gauge-invariant tensor reconstructed from the five scalar helicity
 * functions. Contracts actual polarization vectors, including their phases.
 * ACS convention: s=(p1+p2)^2, t=(p2+p3)^2, u=(p1+p3)^2.
 * The coupling/colour factor is 2 sqrt(alphaEM alphaS^3) d^{abc} Qq.
 */
Complex contract(const Momenta &, const std::array<Vector,4> &,
                 const HelicityFunctions &);
Complex dot(const Vector &, const Vector &);
}}
#endif
