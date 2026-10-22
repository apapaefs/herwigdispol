// -*- C++ -*-
// Native loop-induced photon+jet production. GPL version 3; see COPYING.
#ifndef HERWIG_MEPP2GammaJetLoop_H
#define HERWIG_MEPP2GammaJetLoop_H
#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include <array>

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/** Separate native gg -> gamma g quark-box process, including hard spin
 * correlations. The diagnostic Flat mode is not a physics prediction.
 * Use QCD showers: the neutral hard state has no QED shower partner.
 */
class MEPP2GammaJetLoop: public HwMEBase {
public:
  MEPP2GammaJetLoop();
  unsigned int orderInAlphaS() const override {return 3;}
  unsigned int orderInAlphaEW() const override {return 1;}
  double me2() const override;
  Energy2 scale() const override;
  void getDiagrams() const override;
  Selector<DiagramIndex> diagrams(const DiagramVector &) const override;
  Selector<const ColourLines *> colourGeometries(tcDiagPtr) const override;
  void constructVertex(tSubProPtr) override;
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  static void Init();
protected:
  IBPtr clone() const override {return new_ptr(*this);}
  IBPtr fullclone() const override {return new_ptr(*this);}
private:
  ProductionMatrixElement amplitudes(
    const std::array<Lorentz5Momentum,4> &,
    const std::array<vector<VectorWaveFunction>,4> &) const;
  pair<RhoDMatrix,RhoDMatrix> incomingDensities() const;
  unsigned int amplitudeMode_;
  double scalePreFactor_;
};
}
#endif
