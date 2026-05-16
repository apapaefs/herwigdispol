// -*- C++ -*-
#ifndef HERWIG_MEChargedCurrentDIS_H
#define HERWIG_MEChargedCurrentDIS_H
//
// This is the declaration of the MEChargedCurrentDIS class.
//

#include "DISBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Charged-current DIS matrix element with W exchange, exact longitudinal
 * spin-density support, charged-current chiral NLO factors, and optional
 * spin-only POWHEG real-emission vertices.
 *
 * By default both the incoming and outgoing quarks are treated as massless.
 * The outgoing-quark mass option is retained for heavy-flavour production,
 * especially top channels.
 *
 * @see \ref MEChargedCurrentDISInterfaces "The interfaces"
 * defined for MEChargedCurrentDIS.
 */
class MEChargedCurrentDIS: public DISBase {

public:

  /**
   * Construct with the validated charged-current defaults.
   */
  MEChargedCurrentDIS();

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Add all allowed W-mediated lepton/quark Born diagrams.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   *  Construct the Born spin-correlation vertex for the hard subprocess.
   */
  virtual void constructVertex(tSubProPtr);
  //@}


public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Write the charged-current persistent state.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Read the charged-current persistent state.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   * Charged-current DIS uses the charged-current SimpleDISCut window.
   */
  virtual bool usesChargedCurrentDISWindow() const override { return true; }

  /**
   * Coherent helicity matrix element for \f$\ell q\to W^\pm \to \ell q\f$.
   * @param rhoin Rho matrices for incoming particles
   * @param f1 Fermion on lepton line
   * @param a1 Anti-fermion on lepton line
   * @param f2 Fermion on quark line
   * @param a2 Anti-fermion on quark line
   * @param lorder The order of particles on the lepton line
   * @param qorder The order of particles on the quark line
   * @param me Whether to cache the matrix element for spin correlations.
   */
  double helicityME(const pair<RhoDMatrix,RhoDMatrix> & rhoin,
		    vector<SpinorWaveFunction>    & f1 ,
		    vector<SpinorWaveFunction>    & f2,
		    vector<SpinorBarWaveFunction> & a1 ,
		    vector<SpinorBarWaveFunction> & a2,
		    bool lorder, bool qorder,
		    bool me) const;

  /**
   *  Charged-current analyzing coefficient for hard-radiation kernels.
   */
  virtual double A(tcPDPtr lin, tcPDPtr lout, tcPDPtr qin, tcPDPtr qout,
		   Energy2 scale) const;

  /**
   * Charged-current analyzing power is independent of incoming longitudinal
   * polarizations once the chiral hadron factor is handled separately, so this
   * reduces to the charge-conjugation-dependent unpolarized value.
   */
  virtual double A_pol(tcPDPtr lin, tcPDPtr lout,
                       tcPDPtr qin, tcPDPtr qout,
                       Energy2 scale, double Pl, double Pq) const override;

  /**
   * Charged-current DIS factorises the hadron-spin dependence into an overall
   * chiral prefactor on the quark line, so the NLO quark/gluon projectors must
   * be built from that prefactor rather than from A_pol().
   */
  virtual CollinearBlendWeights collinearBlendWeights(tcPDPtr lin, tcPDPtr lout,
                                                      tcPDPtr qin, tcPDPtr qout,
                                                      Energy2 scale,
                                                      double Pl, double Pq,
                                                      double ell) const override;

  /**
   * Ratio of the mapped to Born charged-current quark-line spin prefactors.
   */
  virtual double qcdcMappedDenominatorRatio(tcPDPtr lin, tcPDPtr lout,
                                            tcPDPtr qin, tcPDPtr qout,
                                            Energy2 scale,
                                            double Pl,
                                            double PqBorn,
                                            double PqMapped) const override;

  /**
   * Charged-current real-emission kernels use the same mapped chiral
   * prefactor as the Born and collinear terms.
   */
  virtual double realEmissionDenominatorFactor(tcPDPtr lin, tcPDPtr lout,
                                               tcPDPtr qin, tcPDPtr qout,
                                               Energy2 scale,
                                               double Pl,
                                               double Pq) const override;

  /**
   * Accepted-emission kernels should follow the mapped x_B/x_p polarization.
   */
  virtual bool useMappedPolarizedEmissionKernel() const override;

  /**
   * Attach an exact spin-only HardVertex to the realised POWHEG 2->3 state.
   * This does not alter generated kinematics or event weights.
   */
  virtual void constructRealEmissionSpinVertex(RealEmissionProcessPtr proc,
                                               bool isCompton) const override;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * Exact helicity amplitude for the realised QCDC POWHEG configuration.
   * The matrix element is used only for the spin-correlation HardVertex.
   */
  ProductionMatrixElement qcdcRealEmissionME(PPtr lin, PPtr qin,
                                             PPtr lout, PPtr qout,
                                             PPtr gout, Energy2 q2) const;

  /**
   * Exact helicity amplitude for the realised BGF POWHEG configuration.
   * The matrix element is used only for the spin-correlation HardVertex.
   */
  ProductionMatrixElement bgfRealEmissionME(PPtr lin, PPtr gin,
                                            PPtr lout, PPtr qout,
                                            PPtr qbout, Energy2 q2) const;

  /**
   * Access the QCD quark-gluon vertex used by the exact POWHEG spin amplitudes.
   */
  AbstractFFVVertexPtr gluonVertex() const;

  /**
   * Charged-current quark-line chiral prefactor on the hadron side.
   */
  double ccHadronSpinFactor(tcPDPtr qin, double Pq) const;

  /**
   * Evaluate the Born matrix element for explicit longitudinal beam and
   * parton polarizations.
   */
  double me2ForPolarizations(double Pl, double Pq) const;

  /**
   * Optional Rivet-weight generation uses the same explicit Born helicity path.
   */
  virtual double rivetWeightBornME2(double Pl, double Pq) const override;

  /**
   * Return the ThePEG propagator option for spacelike W exchange in the
   * helicity-amplitude path.
   */
  int wHelicityPropagatorOption(int legacyOption) const {
    return _useFiniteWidthSpacelikeWPropagator ? 7 : legacyOption;
  }

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEChargedCurrentDIS & operator=(const MEChargedCurrentDIS &) = delete;

private:

  /**
   *  W-fermion vertex used by Born and real-emission helicity amplitudes.
   */
  AbstractFFVVertexPtr _theFFWVertex;

  /**
   *  Heaviest incoming quark flavour allowed in charged-current diagrams.
   */
  unsigned int _maxflavour;

  /**
   *  Option controlling whether outgoing quark masses are retained.
   */
  unsigned int _massopt;

  /**
   * Whether to keep the finite W width for spacelike charged-current exchange.
   */
  bool _useFiniteWidthSpacelikeWPropagator;

  /**
   * Born production matrix element cached for spin correlations.
   */
  ProductionMatrixElement _me;

  /**
   *  Intermediate charged-current boson data.
   */
  //@{
  /**
   *  W+ ParticleData object.
   */
  tcPDPtr _wp;

  /**
   *  W- ParticleData object.
   */
  tcPDPtr _wm;
  //@}
};

}

#endif /* HERWIG_MEChargedCurrentDIS_H */
