// -*- C++ -*-
#ifndef HERWIG_MENeutralCurrentDIS_H
#define HERWIG_MENeutralCurrentDIS_H
//
// This is the declaration of the MENeutralCurrentDIS class.
//

#include "DISBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Neutral-current DIS matrix element with coherent photon/Z exchange,
 * exact longitudinal spin-density support, neutral-current NLO response
 * coefficients, and optional spin-only POWHEG real-emission vertices.
 * This is where the gamma/Z interference and parity bookkeeping live; the
 * common DISBase code only sees compact response factors.
 *
 * For consistency both the incoming and outgoing quarks are assumed to be
 * massless, matching the analytic DIS NLO and real-emission kernels.
 *
 * @see \ref MENeutralCurrentDISInterfaces "The interfaces"
 * defined for MENeutralCurrentDIS.
 */
class MENeutralCurrentDIS: public DISBase {

public:

  /**
   * Construct with the validated neutral-current defaults.
   */
  MENeutralCurrentDIS();

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
   * Add all allowed lepton/quark photon and Z Born diagrams.
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
   * Write the neutral-current persistent state.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Read the neutral-current persistent state.
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
   * Neutral-current DIS uses the neutral-current SimpleDISCut window.
   */
  virtual bool usesChargedCurrentDISWindow() const override { return false; }

  /**
   * Coherent helicity matrix element for \f$\ell q\to \gamma/Z \to \ell q\f$.
   * @param rhoin Rho matrices for incoming particles
   * @param f1 Fermion on lepton line
   * @param a1 Anti-fermion on lepton line
   * @param f2 Fermion on quark line
   * @param a2 Anti-fermion on quark line
   * @param lorder The order of particles on the lepton line
   * @param qorder The order of particles on the quark line
   * @param me Whether to cache the matrix element for spin correlations.
   */
  double helicityME(const pair<RhoDMatrix,RhoDMatrix> & rhoin ,
                    vector<SpinorWaveFunction>    & f1 ,
		    vector<SpinorWaveFunction>    & f2,
		    vector<SpinorBarWaveFunction> & a1 ,
		    vector<SpinorBarWaveFunction> & a2,
		    bool lorder, bool qorder,
		    bool me) const;


  /**
   *  Option for treatment of \f$\gamma/Z\f$ terms.
   */
  inline unsigned int gammaZOption() const {return _gammaZ;}

  /**
   *  Unpolarized analyzing coefficient for the neutral-current MEC kernels.
   */
  virtual double A(tcPDPtr lin, tcPDPtr lout, tcPDPtr qin, tcPDPtr qout,
		   Energy2 scale) const;

   /**
    * Exact analyzing power for longitudinally polarized beams in NC DIS.
    * Implements the retained gamma/Z parity structure and falls back to
    * the photon-only expression when only photon exchange is selected.
    *
    * @param lin   incoming lepton PD
    * @param lout  outgoing lepton PD
    * @param qin   incoming quark  PD
    * @param qout  outgoing quark  PD
    * @param scale Q^2
    * @param Pl    lepton longitudinal polarisation
    * @param Pq    quark  longitudinal polarisation
    * @return      A_pol(Pl,Pq) for NC DIS
    */
   virtual double A_pol(tcPDPtr lin, tcPDPtr lout,
                        tcPDPtr qin, tcPDPtr qout,
                        Energy2 scale, double Pl, double Pq) const override;

  /**
   * Attach an exact spin-only HardVertex to the realised POWHEG 2->3 state.
   * This does not change the generated kinematics or event weight.
   */
  virtual void constructRealEmissionSpinVertex(RealEmissionProcessPtr proc,
                                               bool isCompton) const override;

  /**
   * Exact collinear projectors for NC DIS. These separate the Born pieces
   * that multiply unpolarized kernels (F2/F3-like) from the pieces that
   * multiply polarized kernels (G2/G4-like), and remove the gluon projector
   * for the parity-odd spin-even terms that have no gluon kernel at this
   * order.
   */
  virtual CollinearBlendWeights collinearBlendWeights(tcPDPtr lin, tcPDPtr lout,
                                                      tcPDPtr qin, tcPDPtr qout,
                                                      Energy2 scale,
                                                      double Pl, double Pq,
                                                      double ell) const override;

  /**
   * Exact QCDC mapped/Born parity-even ratio D_m / D_B for NC DIS.
   */
  virtual double qcdcMappedDenominatorRatio(tcPDPtr lin, tcPDPtr lout,
                                            tcPDPtr qin, tcPDPtr qout,
                                            Energy2 scale,
                                            double Pl,
                                            double PqBorn,
                                            double PqMapped) const override;

  /**
   * Exact parity-even denominator for the mapped NC real-emission kernel.
   */
  virtual double realEmissionDenominatorFactor(tcPDPtr lin, tcPDPtr lout,
                                               tcPDPtr qin, tcPDPtr qout,
                                               Energy2 scale,
                                               double Pl,
                                               double mappedPartonPol) const override;

  /**
   * Neutral-current accepted-emission kernels should use the mapped
   * x_B/x_p polarization consistently with the NLO weight.
   */
  virtual bool useMappedPolarizedEmissionKernel() const override;

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

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

private:

  /**
   * Neutral-current coupling coefficients entering the exact Born angular
   * structure and analyzing powers, including the finite-width Z propagator.
   * The N entries multiply the ell term; the D entries multiply the
   * 1 + ell^2 term. Suffixes mark which incoming spin factors are present.
   */
  struct NCCoefficients {
    double N0, Nl, Nq, Nlq;
    double D0, Dl, Dq, Dlq;
  };

  /**
   * Build the neutral-current coefficient set for the given incoming lepton
   * and quark species at momentum transfer q2. Keeping the eta signs here
   * prevents the NLO formulas from scattering particle/antiparticle rules.
   */
  NCCoefficients ncCoefficients(tcPDPtr lin, tcPDPtr qin, Energy2 q2) const;

  /**
   * Evaluate the Born ME with the current kinematics but custom incoming
   * lepton and quark longitudinal polarizations.
   */
  virtual double rivetWeightBornME2(double Pl, double Pq) const override;

  /**
   * Evaluate the Born ME with the current kinematics but custom incoming
   * lepton and quark longitudinal polarizations.
   */
  double me2ForPolarizations(double Pl, double Pq) const;

  /**
   * Evaluate the Born ME with the current kinematics but a custom incoming
   * quark longitudinal polarization and the current incoming lepton
   * polarization.
   */
  double me2ForPartonPolarization(double Pq) const;

  /**
   * Exact Born spin factor Sigma_B for the current NC couplings and a given
   * incoming quark polarization.
   */
  double sigmaBornFactor(tcPDPtr lin, tcPDPtr qin, Energy2 q2,
                         double Pl, double Pq, double ell) const;

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
   * Return the ThePEG propagator option for the spacelike Z boson used in
   * the helicity-amplitude path.
   */
  int zHelicityPropagatorOption() const {
    return _useFiniteWidthSpacelikeZPropagator ? 7 : 1;
  }

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MENeutralCurrentDIS & operator=(const MENeutralCurrentDIS &) = delete;

private:

  /**
   *  Vertices used by the Born and real-emission helicity amplitudes.
   */
  //@{
  /**
   *  Z-fermion vertex.
   */
  AbstractFFVVertexPtr _theFFZVertex;

  /**
   *  Photon-fermion vertex.
   */
  AbstractFFVVertexPtr _theFFPVertex;

  /**
   *  Gluon-fermion vertex for real-emission spin amplitudes.
   */
  AbstractFFVVertexPtr _theFFGVertex;
  //@}

  /**
   *  Intermediate neutral-current boson data.
   */
  //@{
  /**
   *  Z ParticleData object.
   */
  tcPDPtr _z0;

  /**
   *  Photon ParticleData object.
   */
  tcPDPtr _gamma;
  //@}

  /**
   *  Switches controlling flavour coverage and exchanged bosons.
   */
  //@{
  /**
   *  Minimum incoming quark flavour.
   */
  int _minflavour;

  /**
   *  Maximum incoming quark flavour.
   */
  int _maxflavour;

  /**
   *  Whether to include both \f$Z^0\f$ and \f$\gamma\f$ or only one.
   */
  unsigned int _gammaZ;

  /**
   * Whether to keep a finite width for spacelike Z exchange in both the
   * helicity amplitudes and analytic Born/NLO projectors.
   */
  bool _useFiniteWidthSpacelikeZPropagator;
  //@}

  /**
   * Born production matrix element cached for spin correlations.
   */
  ProductionMatrixElement _me;

  /**
   *  Cached electroweak parameters.
   */
  //@{
  /**
   *  \f$\sin\theta_W\f$
   */
  double _sinW;

  /**
   *  \f$\cos\theta_W\f$
   */
  double _cosW;

  /**
   *  The square of the Z mass
   */
  Energy2 _mz2;
  //@}

};

}

#endif /* HERWIG_MENeutralCurrentDIS_H */
