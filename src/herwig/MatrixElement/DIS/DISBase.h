// -*- C++ -*-
#ifndef HERWIG_DISBase_H
#define HERWIG_DISBase_H
//
// This is the declaration of the DISBase class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/Shower/ShowerAlpha.h"
#include "ThePEG/EventRecord/RhoDMatrix.h"
#include <map>
#include <string>

namespace Herwig {

using namespace ThePEG;

/**
 * The DISBase class is the base class for the implementation
 * of DIS type processes including corrections in both the old
 * fashioned matrix element and POWHEG approaches
 *
 * @see \ref DISBaseInterfaces "The interfaces"
 * defined for DISBase.
 */
class DISBase: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  DISBase();

  /**
   * The default constructor.
   */
  virtual ~DISBase();

  /**
   *  Members for the old-fashioned matrix element correction
   */
  //@{
  /**
   *  Has an old fashioned ME correction
   */
  virtual bool hasMECorrection() {return true;}

  /**
   *  Initialize the ME correction
   */
  virtual void initializeMECorrection(RealEmissionProcessPtr, double &,
				      double & );

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual RealEmissionProcessPtr applyHardMatrixElementCorrection(RealEmissionProcessPtr);

  /**
   * Apply the soft matrix element correction
   * @param parent The initial particle in the current branching
   * @param progenitor The progenitor particle of the jet
   * @param fs Whether the emission is initial or final-state
   * @param highestpT The highest pT so far in the shower
   * @param ids ids of the particles produced in the branching
   * @param z The momentum fraction of the branching
   * @param scale the evolution scale of the branching
   * @param pT The transverse momentum of the branching
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(PPtr parent,
				     PPtr progenitor,
				     const bool & fs,
				     const Energy & highestpT,
				     const vector<tcPDPtr> & ids,
				     const double & z,
				     const Energy & scale,
				     const Energy & pT);
  //@}

  /**
   *  Members for the POWHEG stype correction
   */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return Both;}

  /**
   *  Apply the POWHEG style correction
   */
  virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr,
						 ShowerInteraction);
  //@}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given nDim() uniform
   * random numbers in the interval \f$ ]0,1[ \f$. To help the phase space
   * generator, the dSigHatDR should be a smooth function of these
   * numbers, although this is not strictly necessary.
   * @param r a pointer to the first of nDim() consecutive random numbers.
   * @return true if the generation succeeded, otherwise false.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * Optional named event weights propagated through ThePEG/HepMC.
   */
  virtual std::map<std::string,double> generateOptionalWeights();
  //@}


public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DISBase & operator=(const DISBase &) = delete;

protected:

  /**
   * Extract the longitudinal polarisation (for spin-1/2 beams) from a density matrix \rho.
   * Definition: P = (Re \rho_{++} − Re \rho_{--}) / (Re \rho_{++} + Re \rho_{--}).
   * Returns 0 if the denominator would be zero.
   */
  double longPol(const ThePEG::RhoDMatrix& rho) const;

  /**
   * Return incoming rho matrices with the hadron-side longitudinal
   * polarization rebuilt from the same PDF ratio used in NLOWeight().
   * This keeps the Born ME and the analytic NLO terms aligned when
   * the reconstructed NLO parton bins do not carry polarized PDF
   * information reliably.
   */
  std::pair<RhoDMatrix,RhoDMatrix> correctedLongitudinalRhoMatrices() const;

  /**
   * Collinear projector weights for the quark and gluon counterterms.
   * The default implementation reproduces the legacy scalar f_pol blend.
   */
  struct CollinearBlendWeights {
    double qUnpolarized;
    double qPolarized;
    double gUnpolarized;
    double gPolarized;
  };

  /**
   * Project the Born angular structure onto the quark/gluon collinear
   * kernels. Neutral-current DIS overrides this to keep the P_q-even and
   * P_q-odd structures separate.
   */
  virtual CollinearBlendWeights collinearBlendWeights(tcPDPtr lin, tcPDPtr lout,
                                                      tcPDPtr qin, tcPDPtr qout,
                                                      Energy2 scale,
                                                      double Pl, double Pq,
                                                      double ell) const;

  /**
   * Ratio of the mapped QCDC parity-even Born coefficient to the underlying
   * Born one, D_m / D_B. This is unity in the legacy/photon case and is
   * overridden for full neutral-current DIS.
   */
  virtual double qcdcMappedDenominatorRatio(tcPDPtr lin, tcPDPtr lout,
                                            tcPDPtr qin, tcPDPtr qout,
                                            Energy2 scale,
                                            double Pl,
                                            double PqBorn,
                                            double PqMapped) const;

  /**
   * Ratio of the real-emission parity-even denominator to the Born one.
   * This is unity in the legacy/photon case and can be overridden by
   * neutral-current implementations that keep the mapped quark/gluon
   * channels separate.
   */
  virtual double realEmissionDenominatorFactor(tcPDPtr lin, tcPDPtr lout,
                                               tcPDPtr qin, tcPDPtr qout,
                                               Energy2 scale,
                                               double Pl,
                                               double mappedPartonPol) const;

  /**
   * Whether the real-emission kernels should use the mapped incoming-parton
   * polarization rather than the Born one.
   */
  virtual bool useMappedPolarizedEmissionKernel() const;

  /**
   * Event-local neutral-current response coefficients used by the common NLO
   * term assembly. The default DISBase fallback is sufficient for photon-like
   * channels; full neutral-current DIS overrides this to keep gamma/Z parity
   * structures separated.
   */
  struct NeutralCurrentResponse {
    double qOddResponse;
    double gOddResponse;
  };

  /**
   * Extract event-local neutral-current response coefficients. Returns false
   * when the common photon-like fallback should be used.
   */
  virtual bool neutralCurrentResponse(tcPDPtr lin, tcPDPtr lout,
                                      tcPDPtr qin, tcPDPtr qout,
                                      Energy2 scale,
                                      double Pl,
                                      double PqBorn,
                                      double PqMapped,
                                      double ell,
                                      NeutralCurrentResponse &out) const;

  /**
   * Return the alphaS scale used in the POWHEG hardest-emission
   * generation. By default this is the native POWHEG transverse-momentum
   * scale pT^2 = q2*xT^2/4. Optionally it can be forced to Q^2 for direct
   * comparison to fixed-order reference calculations.
  */
  Energy2 powhegEmissionAlphaSScale(Energy2 q2, double xT) const;

  /**
   * Return the PDF numerator scale used in the POWHEG/MEC emission kernels.
   * By default this remains the native mapped emission scale. Optionally it
   * can be forced to Q^2 together with the alphaS scale for direct comparison
   * to fixed-order reference calculations.
   */
  Energy2 powhegEmissionPDFScale(Energy2 q2, Energy2 mappedScale) const;

  /**
   * Return the fixed-order alphaS value used for POWHEG/MEC emission
   * comparisons. When the beam PDF is backed by LHAPDF, prefer the
   * PDF-set alphaS to match fixed-order reference programs such as
   * POLDIS as closely as possible. Otherwise fall back to the model
   * alphaS used by the validated fixed-order NLO weight.
   */
  double powhegEmissionFixedOrderAlphaSValue(Energy2 scale) const;

  /**
   * Return the alphaS value used in the POWHEG/MEC emission kernels.
   * The default path uses the configured ShowerAlpha object, while the
   * optional fixed-order mode delegates to the fixed-order alphaS source,
   * preferring LHAPDF alphaS when the beam PDF provides it.
   */
  double powhegEmissionAlphaSValue(Energy2 scale) const;

  /**
   * Return the overestimate value for the POWHEG emission coupling.
   * In the fixed-order alphaS mode this is chosen from the fixed-order
   * alphaS at a reference scale so that the existing veto-generation
   * structure remains unchanged.
   */
  double powhegEmissionAlphaSOverestimate(Energy2 referenceScale) const;

  /**
   * Return the ratio of the POWHEG emission coupling to its overestimate.
   */
  double powhegEmissionAlphaSRatio(Energy2 scale,
                                   Energy2 referenceScale) const;

  /**
   * Return the configured POWHEG emission comparison mode.
   */
  unsigned int powhegEmissionComparisonMode() const {
    return powhegEmissionComparisonMode_;
  }

  /**
   * Return the maximum number of real-emission attempts in comparison mode.
   */
  unsigned long powhegEmissionComparisonMaxAttempts() const {
    return powhegEmissionComparisonMaxAttempts_;
  }

  /**
   * Attach a spin-only hard vertex to the POWHEG real-emission event
   * record. The default implementation is a no-op so existing validated
   * subtraction and emission generation remain unchanged.
   */
  virtual void constructRealEmissionSpinVertex(RealEmissionProcessPtr proc,
                                               bool isCompton) const;

  /**
   * Structured result for mapped incoming-parton polarization and the rho
   * matrix derived from it for the real-emission spin machinery.
   */
  struct MappedIncomingSpinDensity {
    RhoDMatrix rho;
    double clampedPolarization;
  };

  /**
   * Build the mapped incoming-parton rho matrix used by the real-emission
   * spin machinery. When the mapped polarization is numerically unsafe, fall
   * back to the unpolarized rho matrix for the requested parton.
   */
  MappedIncomingSpinDensity mappedIncomingSpinDensity(tcPDPtr parton,
                                                      double x,
                                                      Energy2 scale,
                                                      const char * clampSource = nullptr) const;

  /**
   * Construct a longitudinal rho matrix for spin-1/2 and spin-1 particles.
   * Any unsupported particle is returned in the unpolarized state.
   */
  RhoDMatrix longitudinalRhoMatrix(tcPDPtr data, double pol) const;

  /**
   * Longitudinal polarisation of the mapped incoming parton at the
   * real-emission momentum fraction x. Returns zero if the required
   * polarized PDF information is unavailable.
   */
  double mappedIncomingLongitudinalPolarization(tcPDPtr parton,
                                                double x,
                                                Energy2 scale) const;

protected:

  /**
   * True for charged-current DIS matrix elements and false for neutral-current
   * ones. Used when selecting the applicable SimpleDISCut window.
   */
  virtual bool usesChargedCurrentDISWindow() const = 0;

  /**
   * Whether the exact spin-only POWHEG real-emission vertex is enabled.
   */
  bool usePOWHEGRealSpinVertex() const { return usePOWHEGRealSpinVertex_; }


  /**
   *  The NLO weight
   */
  double NLOWeight() const;

  /**
   * Signed NLO density ratio before POS/NEG clipping. When overridePolarizations
   * is true, the supplied beam polarizations are used instead of the event setup.
   */
  double NLOWeightRaw(double leptonPolarization,
                      double hadronPolarization,
                      bool overridePolarizations) const;

  /**
   * Effective Born-side parton polarization built from the same polarized-PDF
   * ratio used by the NLO terms.
   */
  double rivetWeightBornPartonPolarization(double hadronPolarization) const;

  /**
   * Born matrix element for a requested lepton and parton polarization. The
   * neutral-current implementation overrides this with the exact NC structure.
   */
  virtual double rivetWeightBornME2(double leptonPolarization,
                                    double partonPolarization) const;

  /**
   * Build correlated helicity weights for Rivet.
   */
  std::map<std::string,double> generateRivetWeights() const;

  /**
   *  Calculate the coefficient A for the correlations
   */
  virtual double A(tcPDPtr lin, tcPDPtr lout, tcPDPtr qin, tcPDPtr qout,
		   Energy2 scale) const =0;

  /**
    * Exact analysing power for longitudinally polarised beams.
    * Implemented in neutral-current ME to follow eq. (4.50); default
    * falls back to the unpolarised A(...).
    * @param lin   incoming lepton PD
    * @param lout  outgoing lepton PD
    * @param qin   incoming quark  PD
    * @param qout  outgoing quark  PD
    * @param scale Q^2
    * @param Pl    lepton longitudinal polarisation
    * @param Pq    quark  longitudinal polarisation
    */
    virtual double A_pol(tcPDPtr lin, tcPDPtr lout,
                         tcPDPtr qin, tcPDPtr qout,
                         Energy2 scale, double /*Pl*/, double /*Pq*/) const {
      return A(lin,lout,qin,qout,scale);
    }

  /**
   *  Members for the matrix element correction
   */
  //@{
  /**
   *  Generate the values of \f$x_p\f$ and \f$z_p\f$
   * @param xp The value of xp, output
   * @param zp The value of zp, output
   */
  double generateComptonPoint(double &xp, double & zp);

  /**
   *  Generate the values of \f$x_p\f$ and \f$z_p\f$
   * @param xp The value of xp, output
   * @param zp The value of zp, output
   */
  double generateBGFPoint(double &xp, double & zp);

  /**
   * Compact azimuthal kernel representation
   * c0 + c1 cos(phi) + c2 cos^2(phi).
   */
  struct AzimuthalKernelCoefficients {
    double c0;
    double c1;
    double c2;

    AzimuthalKernelCoefficients(): c0(0.0), c1(0.0), c2(0.0) {}

    double average() const { return c0 + 0.5 * c2; }
  };

  /**
   *  Return the coefficients for the matrix element piece for
   *  the QCD compton case. The output is the \f$a_i\f$ coefficients to 
   *  give the function as 
   *  \f$a_0+a_1\cos\phi+a_2\sin\phi+a_3\cos^2\phi+a_4\sin^2\phi\f$
   * @param xp \f$x_p\f$
   * @param x2 \f$x_2\f$
   * @param xperp \f$x_\perp\f$
   * @param norm Normalise to the large $l$ value of the ME
   */
  AzimuthalKernelCoefficients ComptonME(double xp, double x2, double xperp,
			                bool norm) const;
  
  /**
   *  Return the coefficients for the matrix element piece for
   *  the QCD compton case. The output is the \f$a_i\f$ coefficients to 
   *  give the function as 
   *  \f$a_0+a_1\cos\phi+a_2\sin\phi+a_3\cos^2\phi+a_4\sin^2\phi\f$
   * @param xp \f$x_p\f$
   * @param x2 \f$x_3\f$
   * @param x3 \f$x_2\f$
   * @param xperp \f$x_\perp\f$
   * @param norm Normalise to the large $l$ value of the ME
   */
  AzimuthalKernelCoefficients BGFME(double xp, double x2, double x3,
                                    double xperp, bool norm) const;
  //@}

  /**
   *  Members for the POWHEG correction
   */
  //@{
  /**
   *  Generate a Compton process
   */
  void generateCompton();

  /**
   *  Generate a BGF process
   */
  void generateBGF();

  enum POWHEGEmissionComparisonMode {
    POWHEGEmissionComparisonModeDefault = 0,
    POWHEGEmissionComparisonModeRealOnly = 1
  };

  /**
   * Prepare the common Born-level state used by the comparison-mode emission
   * samplers.
   */
  void initializePOWHEGEmissionState(RealEmissionProcessPtr born,
                                     PPtr quark[2], PPtr lepton[2],
                                     PPtr & hadron,
                                     unsigned int & iqIn,
                                     unsigned int & iqOut);

  /**
   * Native POWHEG hardest-emission generation with optional Born fallback.
   */
  RealEmissionProcessPtr generateNativePOWHEGHardest(RealEmissionProcessPtr born,
                                                     bool allowBornFallback);

  /**
   * Generate the configured comparison-mode hardest emission.
   */
  RealEmissionProcessPtr generateComparisonModePOWHEGHardest(RealEmissionProcessPtr born);

  /**
   *  Parameters for the matrix element correction
   */
  //@{
  /**
   *  Enchancement factor for ISR
   */
  double initial_;

  /**
   *  Enchancement factor for FSR
   */
  double final_;

  /**
   *   Relative fraction of compton and BGF processes to generate
   */
  double procProb_;

  /**
   *  Integral for compton process
   */
  double comptonInt_;

  /**
   *  Integral for BGF process
   */
  double bgfInt_;
  //@}

  /**
   *  Parameters for the POWHEG correction
   */
  //@{
  /**
   *  Weight for the compton channel
   */
  double comptonWeight_;

  /**
   *  Weight for the BGF channel
   */
  double BGFWeight_;

  /**
   *  Minimum value of \f$p_T\f$
   */
  Energy pTmin_;
  //@}

  /**
   *  Parameters for the point being generated
   */
  //@{
  /**
   *   \f$Q^2\f$
   */
  Energy2 q2_;

  /**
   *  
   */
  double l_;

  /**
   *  Borm momentum fraction
   */
  double xB_;

  /**
   *  Beam particle
   */
  tcBeamPtr beam_;

  /**
   *  Partons
   */
  tcPDPtr partons_[2];

  /**
   *  Leptons
   */
  tcPDPtr leptons_[2];

  /**
   *  PDF object
   */
  tcPDFPtr pdf_;
  /**
   *  Rotation to the Breit frame
   */
  LorentzRotation rot_;

  /**
   *  Lepton momenta
   */
  Lorentz5Momentum pl_[2];

  /**
   *  Quark momenta
   */
  Lorentz5Momentum pq_[2];

  /**
   *  q
   */
  Lorentz5Momentum q_;

  /**
   *  Compton parameters
   */
  Energy pTCompton_;
  bool ComptonISFS_;
  vector<Lorentz5Momentum> ComptonMomenta_;

  /**
   *  BGF parameters
   */
  Energy pTBGF_;
  vector<Lorentz5Momentum> BGFMomenta_;
  //@}

  /**
   *  The coefficient for the correlations
   */
  double acoeff_;

  /**
   *  Coupling
   */
  ShowerAlphaPtr alpha_;

  /**
   *  Gluon particle data object
   */
  PDPtr gluon_;

private:

  /**
   * Cached copy of the active SimpleDISCut window for native DIS generation.
   */
  struct NativeDISWindowDefinition {
    NativeDISWindowDefinition()
      : resolutionAttempted(false), available(false),
        minQ2(ZERO), maxQ2(ZERO), minY(0.0), maxY(1.0),
        minW2(ZERO), maxW2(ZERO) {}

    bool resolutionAttempted;
    bool available;
    Energy2 minQ2;
    Energy2 maxQ2;
    double minY;
    double maxY;
    Energy2 minW2;
    Energy2 maxW2;
  };

  /**
   * Resolve and cache the active SimpleDISCut window for the current DIS
   * process type. Returns false when no unique matching window is available.
   */
  bool resolveNativeDISWindow() const;

  /**
   * Identify the hadron beam and its Bjorken-x for the current Born point.
   */
  bool determineBornHadronAndXB(tcBeamPtr & hadron, double & xB) const;

  /**
   * Tighten the Born cos(theta) interval with the native DIS window.
   */
  bool tightenBornCosThetaWithNativeDISWindow(const HwMEBase::TwoToTwoKinematicsSetup & setup,
                                              double xB,
                                              double & ctmin,
                                              double & ctmax) const;

  /**
   *  The radiative variables
   */
  //@{
  /**
   *  The \f$x_p\f$ or \f$z\f$ real integration variable
   */
  double xp_;
  //@}

  /**
   *  The hadron
   */
  tcBeamPtr hadron_;

  /**
   * Selects a dynamic or fixed factorization scale
   */
  unsigned int scaleOpt_;

  /**
   * The factorization scale 
   */
  Energy muF_;

  /**
   *  Prefactor if variable scale used
   */
  double scaleFact_;

  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Use the same alphaS implementation as the fixed-order NLO code in the
   *  POWHEG/MEC emission kernels.
   */
  bool useFixedOrderAlphaSInPOWHEGEmission_;

  /**
   *  Use Q^2 instead of the native POWHEG pT^2 scale in the POWHEG/MEC
   *  emission kernels, both for the alphaS evaluation and the emission-PDF
   *  numerator scale.
   */
  bool useQ2ScaleInPOWHEGEmission_;

  /**
   *  Attach correlated helicity weights for Rivet analyses.
   */
  bool generateRivetWeights_;

  /**
   *  Tighten the Born generation to the full DIS window before applying the
   *  safety-check passCuts() veto.
   */
  bool useNativeDISWindowGeneration_;

  /**
   *  Comparison-only ladder for bringing the POWHEG hardest-emission path
   *  closer to a fixed-order real-emission generator.
   */
  unsigned int powhegEmissionComparisonMode_;

  /**
   *  Maximum number of real-emission attempts before vetoing the event in a
   *  non-default POWHEG emission comparison mode.
   */
  unsigned long powhegEmissionComparisonMaxAttempts_;

  /**
   *  Attach an exact spin-only HardVertex to POWHEG real-emission states.
   */
  bool usePOWHEGRealSpinVertex_;

  /**
   *  Lepton longitudinal polarization for the current real-emission state.
   */
  double leptonPolarization_;

  /**
   *  Channel-local POWHEG raw variables for the currently selected winner.
   *  These transient generation values are rebuilt for each event.
   */
  double comptonRawXP_, comptonRawZP_, bgfRawXP_, bgfRawZP_;

  /**
   *  Sampling state for the mapped xp variable.
   */
  double xpSamplingRandom_, xpSamplingRho_, xpSamplingRhomin_;

  /**
   *  Cached copy of the active DIS window used by native generation.
   */
  mutable NativeDISWindowDefinition nativeDISWindow_;

  /**
   *  Power for sampling \f$x_p\f$
   */
  double power_;

  /**
   *  Jacobian for \f$x_p\f$ integral
   */
  double jac_;

};

}

#endif /* HERWIG_DISBase_H */
