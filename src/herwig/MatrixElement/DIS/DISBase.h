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
 * Common base for DIS matrix elements. It owns the shared Born/NLO
 * kinematics, the legacy hard/soft matrix-element correction, the
 * POWHEG hardest-emission path, fixed-order comparison controls, and
 * correlated helicity weights. Neutral- and charged-current subclasses
 * provide the current-specific Born and real-emission response factors.
 *
 * In practice this is the shared map of a DIS event: which density was
 * sampled, which scales were used, and how optional helicity weights sit on
 * top of that event. The subclasses should only need to answer what current
 * was exchanged and how its helicity algebra responds.
 *
 * @see \ref DISBaseInterfaces "The interfaces"
 * defined for DISBase.
 */
class DISBase: public HwMEBase {

public:

  /**
   * Construct with the validated release defaults for DIS generation.
   */
  DISBase();

  /**
   * Virtual destructor for the polymorphic matrix-element base.
   */
  virtual ~DISBase();

  /**
   *  Members for the legacy shower matrix-element correction.
   */
  //@{
  /**
   *  This matrix element supplies the legacy DIS hard/soft MEC hooks.
   */
  virtual bool hasMECorrection() {return true;}

  /**
   *  Extract the Born DIS state needed by the legacy MEC sampler.
   */
  virtual void initializeMECorrection(RealEmissionProcessPtr, double &,
				      double & );

  /**
   *  Try to replace the Born DIS process by one accepted hard real emission.
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
   *  Members for the POWHEG-style correction.
   */
  //@{
  /**
   *  DIS can generate both initial- and final-state POWHEG emissions.
   */
  virtual POWHEGType hasPOWHEGCorrection() {return Both;}

  /**
   *  Generate the hardest POWHEG emission or return/veto according to mode.
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
   * Optional named event weights propagated through ThePEG/HepMC. These are
   * multipliers for the already-sampled event, not standalone cross sections.
   */
  virtual std::map<std::string,double> generateOptionalWeights();
  //@}


public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Write the current release-facing persistent layout.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Read the current release-facing persistent layout.
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
   * Polarization inputs for the physical extractor side containing hadron_.
   * A non-unique match is deliberately returned as unpolarized with no
   * difference PDF.
   */
  struct HadronBeamPolarizationData {
    bool matched;
    double longitudinal;
    tcPDFPtr differencePDF;
  };

  HadronBeamPolarizationData hadronBeamPolarizationData() const;

  /**
   * Return incoming rho matrices with the hadron-side longitudinal
   * polarization rebuilt from the same PDF ratio used in NLOWeight().
   * This keeps the Born ME and the analytic NLO terms aligned when
   * the reconstructed NLO parton bins do not carry polarized PDF
   * information reliably.
   */
  std::pair<RhoDMatrix,RhoDMatrix> correctedLongitudinalRhoMatrices() const;

  /**
   * Reduced collinear projectors for the quark and gluon counterterms.
   * The polarized entries multiply Pz*Delta f/q directly and therefore
   * remain finite when the Born difference PDF crosses zero.
   */
  struct CollinearBlendWeights {
    double qUnpolarized;
    double qPolarizedReduced;
    double gUnpolarized;
    double gPolarizedReduced;
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
   * Return the alphaS scale used in the POWHEG hardest-emission
   * generation. By default this is the native POWHEG transverse-momentum
   * scale pT^2 = q2*xT^2/4. Optionally it can be forced to Q^2 for direct
   * comparison to fixed-order reference calculations. The configured
   * ScaleFactor wraps either choice so variations remain a common convention.
  */
  Energy2 powhegEmissionAlphaSScale(Energy2 q2, double xT) const;

  /**
   * Return the PDF numerator scale used in the POWHEG/MEC emission kernels.
   * By default this remains the native mapped emission scale. Optionally it
   * can be forced to Q^2 together with the alphaS scale for direct comparison
   * to fixed-order reference calculations. ScaleFactor is applied after that
   * central-scale choice, matching powhegEmissionAlphaSScale().
   */
  Energy2 powhegEmissionPDFScale(Energy2 q2, Energy2 mappedScale) const;

  /**
   * Apply the optional common lower boundary to a squared perturbative or
   * PDF scale.  The default zero boundary leaves all historical choices
   * unchanged.
   */
  Energy2 flooredScale(Energy2 value) const;

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
   * Positive POS/NEG contribution weight derived from NLOWeightRaw().
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
   * Build correlated helicity weights for Rivet. The weights are normalized to
   * the sampled event density; ThePEG applies the carrier event weight later.
   */
  std::map<std::string,double> generateRivetWeights() const;

  /**
   * Born-level longitudinal analyzing coefficient for the current.
   */
  virtual double A(tcPDPtr lin, tcPDPtr lout, tcPDPtr qin, tcPDPtr qout,
		   Energy2 scale) const =0;

  /**
    * Exact analyzing power for longitudinally polarized beams. Subclasses
    * override when the current carries nontrivial parity structure; the base
    * fallback preserves the historical unpolarized A(...) behaviour.
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
   *  Members for the matrix-element correction.
   */
  //@{
  /**
   * Sample the QCDC/Compton real-emission variables for the MEC path.
   * @param xp The value of xp, output
   * @param zp The value of zp, output
   */
  double generateComptonPoint(double &xp, double & zp);

  /**
   * Sample the BGF real-emission variables for the MEC path.
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
   * Return the compact azimuthal kernel for QCDC/Compton real emission.
   * @param xp \f$x_p\f$
   * @param x2 \f$x_2\f$
   * @param xperp \f$x_\perp\f$
   * @param norm Normalise to the large $l$ value of the ME
   */
  AzimuthalKernelCoefficients ComptonME(double xp, double x2, double xperp,
			                bool norm) const;
  
  /**
   * Return the compact azimuthal kernel for BGF real emission.
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
   *  Members for the POWHEG correction.
   */
  //@{
  /**
   *  Generate and cache a candidate QCDC/Compton POWHEG emission.
   */
  void generateCompton();

  /**
   *  Generate and cache a candidate BGF POWHEG emission.
   */
  void generateBGF();

  enum POWHEGEmissionComparisonMode {
    POWHEGEmissionComparisonModeDefault = 0,
    POWHEGEmissionComparisonModeRealOnly = 1
  };

  /**
   * Prepare the common Born-level state used by the POWHEG emission samplers.
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
   * Generate the configured RealOnly comparison-mode hardest emission.
   */
  RealEmissionProcessPtr generateComparisonModePOWHEGHardest(RealEmissionProcessPtr born);

  /**
   *  Parameters for the matrix-element correction.
   */
  //@{
  /**
   *  Enhancement factor for the ISR overestimate.
   */
  double initial_;

  /**
   *  Enhancement factor for the FSR overestimate.
   */
  double final_;

  /**
   *   Relative sampling fraction for QCDC/Compton versus BGF processes.
   */
  double procProb_;

  /**
   *  Analytic sampling integral for the QCDC/Compton channel.
   */
  double comptonInt_;

  /**
   *  Analytic sampling integral for the BGF channel.
   */
  double bgfInt_;
  //@}

  /**
   *  Parameters for the POWHEG correction.
   */
  //@{
  /**
   *  Overestimate weight for the QCDC/Compton channel.
   */
  double comptonWeight_;

  /**
   *  Overestimate weight for the BGF channel.
   */
  double BGFWeight_;

  /**
   *  Minimum transverse momentum for generated real emissions.
   */
  Energy pTmin_;
  //@}

  /**
   *  Cached state for the point being generated.
   */
  //@{
  /**
   *   \f$Q^2\f$
   */
  Energy2 q2_;

  /**
   * Born DIS angular variable ell = 2/y - 1.
   */
  double l_;

  /**
   *  Born Bjorken momentum fraction.
   */
  double xB_;

  /**
   *  Hadron beam particle for the current DIS point.
   */
  tcBeamPtr beam_;

  /**
   *  Incoming and outgoing Born-side parton data.
   */
  tcPDPtr partons_[2];

  /**
   *  Incoming and outgoing Born-side lepton data.
   */
  tcPDPtr leptons_[2];

  /**
   *  Unpolarized beam PDF object used for flux and emission ratios.
   */
  tcPDFPtr pdf_;
  /**
   *  Rotation to the Breit frame
   */
  LorentzRotation rot_;

  /**
   *  Born lepton momenta in the Breit-frame construction.
   */
  Lorentz5Momentum pl_[2];

  /**
   *  Born parton momenta in the Breit-frame construction.
   */
  Lorentz5Momentum pq_[2];

  /**
   *  Exchanged-current momentum.
   */
  Lorentz5Momentum q_;

  /**
   *  Cached QCDC/Compton candidate from POWHEG generation.
   */
  Energy pTCompton_;
  bool ComptonISFS_;
  vector<Lorentz5Momentum> ComptonMomenta_;

  /**
   *  Cached BGF candidate from POWHEG generation.
   */
  Energy pTBGF_;
  vector<Lorentz5Momentum> BGFMomenta_;
  //@}

  /**
   *  Born analyzing coefficient cached for legacy MEC kernels.
   */
  double acoeff_;

  /**
   *  Shower alphaS object used by the default emission kernels.
   */
  ShowerAlphaPtr alpha_;

  /**
   *  Gluon particle data object.
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
   * Convert the cached DIS window into tighter Born cos(theta) limits.
   */
  bool tightenBornCosThetaWithNativeDISWindow(const HwMEBase::TwoToTwoKinematicsSetup & setup,
                                              double xB,
                                              double & ctmin,
                                              double & ctmax) const;

  /**
   *  NLO radiative variables.
   */
  //@{
  /**
   *  The \f$x_p\f$ real integration variable for NLO contributions.
   */
  double xp_;
  //@}

  /**
   *  Hadron beam cached for the NLO and correlated-weight paths.
   */
  tcBeamPtr hadron_;

  /**
   * Selects a dynamic or fixed factorization scale.
   */
  unsigned int scaleOpt_;

  /**
   * The fixed factorization scale.
   */
  Energy muF_;

  /**
   *  Multiplicative factor for the dynamic scale option.
   */
  double scaleFact_;

  /**
   * Optional common lower boundary for factorization, renormalization, and
   * POWHEG-emission PDF scales.
   */
  Energy minimumScale_;

  /**
   *  Whether to generate the leading, positive NLO, or negative NLO stream.
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
   *  Sampling state for the mapped NLO xp variable.
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
