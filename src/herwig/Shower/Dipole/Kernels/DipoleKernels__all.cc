#line 1 "./DipoleSplittingKernel.cc"
// -*- C++ -*-
//
// DipoleSplittingKernel.cc is a part of Herwig - 
// A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleSplittingKernel class.
//

#include "DipoleSplittingKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/Shower/ShowerHandler.h"

using namespace Herwig;

DipoleSplittingKernel::DipoleSplittingKernel() 
  : HandlerBase(), theScreeningScale(0.0*GeV), 
    thePresamplingPoints(2000), theMaxtry(100000),
    theFreezeGrid(500000),
    theDetuning(1.0),
    theStrictLargeN(false), 
    theFactorizationScaleFactor(1.0),
    theRenormalizationScaleFactor(1.0),
    theRenormalizationScaleFreeze(1.*GeV), 
    theFactorizationScaleFreeze(1.*GeV),
    theVirtualitySplittingScale(false),
    theCMWScheme(0),
    presampling(false) {}

// initialize static variable out of line
double DipoleSplittingKernel::theMaxPDFRatio = 1000000.;

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleSplittingKernel::persistentOutput(PersistentOStream & os) const {
  os << theAlphaS << ounit(theScreeningScale,GeV)
     << theSplittingKinematics << thePDFRatio
     << thePresamplingPoints << theMaxtry << theFreezeGrid << theDetuning
     << theFlavour << theMCCheck << theStrictLargeN
     << theFactorizationScaleFactor
     << theRenormalizationScaleFactor
     << ounit(theRenormalizationScaleFreeze,GeV)
     << ounit(theFactorizationScaleFreeze,GeV)
     << theVirtualitySplittingScale<<theCMWScheme<<theUseThisKernel;
}

void DipoleSplittingKernel::persistentInput(PersistentIStream & is, int) {
  is >> theAlphaS >> iunit(theScreeningScale,GeV) 
     >> theSplittingKinematics >> thePDFRatio
     >> thePresamplingPoints >> theMaxtry >> theFreezeGrid >> theDetuning
     >> theFlavour >> theMCCheck >> theStrictLargeN
     >> theFactorizationScaleFactor
     >> theRenormalizationScaleFactor
     >> iunit(theRenormalizationScaleFreeze,GeV)
     >> iunit(theFactorizationScaleFreeze,GeV)
     >> theVirtualitySplittingScale>>theCMWScheme>>theUseThisKernel;
}

double DipoleSplittingKernel::alphaPDF(const DipoleSplittingInfo& split,
				       Energy optScale,
				       double rScaleFactor,
				       double fScaleFactor) const {

  Energy pt = optScale == ZERO ? split.lastPt() : optScale;

  Energy2 scale = ZERO;
  if ( !virtualitySplittingScale() ) {
    scale = sqr(pt) + sqr(theScreeningScale);
  } else {
    scale = sqr(splittingKinematics()->QFromPt(pt,split)) + sqr(theScreeningScale);
  }


 
  Energy2 fScale = sqr(theFactorizationScaleFactor*fScaleFactor)*scale;
  fScale = max( fScale , sqr(factorizationScaleFreeze()) );
 
  Energy2 rScale = sqr(theRenormalizationScaleFactor*rScaleFactor)*scale;
  rScale = max( rScale , sqr(renormalizationScaleFreeze()) );

  if(split.calcFixedExpansion()){
    fScale = max( sqr(split.fixedScale()) , sqr(factorizationScaleFreeze()) );
    rScale = max( sqr(split.fixedScale()) , sqr(renormalizationScaleFreeze()) );
  }

  double alphas = 1.0;
  double pdf = 1.0;

  // check if we are potentially reweighting and cache evaluations
  bool evaluatePDF = true;
  bool evaluateAlphaS = true;
  bool variations = 
    !ShowerHandler::currentHandler()->showerVariations().empty() &&
    !presampling;
  if ( variations ) {
    
    map<double,double>::const_iterator pit = thePDFCache.find(fScaleFactor);
    evaluatePDF = (pit == thePDFCache.end());
    if ( !evaluatePDF ) {
      pdf = pit->second;
    }
    map<double,double>::const_iterator ait = theAlphaSCache.find(rScaleFactor);
    evaluateAlphaS = (ait == theAlphaSCache.end());
    if ( !evaluateAlphaS ) {
      alphas = ait->second;
    }
  }

  if ( evaluateAlphaS ){
    if (theCMWScheme==0||split.calcFixedExpansion()) {
      alphas = alphaS()->value(rScale);
    }else if(theCMWScheme==1){
      
      alphas = alphaS()->value(rScale);
      alphas *=1.+(3.*(67./18.-1./6.*sqr(Constants::pi))
                   -5./9.*alphaS()->Nf(rScale))*
               alphas/2./Constants::pi;
      
    }else if(theCMWScheme==2){
      double kg=exp(-(67.-3.*sqr(Constants::pi)-10/3*alphaS()->Nf(rScale))
                    /(33.-2.*alphaS()->Nf(rScale)));
      Energy2 cmwscale2=max(kg*rScale, sqr(renormalizationScaleFreeze()) );
      alphas = alphaS()->value(cmwscale2);
      
    }else{
      throw Exception()
      << "This CMW-Scheme is not implemented."
      << Exception::abortnow;
    
    }
  }
  if ( evaluatePDF ) {
    if ( split.index().initialStateEmitter() ) {
      assert(pdfRatio());
      pdf *= 
	split.lastEmitterZ() * 
	(*pdfRatio())(split.index().emitterPDF(), fScale,
		      split.index().emitterData(),split.emitterData(),
		      split.emitterX(),split.lastEmitterZ());
    }

    if ( split.index().initialStateSpectator() ) {
      assert(pdfRatio());
      pdf *= 
	split.lastSpectatorZ() * 
	(*pdfRatio())(split.index().spectatorPDF(), fScale,
		      split.index().spectatorData(),split.spectatorData(),
		      split.spectatorX(),split.lastSpectatorZ());
    }
  }

  if ( evaluatePDF && variations ) {
    thePDFCache[fScaleFactor] = min(pdf,theMaxPDFRatio);
  }

  if ( evaluateAlphaS && variations ) {
    theAlphaSCache[rScaleFactor] = alphas;
  }

  double ret = min(pdf,theMaxPDFRatio)*
               (split.calcFixedExpansion()?
                1.:(alphas / (2.*Constants::pi)));

  if ( ret < 0. )
    ret = 0.;

  return ret;

}

void DipoleSplittingKernel::accept(const DipoleSplittingInfo& split,
				   double, double,
				   map<string,double>& weights) const {
  if ( ShowerHandler::currentHandler()->showerVariations().empty() )
    return;
  double reference = alphaPDF(split);
  assert(reference > 0.);
  for ( map<string,ShowerVariation>::const_iterator var =
	  ShowerHandler::currentHandler()->showerVariations().begin();
	var != ShowerHandler::currentHandler()->showerVariations().end(); ++var ) {
    if ( ( ShowerHandler::currentHandler()->firstInteraction() 
	   && var->second.firstInteraction ) ||
	 ( !ShowerHandler::currentHandler()->firstInteraction() 
           && var->second.secondaryInteractions ) ) {
      double varied = alphaPDF(split,ZERO,
			       var->second.renormalizationScaleFactor,
			       var->second.factorizationScaleFactor);
      if ( varied != reference ) {
	map<string,double>::iterator wi = weights.find(var->first);
	if ( wi != weights.end() )
	  wi->second *= varied/reference;
	else
	  weights[var->first] = varied/reference;
      }
    }
  }
}

void DipoleSplittingKernel::veto(const DipoleSplittingInfo& split,
				 double p, double r,
				 map<string,double>& weights) const {
  if ( ShowerHandler::currentHandler()->showerVariations().empty() )
    return;
  double reference = alphaPDF(split);
  // this is dangerous, but we have no other choice currently -- need to
  // carefully check for the effects; the assumption is that if the central
  // one ius zero, then so will be the variations.
  if ( reference == 0.0 )
    return;
  for ( map<string,ShowerVariation>::const_iterator var =
	  ShowerHandler::currentHandler()->showerVariations().begin();
	var != ShowerHandler::currentHandler()->showerVariations().end(); ++var ) {
    if ( ( ShowerHandler::currentHandler()->firstInteraction() 
           && var->second.firstInteraction ) ||
	 ( !ShowerHandler::currentHandler()->firstInteraction() 
           && var->second.secondaryInteractions ) ) {
      double varied = alphaPDF(split,ZERO,
			       var->second.renormalizationScaleFactor,
			       var->second.factorizationScaleFactor);
      if ( varied != reference ) {
	map<string,double>::iterator wi = weights.find(var->first);
	if ( wi != weights.end() )
	  wi->second *= (r - varied*p/reference) / (r-p);
	else
	  weights[var->first] = (r - varied*p/reference) / (r-p);
      }
    }
  }
}

AbstractClassDescription<DipoleSplittingKernel> 
DipoleSplittingKernel::initDipoleSplittingKernel;
// Definition of the static class description member.

void DipoleSplittingKernel::Init() {

  static ClassDocumentation<DipoleSplittingKernel> documentation
    ("DipoleSplittingKernel is the base class for all kernels "
     "used within the dipole shower.");

  static Reference<DipoleSplittingKernel,AlphaSBase> interfaceAlphaS
    ("AlphaS",
     "The strong coupling to be used by this splitting kernel.",
     &DipoleSplittingKernel::theAlphaS, false, false, true, true, false);


  static Parameter<DipoleSplittingKernel,Energy> interfaceScreeningScale
    ("ScreeningScale",
     "A colour screening scale",
     &DipoleSplittingKernel::theScreeningScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


  static Reference<DipoleSplittingKernel,DipoleSplittingKinematics> 
     interfaceSplittingKinematics
    ("SplittingKinematics",
     "The splitting kinematics to be used by this splitting kernel.",
     &DipoleSplittingKernel::theSplittingKinematics, false, false, true, false, false);


  static Reference<DipoleSplittingKernel,PDFRatio> interfacePDFRatio
    ("PDFRatio",
     "Set the optional PDF ratio object to evaluate this kernel",
     &DipoleSplittingKernel::thePDFRatio, false, false, true, true, false);

  static Parameter<DipoleSplittingKernel,unsigned long> interfacePresamplingPoints
    ("PresamplingPoints",
     "The number of points used to presample this kernel.",
     &DipoleSplittingKernel::thePresamplingPoints, 2000, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<DipoleSplittingKernel,unsigned long> interfaceMaxtry
    ("Maxtry",
     "The maximum number of attempts to generate a splitting.",
     &DipoleSplittingKernel::theMaxtry, 10000, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<DipoleSplittingKernel,unsigned long> interfaceFreezeGrid
    ("FreezeGrid",
     "",
     &DipoleSplittingKernel::theFreezeGrid, 500000, 1, 0,
     false, false, Interface::lowerlim);

  static Reference<DipoleSplittingKernel,ParticleData> interfaceFlavour
    ("Flavour",
     "Set the flavour to be produced if ambiguous.",
     &DipoleSplittingKernel::theFlavour, false, false, true, true, false);

  static Reference<DipoleSplittingKernel,DipoleMCCheck> interfaceMCCheck
    ("MCCheck",
     "[debug option] MCCheck",
     &DipoleSplittingKernel::theMCCheck, false, false, true, true, false);

  interfaceMCCheck.rank(-1);

  static Switch<DipoleSplittingKernel,bool> interfaceStrictLargeN
    ("StrictLargeN",
     "Work in a strict large-N limit.",
     &DipoleSplittingKernel::theStrictLargeN, false, false, false);
  static SwitchOption interfaceStrictLargeNYes
    (interfaceStrictLargeN,
     "Yes",
     "Replace C_F -> C_A/2 where present",
     true);
  static SwitchOption interfaceStrictLargeNNo
    (interfaceStrictLargeN,
     "No",
     "Keep C_F=4/3",
     false);

  interfaceStrictLargeN.rank(-2);

  static Switch<DipoleSplittingKernel,unsigned int> interfaceCMWScheme
    ("CMWScheme",
     "Use the CMW Scheme related Kg expression to the splitting",
    &DipoleSplittingKernel::theCMWScheme, 0, false, false);
  static SwitchOption interfaceCMWSchemeNo
    (interfaceCMWScheme,"No","No CMW-Scheme", 0);
  static SwitchOption interfaceCMWSchemeLinear
  (interfaceCMWScheme,"Linear",
   "Linear CMW multiplication: alpha_s(q) -> alpha_s(q)(1+K_g*alpha_s(q)/2pi )",1);
  static SwitchOption interfaceCMWSchemeFactor
  (interfaceCMWScheme,"Factor",
   "Use factor in alpha_s argument: alpha_s(q) -> alpha_s(k_g*q) with  kfac=exp(-(67-3pi^2-10/3*Nf)/(33-2Nf)) ",2);

  static Parameter<DipoleSplittingKernel,double> interfaceFactorizationScaleFactor
    ("FactorizationScaleFactor",
     "The factorization scale factor.",
     &DipoleSplittingKernel::theFactorizationScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  interfaceFactorizationScaleFactor.rank(-2);

  static Parameter<DipoleSplittingKernel,double> interfaceRenormalizationScaleFactor
    ("RenormalizationScaleFactor",
     "The renormalization scale factor.",
     &DipoleSplittingKernel::theRenormalizationScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  interfaceRenormalizationScaleFactor.rank(-2);

  static Parameter<DipoleSplittingKernel,Energy> interfaceRenormalizationScaleFreeze
    ("RenormalizationScaleFreeze",
     "The freezing scale for the renormalization scale.",
     &DipoleSplittingKernel::theRenormalizationScaleFreeze, 
      GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<DipoleSplittingKernel,Energy> interfaceFactorizationScaleFreeze
    ("FactorizationScaleFreeze",
     "The freezing scale for the factorization scale.",
     &DipoleSplittingKernel::theFactorizationScaleFreeze, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Switch<DipoleSplittingKernel,bool> interfaceVirtualitySplittingScale
    ("VirtualitySplittingScale",
     "Use the virtuality as the splitting scale.",
     &DipoleSplittingKernel::theVirtualitySplittingScale, false, false, false);
  static SwitchOption interfaceVirtualitySplittingScaleYes
    (interfaceVirtualitySplittingScale,
     "Yes",
     "Use vrituality.",
     true);
  static SwitchOption interfaceVirtualitySplittingScaleNo
    (interfaceVirtualitySplittingScale,
     "No",
     "Use transverse momentum.",
     false);

  static Parameter<DipoleSplittingKernel,double> interfaceDetuning
  ("Detuning",
   "A value to detune the overestimate kernel.",
   &DipoleSplittingKernel::theDetuning, 1.0, 1.0, 0,
   false, false, Interface::lowerlim);
  
  
  static Switch<DipoleSplittingKernel,bool> interfaceUseThisKernel
  ("UseKernel",
   "Turn On and of the Kernel.",
   &DipoleSplittingKernel::theUseThisKernel, true, false, false);
  static SwitchOption interfaceUseThisKernelYes
  (interfaceUseThisKernel,
   "Yes",
   "Use this Kernel.",
   true);
  static SwitchOption interfaceUseThisKernelNo
  (interfaceUseThisKernel,
   "No",
   "Dont use this Kernel.",
   false);
  

}

#line 1 "./ColourMatrixElementCorrection.cc"
// -*- C++ -*-
//
// ColourMatrixElementCorrection.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourMatrixElementCorrection class.
//

#include "ColourMatrixElementCorrection.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/Shower/Dipole/DipoleShowerHandler.h"

using namespace Herwig;

ColourMatrixElementCorrection::ColourMatrixElementCorrection():lambda(1.0),negCMECScaling(1.0) {}

IBPtr ColourMatrixElementCorrection::clone() const {
  return new_ptr(*this);
}

IBPtr ColourMatrixElementCorrection::fullclone() const {
  return new_ptr(*this);
}

double ColourMatrixElementCorrection::cmec(const DipoleSplittingInfo& dsplit) const {
  // Do not calculate CMECs if the subleading Nc shower has stopped 
  if ( !(currentHandler()->continueSubleadingNc()) ) {
    return 1.;
  }
  const PPtr em = dsplit.emitter();
  const PPtr sp = dsplit.spectator();
  const tcPDPtr emis = dsplit.emissionData();
  // Get the dictionary from particle pointers to herwig particle numbers
  const map<PPtr,size_t>& theDictionary = currentHandler()->particleIndices();
  const cPDVector& theParticles = currentHandler()->particlesAfter();
  // Use the dictionary to find
  // - emitter index
  // - spectator index
  // - emission ParticleID
  const std::tuple<size_t,size_t,long> ikemission = std::make_tuple(theDictionary.at(em),
								    theDictionary.at(sp),
								    emis->id());
  double factor = currentHandler()->densityOperator().colourMatrixElementCorrection(ikemission,theParticles);

  return factor > 0.0 ? factor : negCMECScaling*factor;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ColourMatrixElementCorrection::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void ColourMatrixElementCorrection::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ColourMatrixElementCorrection,Herwig::DipoleSplittingReweight>
  describeHerwigColourMatrixElementCorrection("Herwig::ColourMatrixElementCorrection", 
					      "HwDipoleShower.so");

void ColourMatrixElementCorrection::Init() {

  static ClassDocumentation<ColourMatrixElementCorrection> documentation
    ("There is no documentation for the ColourMatrixElementCorrection class");

}

#line 1 "./FFMqx2qgxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMqx2qgxDipoleKernel class.
//

#include "FFMqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFMqx2qgxDipoleKernel::FFMqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr FFMqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFMqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFMqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 7  &&
    // 2012-05-01
    abs(ind.emitterData()->id()) == abs(flavour()->id()) &&
    !( ind.emitterData()->mass() == ZERO &&
       ind.spectatorData()->mass() == ZERO ) &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator() &&
    !ind.incomingDecayEmitter() && !ind.incomingDecaySpectator();
}

bool FFMqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {
  
  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    abs(sk.emitter(b)->id()) < 7 &&
    abs(sk.emitter(b)->mass()) == abs(emitter(a)->mass()) &&
    abs(sk.spectator(b)->mass()) == abs(spectator(a)->mass());

}

// 2012-05-01
tcPDPtr FFMqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  assert(flavour());
  assert(abs(flavour()->id())<7);
  return ind.emitterData()->id() > 0 ?
    (tcPDPtr) flavour() : (tcPDPtr) flavour()->CC();
}

tcPDPtr FFMqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFMqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFMqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);

  // Sudakov parameterisation variables,
  // needed to calculate y.
  double z = split.lastZ();
  Energy pt = split.lastPt();

  // Construct mass squared variables
  // Note for q->qg can use the emitterMass
  // (i.e. mass of emitter before splitting = mass of emitter after)
  Energy2 Qijk = sqr(split.scale());
  Energy2 mi2 = sqr(split.emitterMass());
  Energy2 mk2 = sqr(split.spectatorMass());
  Energy2 sbar = Qijk - mi2 - mk2;

  // Calculate y
  double y = (sqr(pt) + sqr(1.-z)*mi2) / sbar / z / (1.-z);

  // zi, used in dipole splitting kernel
  double zi = split.lastSplittingParameters()[0];

  double vijk = sqrt( sqr(2.*mk2 + sbar*(1.-y)) - 4.*mk2*Qijk ) / sbar / (1.-y);
  double vtilde = sqrt( sqr(Qijk) + sqr(mi2) + sqr(mk2)
                        - 2.*(mi2*Qijk + mk2*Qijk + mi2*mk2) ) / sbar;
  
  ret *= (!strictLargeN() ? 4./3. : 3./2.)*
    ( 2./(1.-zi*(1.-y)) - vtilde/vijk*( 1. + zi + 2.*mi2/sbar/y ) );
  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FFMqx2qgxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {
  
  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr FFMqx2qgxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {

  // Need variables for the AP kernels
  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();
  Energy mi = dInfo.emitterMass();

  // Altarelli-Parisi spin-indexed kernels:
  Energy den = sqrt(sqr(mi)*sqr(1.-z) + sqr(pt));
  double v_AP_ppp = pt / den / sqrt(1.-z);
  double v_AP_ppm = - z * v_AP_ppp ;
  double v_AP_pmp = mi*(1.-z)*sqrt(1.-z) / den ;

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = v_AP_pmp;

  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -1 or -1/2, 1=+1/2, 2 = +1
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,1,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(1,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm;
  (*kernelPhiDep)(1,0,2) = v_AP_pmp;
  (*kernelPhiDep)(0,1,2) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FFMqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFMqx2qgxDipoleKernel> FFMqx2qgxDipoleKernel::initFFMqx2qgxDipoleKernel;
// Definition of the static class description member.

void FFMqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<FFMqx2qgxDipoleKernel> documentation
    ("FFMqx2qgxDipoleKernel");

}

#line 1 "./FFMgx2ggxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMgx2ggxDipoleKernel class.
//

#include "FFMgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFMgx2ggxDipoleKernel::FFMgx2ggxDipoleKernel() 
  : DipoleSplittingKernel(){}

IBPtr FFMgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFMgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFMgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() != ZERO &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator() &&
    !ind.incomingDecayEmitter() && !ind.incomingDecaySpectator();

}

bool FFMgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {
    
  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    abs(spectator(a)->mass()) == abs(sk.spectator(b)->mass());

}

tcPDPtr FFMgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFMgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFMgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFMgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);

  // Sudakov parameterisation variables,
  // needed to calculate y.
  double z = split.lastZ();
  Energy pt = split.lastPt();
  
  // Construct mass squared variables
  Energy2 Qijk = sqr(split.scale());
  Energy2 mk2 = sqr(split.spectatorMass());
  Energy2 sbar = Qijk - mk2;

  // Calculate y
  double y = sqr(pt) / sbar / z / (1.-z);

  double vijk = sqrt( sqr(2.*mk2 + sbar*(1.-y)) - 4.*mk2*Qijk ) / sbar / (1.-y);
  double viji = 1.;

  // zi, used in dipole splitting kernel
  double zi = split.lastSplittingParameters()[0];

  double zip = 0.5*(1.+viji*vijk);
  double zim = 0.5*(1.-viji*vijk);
  
  // how to choose kappa?
  double kappa = 0.;

  double S1=1./(1.-zi*(1.-y));
  double S2=1./(1.-(1.-zi)*(1.-y));
  double NS=(zi*(1.-zi)-(1.-kappa)*zip*zim-2.)/vijk;
  
  if( theAsymmetryOption == 0 ){
    ret *= 3.*( S1 + 0.5 * NS);
  }else if ( theAsymmetryOption == 1 ){
    ret *= 3.*zi*( S1 +S2 + NS );
  }else{
    ret *= 3.*0.5*( S1 + S2 + NS );
  }

  return ret > 0. ? ret : 0.;  
}

vector< pair<int, Complex> >
FFMgx2ggxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  // Need variables for the AP kernels
  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_pmp = (1.-z)*sqrt( (1.-z)/z );

  //double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp)) - 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);

  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp))/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );

  return distPhiDep;
}

DecayMEPtr FFMgx2ggxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  // Need variables for the AP kernels
  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_pmp = (1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(2,2,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,2,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,2) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,2) = 0;
  (*kernelPhiDep)(2,0,0) = 0;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMgx2ggxDipoleKernel::persistentOutput(PersistentOStream & os) const {
  os << theAsymmetryOption;
}

void FFMgx2ggxDipoleKernel::persistentInput(PersistentIStream & is, int) {
  is >> theAsymmetryOption;
}

ClassDescription<FFMgx2ggxDipoleKernel> FFMgx2ggxDipoleKernel::initFFMgx2ggxDipoleKernel;
// Definition of the static class description member.

void FFMgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<FFMgx2ggxDipoleKernel> documentation
    ("FFMgx2ggxDipoleKernel");
  
  static Parameter<FFMgx2ggxDipoleKernel,int> interfacetheAsymmetryOption
  ("AsymmetryOption",
   "The asymmetry option for final state gluon spliitings.",
   &FFMgx2ggxDipoleKernel::theAsymmetryOption, 1, 0, 0,
   false, false, Interface::lowerlim);

}
#line 1 "./FFMgx2qqxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMgx2qqxDipoleKernel class.
//

#include "FFMgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFMgx2qqxDipoleKernel::FFMgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr FFMgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFMgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFMgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    !( ind.spectatorData()->mass() == ZERO &&
       flavour()->mass() == ZERO ) &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator() &&
    !ind.incomingDecayEmitter() && !ind.incomingDecaySpectator();
}

bool FFMgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {
  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    emitter(a)->id() == sk.emitter(b)->id() &&
    abs(sk.spectator(b)->mass()) == abs(spectator(a)->mass());
}

tcPDPtr FFMgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6);
  return flavour();
}

tcPDPtr FFMgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6);
  return flavour()->CC();
}

tcPDPtr FFMgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFMgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);

  // Sudakov parameterisation variables,
  // needed to calculate y.
  double z = split.lastZ();
  Energy pt = split.lastPt();

  // Construct mass squared variables
  Energy2 Qijk = sqr(split.scale());
  Energy2 mi2 = sqr(split.emitterData()->mass());
  Energy2 mj2 = mi2;
  Energy2 mk2 = sqr(split.spectatorMass());
  Energy2 sbar = Qijk - mi2 - mj2 - mk2;

  // Calculate y
  double y = (sqr(pt) + sqr(1.-z)*mi2 + sqr(z)*mj2) / sbar / z / (1.-z);
  
  // zi, used in dipole splitting kernel
  double zi = split.lastSplittingParameters()[0];
  
  double vijk = sqrt( sqr(2.*mk2 + sbar*(1.-y)) - 4.*mk2*Qijk ) / sbar / (1.-y);
  double viji = sqrt( sqr(sbar*y) - 4.*sqr(mi2) ) / (sbar*y + 2.*mi2);

  double zip = 0.5*(1.+viji*vijk);
  double zim = 0.5*(1.-viji*vijk);

  // how to choose kappa??
  double kappa = 0.;

  ret *= 0.25 / vijk *
    ( 1. - 2.*( zi*(1.-zi) - (1.-kappa)*zip*zim 
                - kappa*mi2 / ( 2.*mi2 + (Qijk - 2.*mi2 - mk2)*y) ) );
    
  return ret > 0. ? ret : 0.;
  
}

vector< pair<int, Complex> >
FFMgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  // Need variables for the AP kernels
  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();
  Energy2 mi2 = sqr(dInfo.emitterData()->mass());

  // Altarelli-Parisi spin-indexed kernels:
  double ratio = mi2 / ( mi2 + sqr(pt) );
  double root = sqrt(1.-ratio);
  double v_AP_ppp = sqrt(ratio);
  double v_AP_ppm = z*root;
  double v_AP_pmp = -(1.-z)*root;

  //double v_AP_mmm = v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp)) + 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);
  
  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp) )/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );

  return distPhiDep;
}

DecayMEPtr FFMgx2qqxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {
  
  // Need variables for the AP kernels
  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();
  Energy2 mi2 = sqr(dInfo.emitterData()->mass());

  // Altarelli-Parisi spin-indexed kernels:
  double ratio = mi2 / ( mi2 + sqr(pt) );
  double root = sqrt(1.-ratio);
  double v_AP_ppp = sqrt(ratio);
  double v_AP_ppm = z*root;
  double v_AP_pmp = -(1.-z)*root;

  double v_AP_mmm = v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm;
  (*kernelPhiDep)(2,1,1) = v_AP_ppp;
  (*kernelPhiDep)(0,0,1) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,1,1) = 0.;
  (*kernelPhiDep)(2,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FFMgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFMgx2qqxDipoleKernel> 
FFMgx2qqxDipoleKernel::initFFMgx2qqxDipoleKernel;
// Definition of the static class description member.

void FFMgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FFMgx2qqxDipoleKernel> documentation
    ("FFMgx2qqxDipoleKernel");

}

#line 1 "./FFqx2qgxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFqx2qgxDipoleKernel class.
//

#include "FFqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFqx2qgxDipoleKernel::FFqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr FFqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}
#ifndef NDEBUG
bool FFqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
#else
bool FFqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& ,
#endif
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    abs(sk.emitter(b)->id()) < 6 &&
    sk.emitter(b)->mass() == ZERO;
       

}

tcPDPtr FFqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  return ind.emitterData();
}

tcPDPtr FFqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double y = sqr(split.lastPt() / split.scale()) / (z*(1.-z));

  ret *= (!strictLargeN() ? 4./3. : 3./2.)*( 2./(1.-z*(1.-y)) - (1.+z) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FFqx2qgxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {
  
  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr FFqx2qgxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {

  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt( 1./(1.-z) );
  double v_AP_ppm = -z/sqrt(1.-z);

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -1 or -1/2, 1=+1/2, 2 = +1
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,1,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(1,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = 0.;
  (*kernelPhiDep)(1,0,2) = 0.;
  (*kernelPhiDep)(0,1,2) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FFqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFqx2qgxDipoleKernel> FFqx2qgxDipoleKernel::initFFqx2qgxDipoleKernel;
// Definition of the static class description member.

void FFqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<FFqx2qgxDipoleKernel> documentation
    ("FFqx2qgxDipoleKernel");

}

#line 1 "./FFgx2ggxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFgx2ggxDipoleKernel class.
//

#include "FFgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFgx2ggxDipoleKernel::FFgx2ggxDipoleKernel() 
  : DipoleSplittingKernel(){}

IBPtr FFgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

#ifndef NDEBUG
bool FFgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
#else
bool FFgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& ,
#endif
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g;
       

}

tcPDPtr FFgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double y = sqr(split.lastPt() / split.scale()) / (z*(1.-z));

  double S1=1./(1.-z*(1.-y));
  double S2=1./(1.-(1.-z)*(1.-y));
  double NS=(-2 + z*(1.-z));
  
  if( theAsymmetryOption == 0 ){
    ret *= 3.*( S1 + 0.5 * NS);
  }else if ( theAsymmetryOption == 1 ){
    ret *= 3.*z*( S1 +S2 + NS );
  }else{
    ret *= 3.*0.5*( S1 + S2 + NS );
  }
  
  return ret > 0. ? ret : 0.;
}

vector< pair<int, Complex> >
FFgx2ggxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_pmp = (1.-z)*sqrt( (1.-z)/z );

  //double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp)) - 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);

  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp))/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );

  return distPhiDep;
}

DecayMEPtr FFgx2ggxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_pmp = (1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(2,2,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,2,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,2) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,2) = 0;
  (*kernelPhiDep)(2,0,0) = 0;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFgx2ggxDipoleKernel::persistentOutput(PersistentOStream & os) const {
  os << theAsymmetryOption;
}

void FFgx2ggxDipoleKernel::persistentInput(PersistentIStream & is, int) {
  is >> theAsymmetryOption;
}

ClassDescription<FFgx2ggxDipoleKernel> FFgx2ggxDipoleKernel::initFFgx2ggxDipoleKernel;
// Definition of the static class description member.

void FFgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<FFgx2ggxDipoleKernel> documentation
    ("FFgx2ggxDipoleKernel");

  static Parameter<FFgx2ggxDipoleKernel,int> interfacetheAsymmetryOption
    ("AsymmetryOption",
     "The asymmetry option for final state gluon spliitings.",
     &FFgx2ggxDipoleKernel::theAsymmetryOption, 1, 0, 0,
     false, false, Interface::lowerlim);
}
#line 1 "./FFgx2qqxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFgx2qqxDipoleKernel class.
//

#include "FFgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFgx2qqxDipoleKernel::FFgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr FFgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() == ZERO &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

#ifndef NDEBUG
bool FFgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
#else
bool FFgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& ,
#endif
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    sk.emitter(b)->mass() == ZERO;
       
}


tcPDPtr FFgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr FFgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour()->CC();
}

tcPDPtr FFgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();

  ret *= .25 * ( 1. - 2.*z*(1.-z) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FFgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppm = z;
  double v_AP_pmp = -(1.-z);

  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppm) + sqr(v_AP_pmp)) + 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);
  
  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*( sqr(v_AP_ppm) + sqr(v_AP_pmp) )/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );

  return distPhiDep;
}

DecayMEPtr FFgx2qqxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {
  
  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppm = z;
  double v_AP_pmp = -(1.-z);

  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());
  
  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = 0.;
  (*kernelPhiDep)(2,1,1) = 0.;
  (*kernelPhiDep)(0,0,1) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,1,1) = 0.;
  (*kernelPhiDep)(2,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FFgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFgx2qqxDipoleKernel> FFgx2qqxDipoleKernel::initFFgx2qqxDipoleKernel;
// Definition of the static class description member.

void FFgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FFgx2qqxDipoleKernel> documentation
    ("FFgx2qqxDipoleKernel");

}

#line 1 "./FIqx2qgxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIqx2qgxDipoleKernel class.
//

#include "FIqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIqx2qgxDipoleKernel::FIqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr FIqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6 &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    abs(sk.emitter(b)->id()) < 6 &&
    sk.emitter(b)->mass() == ZERO &&
    a.spectatorPDF() == b.spectatorPDF();
       

}


tcPDPtr FIqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  return ind.emitterData();
}

tcPDPtr FIqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double x = 1. / ( 1. + sqr(split.lastPt()/split.scale()) / (z*(1.-z)) );

  ret *= (!strictLargeN() ? 4./3. : 3./2.) * ( 2./(1.-z+(1.-x)) -(1.+z) + (1.-x)*(1.+3.*x*z) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FIqx2qgxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {

  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr FIqx2qgxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {

  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt( 1./(1.-z) );
  double v_AP_ppm = -z/sqrt(1.-z);

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -1 or -1/2, 1=+1/2, 2 = +1
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,1,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(1,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = 0.;
  (*kernelPhiDep)(1,0,2) = 0.;
  (*kernelPhiDep)(0,1,2) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIqx2qgxDipoleKernel> FIqx2qgxDipoleKernel::initFIqx2qgxDipoleKernel;
// Definition of the static class description member.

void FIqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<FIqx2qgxDipoleKernel> documentation
    ("FIqx2qgxDipoleKernel");

}

#line 1 "./FIgx2ggxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIgx2ggxDipoleKernel class.
//

#include "FIgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIgx2ggxDipoleKernel::FIgx2ggxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr FIgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    a.spectatorPDF() == b.spectatorPDF();

}

tcPDPtr FIgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double x = 1. / ( 1. + sqr(split.lastPt()/split.scale()) / (z*(1.-z)) );

  double S1=1./(1.-z+(1.-x));
  double S2=1./(1.-(1.-z)+(1.-x));
  double NS=(-2 + z*(1.-z)+(1.-x)*(1.+x*z*(1.-z)));
  
  if( theAsymmetryOption == 0 ){
    ret *= 3.*( S1 + 0.5 * NS);
  }else if ( theAsymmetryOption == 1 ){
    ret *= 3.*z*( S1 +S2 + NS );
  }else{
    ret *= 3.*0.5*( S1 + S2 + NS );
  }
  
  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FIgx2ggxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_pmp = (1.-z)*sqrt( (1.-z)/z );

  //double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;
  
  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp)) - 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);
  
  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp))/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );
  
  return distPhiDep;
}

DecayMEPtr FIgx2ggxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_pmp = (1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(2,2,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,2,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,2) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,2) = 0;
  (*kernelPhiDep)(2,0,0) = 0;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIgx2ggxDipoleKernel::persistentOutput(PersistentOStream & os) const {

  os << theAsymmetryOption;
}

void FIgx2ggxDipoleKernel::persistentInput(PersistentIStream & is, int) {
  is >> theAsymmetryOption;
}

ClassDescription<FIgx2ggxDipoleKernel> FIgx2ggxDipoleKernel::initFIgx2ggxDipoleKernel;
// Definition of the static class description member.

void FIgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<FIgx2ggxDipoleKernel> documentation
    ("FIgx2ggxDipoleKernel");

  static Parameter<FIgx2ggxDipoleKernel,int> interfacetheAsymmetryOption
  ("AsymmetryOption",
   "The asymmetry option for final state gluon spliitings.",
   &FIgx2ggxDipoleKernel::theAsymmetryOption, 1, 0, 0,
   false, false, Interface::lowerlim);
}

#line 1 "./FIgx2qqxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIgx2qqxDipoleKernel class.
//

#include "FIgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIgx2qqxDipoleKernel::FIgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr FIgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() == ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    // sk.emitter(b)->mass() == ZERO &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr FIgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr FIgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour()->CC();
}

tcPDPtr FIgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();

  ret *= .25 * (1.-2.*z*(1.-z));

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FIgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppm = z;
  double v_AP_pmp = -(1.-z);

  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppm) + sqr(v_AP_pmp)) + 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);
  
  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*( sqr(v_AP_ppm) + sqr(v_AP_pmp) )/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );

  return distPhiDep;
}

DecayMEPtr FIgx2qqxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {
  
  double z = dInfo.lastZ();
    
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppm = z;
  double v_AP_pmp = -(1.-z);

  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = 0.;
  (*kernelPhiDep)(2,1,1) = 0.;
  (*kernelPhiDep)(0,0,1) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,1,1) = 0.;
  (*kernelPhiDep)(2,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIgx2qqxDipoleKernel> FIgx2qqxDipoleKernel::initFIgx2qqxDipoleKernel;
// Definition of the static class description member.

void FIgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FIgx2qqxDipoleKernel> documentation
    ("FIgx2qqxDipoleKernel");

}

#line 1 "./IFqx2qgxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFqx2qgxDipoleKernel class.
//

#include "IFqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFqx2qgxDipoleKernel::IFqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IFqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    emitter(a) == sk.emitter(b) &&
    emission(a) == sk.emission(b) &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  return ind.emitterData();
}

tcPDPtr IFqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());

  double rho = 1. - 4.*ratio*z*(1.-z)/sqr(1.-z+ratio);
  double x = 0.5*((1.-z+ratio)/ratio)*(1.-sqrt(rho));
  double u = 0.5*((1.-z+ratio)/(1.-z))*(1.-sqrt(rho));

  ret *= (!strictLargeN() ? 4./3. : 3./2.) * ( 2./(1.-x+u) - (1.+x) + u*(1.+3.*x*(1.-u) ) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IFqx2qgxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {

  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr IFqx2qgxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt( 1./(1.-z) );
  double v_AP_ppm = -z/sqrt(1.-z);

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half, PDT::Spin1Half, PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,1,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(1,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = 0.;
  (*kernelPhiDep)(1,0,2) = 0.;
  (*kernelPhiDep)(0,1,2) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;
  
  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFqx2qgxDipoleKernel> IFqx2qgxDipoleKernel::initIFqx2qgxDipoleKernel;
// Definition of the static class description member.

void IFqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<IFqx2qgxDipoleKernel> documentation
    ("IFqx2qgxDipoleKernel");

}

#line 1 "./IFqx2gqxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFqx2gqxDipoleKernel class.
//

#include "IFqx2gqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFqx2gqxDipoleKernel::IFqx2gqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IFqx2gqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFqx2gqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFqx2gqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFqx2gqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    a.emitterData() == b.emitterData() &&
    emitter(a) == sk.emitter(b) &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFqx2gqxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFqx2gqxDipoleKernel::emission(const DipoleIndex& ind) const {
  return ind.emitterData()->CC();
}

tcPDPtr IFqx2gqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFqx2gqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());

  double rho = 1. - 4.*ratio*z*(1.-z)/sqr(1.-z+ratio);
  double x = 0.5*((1.-z+ratio)/ratio)*(1.-sqrt(rho));

  ret *= .5 * ( 1.-2.*x*(1.-x)  );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IFqx2gqxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {

  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr IFqx2gqxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppm = z;
  double v_AP_mpm = (1.-z);
  
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_pmp = -v_AP_mpm;
    
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1, PDT::Spin1Half, PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = 0.;
  (*kernelPhiDep)(2,1,1) = 0.;
  (*kernelPhiDep)(0,0,1) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,1,1) = 0.;
  (*kernelPhiDep)(2,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFqx2gqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFqx2gqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFqx2gqxDipoleKernel> IFqx2gqxDipoleKernel::initIFqx2gqxDipoleKernel;
// Definition of the static class description member.

void IFqx2gqxDipoleKernel::Init() {

  static ClassDocumentation<IFqx2gqxDipoleKernel> documentation
    ("IFqx2gqxDipoleKernel");

}

#line 1 "./IFgx2ggxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFgx2ggxDipoleKernel class.
//

#include "IFgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFgx2ggxDipoleKernel::IFgx2ggxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IFgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());

  double rho = 1. - 4.*ratio*z*(1.-z)/sqr(1.-z+ratio);
  double x = 0.5*((1.-z+ratio)/ratio)*(1.-sqrt(rho));
  double u = 0.5*((1.-z+ratio)/(1.-z))*(1.-sqrt(rho));

  ret *= 3. * ( 1./(1.-x+u) + (1.-x)/x - 1. + x*(1.-x) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IFgx2ggxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_mpm = -(1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  //double v_AP_mmp = -v_AP_ppm;
  double v_AP_pmp = -v_AP_mpm;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_mpm)) - 2.*abs(rho(0,2))*(v_AP_ppp*v_AP_pmp + v_AP_mmm*v_AP_mpm);

  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_mpm))/max ) );
  distPhiDep.push_back( make_pair(2, rho(0,2)*(v_AP_mmm*v_AP_mpm + v_AP_pmp*v_AP_ppp)/max ) );
  distPhiDep.push_back( make_pair(-2, rho(2,0)*(v_AP_ppp*v_AP_pmp + v_AP_mpm*v_AP_mmm)/max) );
  
  return distPhiDep;
}

DecayMEPtr IFgx2ggxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_mpm = -(1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_pmp = -v_AP_mpm;
    
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1, PDT::Spin1, PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(2,2,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,2,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,2) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,2) = 0;
  (*kernelPhiDep)(2,0,0) = 0;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFgx2ggxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFgx2ggxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFgx2ggxDipoleKernel> IFgx2ggxDipoleKernel::initIFgx2ggxDipoleKernel;
// Definition of the static class description member.

void IFgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<IFgx2ggxDipoleKernel> documentation
    ("IFgx2ggxDipoleKernel");

}

#line 1 "./IFgx2qqxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFgx2qqxDipoleKernel class.
//

#include "IFgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFgx2qqxDipoleKernel::IFgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IFgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() == ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    flavour() == sk.flavour() &&
    abs(flavour()->id()) < 6 &&
    flavour()->mass() == ZERO &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IFgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IFgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());

  double rho = 1. - 4.*ratio*z*(1.-z)/sqr(1.-z+ratio);
  double x = 0.5*((1.-z+ratio)/ratio)*(1.-sqrt(rho));

  ret *= 0.5 * (!strictLargeN() ? 4./3. : 3./2.) * ( x + 2.*(1.-x)/x );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IFgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt(1./z);
  double v_AP_mpm = (1.-z)/sqrt(z);

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_pmp = -v_AP_mpm;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = sqr(v_AP_ppp) + sqr(v_AP_mpm)
    - 2.*abs(rho(0,2))*(v_AP_pmp*v_AP_ppp + v_AP_mmm*v_AP_mpm);

  distPhiDep.push_back(make_pair( 0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_mpm))/max ) );
  distPhiDep.push_back(make_pair( 2, rho(0,2)*(v_AP_pmp*v_AP_ppp + v_AP_mmm*v_AP_mpm )/max ) );
  distPhiDep.push_back(make_pair( -2, rho(2,0)*(v_AP_ppp*v_AP_pmp + v_AP_mpm*v_AP_mmm )/max ) );

  return distPhiDep;
}

DecayMEPtr IFgx2qqxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt(1./z);
  double v_AP_mpm = (1.-z)/sqrt(z);
  
  double v_AP_mmm = -v_AP_ppp;
  double v_AP_pmp = -v_AP_mpm;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half, PDT::Spin1, PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -1 or -1/2, 1=+1/2, 2 = +1
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,2,1) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,1) = 0.;
  (*kernelPhiDep)(1,2,0) = 0.;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(1,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,1) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;            
  
  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFgx2qqxDipoleKernel> IFgx2qqxDipoleKernel::initIFgx2qqxDipoleKernel;
// Definition of the static class description member.

void IFgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<IFgx2qqxDipoleKernel> documentation
    ("IFgx2qqxDipoleKernel");

}

#line 1 "./IIqx2qgxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIqx2qgxDipoleKernel class.
//

#include "IIqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IIqx2qgxDipoleKernel::IIqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IIqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IIqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IIqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool IIqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    emitter(a) == sk.emitter(b) &&
    emission(a) == sk.emission(b) &&
    a.emitterPDF() == b.emitterPDF() &&
    a.spectatorData() == b.spectatorData() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr IIqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  return ind.emitterData();
}

tcPDPtr IIqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IIqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IIqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());
  double x = z*(1.-z)/(1.-z+ratio);

  ret *= (!strictLargeN() ? 4./3. : 3./2.) * ( (1.+sqr(x))/(1.-x) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IIqx2qgxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {

  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr IIqx2qgxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt( 1./(1.-z) );
  double v_AP_ppm = -z/sqrt(1.-z);

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half, PDT::Spin1Half, PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,1,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(1,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = 0.;
  (*kernelPhiDep)(1,0,2) = 0.;
  (*kernelPhiDep)(0,1,2) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;
  
  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IIqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IIqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IIqx2qgxDipoleKernel> IIqx2qgxDipoleKernel::initIIqx2qgxDipoleKernel;
// Definition of the static class description member.

void IIqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<IIqx2qgxDipoleKernel> documentation
    ("IIqx2qgxDipoleKernel");

}

#line 1 "./IIqx2gqxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIqx2gqxDipoleKernel class.
//

#include "IIqx2gqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IIqx2gqxDipoleKernel::IIqx2gqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IIqx2gqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IIqx2gqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IIqx2gqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool IIqx2gqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    a.emitterData() == b.emitterData() &&
    emitter(a) == sk.emitter(b) &&
    a.emitterPDF() == b.emitterPDF() &&
    a.spectatorData() == b.spectatorData() &&
    a.spectatorPDF() == b.spectatorPDF();

}

  
tcPDPtr IIqx2gqxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IIqx2gqxDipoleKernel::emission(const DipoleIndex& ind) const {
  return ind.emitterData()->CC();
}

tcPDPtr IIqx2gqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IIqx2gqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());

  double x = z*(1.-z)/(1.-z+ratio);

  ret *= .5 * ( 1.-2.*x*(1.-x) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IIqx2gqxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {

  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr IIqx2gqxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppm = z;
  double v_AP_mpm = (1.-z);
  
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_pmp = -v_AP_mpm;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1, PDT::Spin1Half, PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = 0.;
  (*kernelPhiDep)(2,1,1) = 0.;
  (*kernelPhiDep)(0,0,1) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,1,1) = 0.;
  (*kernelPhiDep)(2,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IIqx2gqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IIqx2gqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IIqx2gqxDipoleKernel> IIqx2gqxDipoleKernel::initIIqx2gqxDipoleKernel;
// Definition of the static class description member.

void IIqx2gqxDipoleKernel::Init() {

  static ClassDocumentation<IIqx2gqxDipoleKernel> documentation
    ("IIqx2gqxDipoleKernel");

}

#line 1 "./IIgx2ggxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIgx2ggxDipoleKernel class.
//

#include "IIgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IIgx2ggxDipoleKernel::IIgx2ggxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IIgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IIgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IIgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool IIgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    a.emitterPDF() == b.emitterPDF() &&
    a.spectatorData() == b.spectatorData() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr IIgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IIgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IIgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IIgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());

  double x = z*(1.-z)/(1.-z+ratio);

  ret *= 3. * ( x/(1.-x) + (1.-x)/x +x*(1.-x) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IIgx2ggxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_mpm = -(1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  //double v_AP_mmp = -v_AP_ppm;
  double v_AP_pmp = -v_AP_mpm;
  
  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_mpm))
    - 2.*abs(rho(0,2))*(v_AP_ppp*v_AP_pmp + v_AP_mmm*v_AP_mpm);

  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*
                                  (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_mpm))/max ) );
  distPhiDep.push_back( make_pair(2, rho(0,2)*(v_AP_mmm*v_AP_mpm + v_AP_pmp*v_AP_ppp)/max ) );
  distPhiDep.push_back( make_pair(-2, rho(2,0)*(v_AP_ppp*v_AP_pmp + v_AP_mpm*v_AP_mmm)/max) );

  return distPhiDep;
}

DecayMEPtr IIgx2ggxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_mpm = -(1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_pmp = -v_AP_mpm;
    
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1, PDT::Spin1, PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(2,2,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,2,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,2) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,2) = 0;
  (*kernelPhiDep)(2,0,0) = 0;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IIgx2ggxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IIgx2ggxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IIgx2ggxDipoleKernel> IIgx2ggxDipoleKernel::initIIgx2ggxDipoleKernel;
// Definition of the static class description member.

void IIgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<IIgx2ggxDipoleKernel> documentation
    ("IIgx2ggxDipoleKernel");

}

#line 1 "./IIgx2qqxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIgx2qqxDipoleKernel class.
//

#include "IIgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IIgx2qqxDipoleKernel::IIgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IIgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IIgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IIgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() == ZERO &&
    ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool IIgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    flavour() == sk.flavour() &&
    abs(flavour()->id()) < 6 &&
    flavour()->mass() == ZERO &&
    a.emitterPDF() == b.emitterPDF() &&
    a.spectatorData() == b.spectatorData() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr IIgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IIgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IIgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IIgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());

  double x = z*(1.-z)/(1.-z+ratio);

  ret *= 0.5 * (!strictLargeN() ? 4./3. : 3./2.) * ( 1./x +sqr(1.-x)/x );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IIgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {
  
  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt(1./z);
  double v_AP_mpm = (1.-z)/sqrt(z);

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_pmp = -v_AP_mpm;
    
  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;  
  double max = sqr(v_AP_ppp) + sqr(v_AP_mpm)
    - 2.*abs(rho(0,2))*(v_AP_pmp*v_AP_ppp + v_AP_mmm*v_AP_mpm);

  distPhiDep.push_back(make_pair( 0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_mpm))/max ) );
  distPhiDep.push_back(make_pair( 2, rho(0,2)*(v_AP_pmp*v_AP_ppp + v_AP_mmm*v_AP_mpm )/max ) );
  distPhiDep.push_back(make_pair( -2, rho(2,0)*(v_AP_ppp*v_AP_pmp + v_AP_mpm*v_AP_mmm )/max ) );

  return distPhiDep;
}

DecayMEPtr IIgx2qqxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt(1./z);
  double v_AP_mpm = (1.-z)/sqrt(z);
  
  double v_AP_mmm = -v_AP_ppp;
  double v_AP_pmp = -v_AP_mpm;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half, PDT::Spin1, PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -1 or -1/2, 1=+1/2, 2 = +1
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,2,1) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,1) = 0.;
  (*kernelPhiDep)(1,2,0) = 0.;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(1,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,1) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;            
  
  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IIgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IIgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IIgx2qqxDipoleKernel> IIgx2qqxDipoleKernel::initIIgx2qqxDipoleKernel;
// Definition of the static class description member.

void IIgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<IIgx2qqxDipoleKernel> documentation
    ("IIgx2qqxDipoleKernel");

}

#line 1 "./FIMqx2qgxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMqx2qgxDipoleKernel class.
//

#include "FIMqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMqx2qgxDipoleKernel::FIMqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr FIMqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 7 &&
    abs(ind.emitterData()->id())==abs(flavour()->id()) &&
    ind.emitterData()->mass() != ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIMqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    abs(sk.emitter(b)->id()) < 7 &&
    sk.emitter(b)->mass() == emitter(a)->mass() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr FIMqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  assert(flavour());
  assert(abs(flavour()->id())<7 && flavour()->mass() != ZERO);
  return ind.emitterData()->id() > 0 ?
    (tcPDPtr) flavour() : (tcPDPtr) flavour()->CC();
}

tcPDPtr FIMqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIMqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

// TODO
// split.scale() should be sqrt(sbar) = sqrt( Mi2 - Q2 ) !!!
double FIMqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  // Mi=mi=mQ, m=0, Mj=mj=0
  Energy2 mQ2 = sqr(split.emitterMass());

  double z = split.lastZ();
  double x = 1. / ( 1. + 
		    ( sqr(split.lastPt()) + sqr(1.-z)*mQ2 ) /
		    ( z*(1.-z) * sqr(split.scale()) ) );

  // Simon has extra terms
  ret *= (!strictLargeN() ? 4./3. : 3./2.) *
    ( 2./(1.-z+(1.-x)) -(1.+z) - mQ2/sqr(split.scale()) * 2.*x/(1.-x) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FIMqx2qgxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {

  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr FIMqx2qgxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {

  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();
  Energy mi = dInfo.emitterMass();
  
  // Altarelli-Parisi spin-indexed kernels:
  Energy den = sqrt(sqr(mi)*sqr(1.-z) + sqr(pt));
  double v_AP_ppp = pt / den / sqrt(1.-z);
  double v_AP_ppm = - z * v_AP_ppp ;
  double v_AP_pmp = mi*(1.-z)*sqrt(1.-z) / den ;

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = v_AP_pmp;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -1 or -1/2, 1=+1/2, 2 = +1
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,1,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(1,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm;
  (*kernelPhiDep)(1,0,2) = v_AP_pmp;
  (*kernelPhiDep)(0,1,2) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIMqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIMqx2qgxDipoleKernel> FIMqx2qgxDipoleKernel::initFIMqx2qgxDipoleKernel;
// Definition of the static class description member.

void FIMqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<FIMqx2qgxDipoleKernel> documentation
    ("FIMqx2qgxDipoleKernel");

}

#line 1 "./FIMgx2qqxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMgx2qqxDipoleKernel class.
//

#include "FIMgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMgx2qqxDipoleKernel::FIMgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr FIMgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() != ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIMgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    emitter(a)->mass() == sk.emitter(b)->mass() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr FIMgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() != ZERO);
  return flavour();
}

tcPDPtr FIMgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() != ZERO);
  return flavour()->CC();
}

tcPDPtr FIMgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

// TODO
// assure split.scale() is sqrt(sbar)
double FIMgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  // mi=m=mQ, Mi=0, mj=Mj=0
  Energy2 mQ2 = sqr(split.emitterData()->mass());

  double z = split.lastZ();
  double x = 1./ ( 1. +
		   ( sqr(split.lastPt()) + mQ2 ) /
		   ( z*(1.-z) * sqr(split.scale()) ) );

  double muQ2 = x * mQ2/sqr(split.scale());

  double zm = .5 * ( 1. - sqrt( 1. - 4.*muQ2/(1.-x) ) );
  double zp = .5 * ( 1. + sqrt( 1. - 4.*muQ2/(1.-x) ) );

  ret *= .25 * (1.-2.*(zp-z)*(z-zm));

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FIMgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();
  Energy2 mi2 = sqr(dInfo.emitterData()->mass());
  
  // Altarelli-Parisi spin-indexed kernels:
  double ratio = mi2 / ( mi2 + sqr(pt) );
  double root = sqrt(1.-ratio);
  double v_AP_ppp = sqrt(ratio);
  double v_AP_ppm = z*root;
  double v_AP_pmp = -(1.-z)*root;

  //double v_AP_mmm = v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp)) + 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);
  
  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp) )/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );

  return distPhiDep;
}

DecayMEPtr FIMgx2qqxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {
  
  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();
  Energy2 mi2 = sqr(dInfo.emitterData()->mass());

  // Altarelli-Parisi spin-indexed kernels:
  double ratio = mi2 / ( mi2 + sqr(pt) );
  double root = sqrt(1.-ratio);
  double v_AP_ppp = sqrt(ratio);
  double v_AP_ppm = z*root;
  double v_AP_pmp = -(1.-z)*root;

  double v_AP_mmm = v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm;
  (*kernelPhiDep)(2,1,1) = v_AP_ppp;
  (*kernelPhiDep)(0,0,1) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,1,1) = 0.;
  (*kernelPhiDep)(2,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIMgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIMgx2qqxDipoleKernel> FIMgx2qqxDipoleKernel::initFIMgx2qqxDipoleKernel;
// Definition of the static class description member.

void FIMgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FIMgx2qqxDipoleKernel> documentation
    ("FIMgx2qqxDipoleKernel");

}

#line 1 "./IFMqx2qgxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMqx2qgxDipoleKernel class.
//

#include "IFMqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMqx2qgxDipoleKernel::IFMqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IFMqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFMqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFMqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() != ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFMqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    emitter(a) == sk.emitter(b) &&
    emission(a) == sk.emission(b) &&
    spectator(a)->mass() == sk.spectator(b)->mass() &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFMqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  return ind.emitterData();
}

tcPDPtr IFMqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFMqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFMqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  Energy pt = split.lastPt();
  double ratio = sqr(pt/split.scale());
  double muk2 = sqr(split.spectatorMass()/split.scale());

// Calculate x and u
    double rho = 1. - 4.*ratio*(1.-muk2)*z*(1.-z)/sqr(1.-z+ratio);
    double x = 0.5*((1.-z+ratio)/(ratio*(1.-muk2))) * (1. - sqrt(rho));
    double u = x*ratio / (1.-z);
  
// 19/01/2017 - SW: Removed finite term as its effect on
// the ratio of -ve/+ve kernels is small
    ret *= (!strictLargeN() ? 4./3. : 3./2.) * ( 2./(1.-x+u) - (1.+x) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IFMqx2qgxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {

  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr IFMqx2qgxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt( 1./(1.-z) );
  double v_AP_ppm = -z/sqrt(1.-z);

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half, PDT::Spin1Half, PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,1,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(1,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = 0.;
  (*kernelPhiDep)(1,0,2) = 0.;
  (*kernelPhiDep)(0,1,2) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;
  
  return kernelPhiDep;
}


// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFMqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFMqx2qgxDipoleKernel> IFMqx2qgxDipoleKernel::initIFMqx2qgxDipoleKernel;
// Definition of the static class description member.

void IFMqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<IFMqx2qgxDipoleKernel> documentation
    ("IFMqx2qgxDipoleKernel");

}

#line 1 "./IFMqx2gqxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMqx2gqxDipoleKernel class.
//

#include "IFMqx2gqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMqx2gqxDipoleKernel::IFMqx2gqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IFMqx2gqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFMqx2gqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFMqx2gqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() != ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFMqx2gqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    a.emitterData() == b.emitterData() &&
    emitter(a) == sk.emitter(b) &&
    spectator(a)->mass() == sk.spectator(b)->mass() &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFMqx2gqxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFMqx2gqxDipoleKernel::emission(const DipoleIndex& ind) const {
  return ind.emitterData()->CC();
}

tcPDPtr IFMqx2gqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFMqx2gqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  Energy pt = split.lastPt();
  double ratio = sqr(pt/split.scale());
  double muk2 = sqr(split.spectatorMass()/split.scale());
  
  // Calculate x
  double rho = 1. - 4.*ratio*(1.-muk2)*z*(1.-z)/sqr(1.-z+ratio);
  double x = 0.5*((1.-z+ratio)/(ratio*(1.-muk2))) * (1. - sqrt(rho));

  ret *= .5 * ( 1.-2.*x*(1.-x)  );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IFMqx2gqxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {

  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr IFMqx2gqxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppm = z;
  double v_AP_mpm = (1.-z);

  double v_AP_mmp = -v_AP_ppm;
  double v_AP_pmp = -v_AP_mpm;
    
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1, PDT::Spin1Half, PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = 0.;
  (*kernelPhiDep)(2,1,1) = 0.;
  (*kernelPhiDep)(0,0,1) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,1,1) = 0.;
  (*kernelPhiDep)(2,0,0) = 0.;

  return kernelPhiDep;
}


// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMqx2gqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFMqx2gqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFMqx2gqxDipoleKernel> IFMqx2gqxDipoleKernel::initIFMqx2gqxDipoleKernel;
// Definition of the static class description member.

void IFMqx2gqxDipoleKernel::Init() {

  static ClassDocumentation<IFMqx2gqxDipoleKernel> documentation
    ("IFMqx2gqxDipoleKernel");

}

#line 1 "./IFMgx2ggxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMgx2ggxDipoleKernel class.
//

#include "IFMgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMgx2ggxDipoleKernel::IFMgx2ggxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IFMgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFMgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFMgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() != ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFMgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    a.emitterPDF() == b.emitterPDF() &&
    sk.spectator(b)->mass() == spectator(a)->mass();

}


tcPDPtr IFMgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFMgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFMgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFMgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  Energy pt = split.lastPt();
  double ratio = sqr(pt/split.scale());
  double muk2 = sqr(split.spectatorMass()/split.scale());

	// Calculate x and u
    double rho = 1. - 4.*ratio*(1.-muk2)*z*(1.-z)/sqr(1.-z+ratio);
    double x = 0.5*((1.-z+ratio)/(ratio*(1.-muk2))) * (1. - sqrt(rho));
    double u = x*ratio / (1.-z);

// NOTE - The definition of muk used in the kinematics differs from that in CS

    double muk2CS = x*muk2;
    ret *= 3. * ( 1./(1.-x+u) - 1. + x*(1.-x) + (1.-x)/x - muk2CS*u/(x*(1.-u)) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IFMgx2ggxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_mpm = -(1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  //double v_AP_mmp = -v_AP_ppm;
  double v_AP_pmp = -v_AP_mpm;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_mpm)) - 2.*abs(rho(0,2))*(v_AP_ppp*v_AP_pmp + v_AP_mmm*v_AP_mpm);

  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_mpm))/max ) );
  distPhiDep.push_back( make_pair(2, rho(0,2)*(v_AP_mmm*v_AP_mpm + v_AP_pmp*v_AP_ppp)/max ) );
  distPhiDep.push_back( make_pair(-2, rho(2,0)*(v_AP_ppp*v_AP_pmp + v_AP_mpm*v_AP_mmm)/max) );
  
  return distPhiDep;
}

DecayMEPtr IFMgx2ggxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_mpm = -(1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_pmp = -v_AP_mpm;
    
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1, PDT::Spin1, PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(2,2,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,2,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,2) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,2) = 0;
  (*kernelPhiDep)(2,0,0) = 0;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMgx2ggxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFMgx2ggxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFMgx2ggxDipoleKernel> IFMgx2ggxDipoleKernel::initIFMgx2ggxDipoleKernel;
// Definition of the static class description member.

void IFMgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<IFMgx2ggxDipoleKernel> documentation
    ("IFMgx2ggxDipoleKernelv");

}

#line 1 "./IFMgx2qqxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMgx2qqxDipoleKernel class.
//

#include "IFMgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMgx2qqxDipoleKernel::IFMgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr IFMgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFMgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFMgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() != ZERO &&
    flavour()->mass() == ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFMgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    flavour() == sk.flavour() &&
    abs(flavour()->id()) < 6 &&
    flavour()->mass() == ZERO &&
    spectator(a)->mass() == sk.spectator(b)->mass() &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFMgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IFMgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IFMgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFMgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  Energy pt = split.lastPt();
  double ratio = sqr(pt/split.scale());
  double muk2 = sqr(split.spectatorMass()/split.scale());
  
// Calculate x and u
  double rho = 1. - 4.*ratio*(1.-muk2)*z*(1.-z)/sqr(1.-z+ratio);
  double x = 0.5*((1.-z+ratio)/(ratio*(1.-muk2))) * (1. - sqrt(rho));
  double u = x*ratio / (1.-z);

  // NOTE - The definition of muk used in the kinematics differs from that in CS
    double muk2CS = x*muk2;
    ret *= 0.5 * (!strictLargeN() ? 4./3. : 3./2.) *
      ( x + 2.*(1.-x)/x - 2.*muk2CS/x*u/(1.-u) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IFMgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt(1./z);
  double v_AP_mpm = (1.-z)/sqrt(z);

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_pmp = -v_AP_mpm;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = sqr(v_AP_ppp) + sqr(v_AP_mpm)
    - 2.*abs(rho(0,2))*(v_AP_pmp*v_AP_ppp + v_AP_mmm*v_AP_mpm);

  distPhiDep.push_back(make_pair( 0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_mpm))/max ) );
  distPhiDep.push_back(make_pair( 2, rho(0,2)*(v_AP_pmp*v_AP_ppp + v_AP_mmm*v_AP_mpm )/max ) );
  distPhiDep.push_back(make_pair( -2, rho(2,0)*(v_AP_ppp*v_AP_pmp + v_AP_mpm*v_AP_mmm )/max ) );

  return distPhiDep;
}

DecayMEPtr IFMgx2qqxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt(1./z);
  double v_AP_mpm = (1.-z)/sqrt(z);
  
  double v_AP_mmm = -v_AP_ppp;
  double v_AP_pmp = -v_AP_mpm;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half, PDT::Spin1, PDT::Spin1Half)));  
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -1 or -1/2, 1=+1/2, 2 = +1
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,2,1) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,1) = 0.;
  (*kernelPhiDep)(1,2,0) = 0.;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(1,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,1) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;            
  
  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFMgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFMgx2qqxDipoleKernel> IFMgx2qqxDipoleKernel::initIFMgx2qqxDipoleKernel;
// Definition of the static class description member.

void IFMgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<IFMgx2qqxDipoleKernel> documentation
    ("IFMgx2qqxDipoleKernel");

}

#line 1 "./FIMDecaygx2qqxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMDecaygx2qqxDipoleKernel class.
//

#include "FIMDecaygx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMDecaygx2qqxDipoleKernel::FIMDecaygx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IBPtr FIMDecaygx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMDecaygx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMDecaygx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.incomingDecaySpectator() && !ind.incomingDecayEmitter() &&
    ind.emitterData()->id() == ParticleID::g &&
    !(ind.spectatorData()->mass() == ZERO) &&
    // Initial state here refers to the entire event
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool FIMDecaygx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    emitter(a)->id() == sk.emitter(b)->id() &&
    abs(sk.spectator(b)->mass()) == abs(spectator(a)->mass());

}


tcPDPtr FIMDecaygx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6);
  return flavour();
}

tcPDPtr FIMDecaygx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6);
  return flavour()->CC();
}

tcPDPtr FIMDecaygx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIMDecaygx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);

  // Sudakov parameterisation variables,
  // needed to calculate y.
  double z = split.lastZ();
  Energy pt = split.lastPt();

  // Construct mass squared variables
  // Note for q->qg can use the emitterMass
  // (i.e. mass of emitter before splitting = mass of emitter after)
  Energy2 Qijk = sqr(split.scale());
  Energy2 mi2 = sqr(split.emitterData()->mass());
  Energy2 mj2 = mi2;
  Energy2 mk2 = sqr(split.recoilMass());
  Energy2 sbar = Qijk - mi2 - mj2 - mk2;

  // Calculate y
  double y = (sqr(pt) + sqr(1.-z)*mi2 + sqr(z)*mj2) / sbar / z / (1.-z);

  if( sqr(2.*mk2+sbar*(1.-y)) - 4.*mk2*Qijk < ZERO ){
    generator()->logWarning( Exception()
    << "error in FIMDecayqx2qgxDipoleKernel::evaluate -- " <<
    "mk2 " << mk2/GeV2 << "  mi2 " << mi2/GeV2 << "  y " << y << Exception::warning );
    return 0.0;
  }

  // zi, used in dipole splitting kernel
  double zi = split.lastSplittingParameters()[0];

  double vijk = sqrt( sqr(2.*mk2 + sbar*(1.-y)) - 4.*mk2*Qijk ) / sbar / (1.-y);
  double viji = sqrt( sqr(sbar*y) - 4.*sqr(mi2) ) / (sbar*y + 2.*mi2);
  
  double zip = 0.5*(1.+viji*vijk);
  double zim = 0.5*(1.-viji*vijk);

  // how to choose kappa?
  double kappa = 0.;

  ret *= 0.25 / vijk
    * ( 1. - 2.*( zi*(1.-zi) - (1.-kappa)*zip*zim - kappa*mi2/(2.*mi2 + sbar*y) ) );
  
  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FIMDecaygx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo,
                                        const RhoDMatrix& rho) const {

  // Need variables for the AP kernels
  double z = dInfo.lastSplittingParameters()[0];
  Energy pt = dInfo.lastPt();
  Energy2 mi2 = sqr(dInfo.emitterData()->mass());

  // Altarelli-Parisi spin-indexed kernels:
  double ratio = mi2 / ( mi2 + sqr(pt) );
  double root = sqrt(1.-ratio);
  double v_AP_ppp = sqrt(ratio);
  double v_AP_ppm = z*root;
  double v_AP_pmp = -(1.-z)*root;

  //  double v_AP_mmm = v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp)) + 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);
  
  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp) )/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );

  return distPhiDep;
}


DecayMEPtr FIMDecaygx2qqxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {
  
  // Need variables for the AP kernels
  double z = dInfo.lastSplittingParameters()[0];
  Energy pt = dInfo.lastPt();
  Energy2 mi2 = sqr(dInfo.emitterData()->mass());
  
  // Altarelli-Parisi spin-indexed kernels:
  double ratio = mi2 / ( mi2 + sqr(pt) );
  double root = sqrt(1.-ratio);
  double v_AP_ppp = sqrt(ratio);
  double v_AP_ppm = z*root;
  double v_AP_pmp = -(1.-z)*root;

  double v_AP_mmm = v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm;
  (*kernelPhiDep)(2,1,1) = v_AP_ppp;;
  (*kernelPhiDep)(0,0,1) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,1,1) = 0.;
  (*kernelPhiDep)(2,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMDecaygx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIMDecaygx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIMDecaygx2qqxDipoleKernel> FIMDecaygx2qqxDipoleKernel::initFIMDecaygx2qqxDipoleKernel;
// Definition of the static class description member.

void FIMDecaygx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FIMDecaygx2qqxDipoleKernel> documentation
    ("FIMDecaygx2qqxDipoleKernel");

}
#line 1 "./FIMDecayqx2qgxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMDecayqx2qgxDipoleKernel class.
//

#include "FIMDecayqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMDecayqx2qgxDipoleKernel::FIMDecayqx2qgxDipoleKernel() : DipoleSplittingKernel() {}

IBPtr FIMDecayqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMDecayqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMDecayqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.incomingDecaySpectator() && !ind.incomingDecayEmitter() &&
    abs(ind.emitterData()->id()) < 7  &&
    // This line matches to the kernel declared in a .in file for the given emitter flavour
    abs(ind.emitterData()->id()) == abs(flavour()->id()) &&
    !(ind.spectatorData()->mass() == ZERO) && 
    // Initial state here refers to the entire event
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool FIMDecayqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
						     const DipoleSplittingKernel& sk,
						     const DipoleIndex& b) const {
  
  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    abs(sk.emitter(b)->id()) < 7 &&
    abs(sk.emitter(b)->mass()) == abs(emitter(a)->mass()) &&
    abs(sk.spectator(b)->mass()) == abs(spectator(a)->mass());

}

tcPDPtr FIMDecayqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {

  assert(flavour());
  assert(abs(flavour()->id()) < 7);

  return ind.emitterData()->id() > 0 ?
    (tcPDPtr) flavour() : (tcPDPtr) flavour()->CC();
}


tcPDPtr FIMDecayqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}


tcPDPtr FIMDecayqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}


double FIMDecayqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);

  // Sudakov parameterisation variables,
  // needed to calculate y.
  double z = split.lastZ();
  Energy pt = split.lastPt();

  // Construct mass squared variables
  // Note for q->qg can use the emitterMass
  // (i.e. mass of emitter before splitting = mass of emitter after)
  Energy2 Qijk = sqr(split.scale());
  Energy2 mi2 = sqr(split.emitterMass());
  Energy2 mk2 = sqr(split.recoilMass());
  Energy2 sbar = Qijk - mi2 - mk2;

  // Note this should be the same as Qijk
  Energy2 ma2 = sqr(split.spectatorMass());

  // Calculate y
  double y = (sqr(pt) + sqr(1.-z)*mi2) / sbar / z / (1.-z);

  if( sqr(2.*mk2+sbar*(1.-y)) - 4.*mk2*Qijk < ZERO ){
    generator()->logWarning( Exception()
    << "error in FIMDecayqx2qgxDipoleKernel::evaluate -- " <<
    "mk2 " << mk2/GeV2 << "  mi2 " << mi2/GeV2 << "  y " << y << Exception::warning );
    return 0.0;
  }

  // zi, used in dipole splitting kernel
  double zi = split.lastSplittingParameters()[0];

  double vijk = sqrt( sqr(2.*mk2 + sbar*(1.-y)) - 4.*mk2*Qijk ) / sbar / (1.-y);
  double vtilde = sqrt( sqr(Qijk) + sqr(mi2) + sqr(mk2)
                        - 2.*(mi2*Qijk + mk2*Qijk + mi2*mk2) ) / sbar;

  ret *=
    (!strictLargeN() ? 4./3. : 3./2.)
    * ( ( 2.*( 2.*mi2/sbar + 2.*y + 1. ) / ((1.+y) - zi*(1.-y))
          - (vtilde/vijk) * ( (1.+zi) + 2.*mi2/y/sbar ) )
	+ y/(1. - zi*(1.-y)) * ( 2.*(2.*mi2/sbar + 2.*y + 1. ) / ((1.+y) - zi*(1.-y))
                             - (vtilde/vijk) * ( 2. + 2.*ma2/((1. - zi*(1.-y))*sbar) ) ) );
  
  return ret > 0. ? ret : 0.;
  
}

vector< pair<int, Complex> >
FIMDecayqx2qgxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {

  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.  
  return {{ {0, 1.} }};
}

DecayMEPtr FIMDecayqx2qgxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {

  // Need variables for the AP kernels
  double z = dInfo.lastSplittingParameters()[0];
  Energy pt = dInfo.lastPt();
  Energy mi = dInfo.emitterMass();

  // Altarelli-Parisi spin-indexed kernels:
  Energy den = sqrt(sqr(mi)*sqr(1.-z) + sqr(pt));
  double v_AP_ppp = pt / den / sqrt(1.-z);
  double v_AP_ppm = - z * v_AP_ppp ;
  double v_AP_pmp = mi*(1.-z)*sqrt(1.-z) / den ;

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = v_AP_pmp;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -1 or -1/2, 1=+1/2, 2 = +1
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,1,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(1,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm;
  (*kernelPhiDep)(1,0,2) = v_AP_pmp;
  (*kernelPhiDep)(0,1,2) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMDecayqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIMDecayqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIMDecayqx2qgxDipoleKernel> FIMDecayqx2qgxDipoleKernel::initFIMDecayqx2qgxDipoleKernel;
// Definition of the static class description member.

void FIMDecayqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<FIMDecayqx2qgxDipoleKernel> documentation
    ("FIMDecayqx2qgxDipoleKernel");

}

#line 1 "./FIMDecaygx2ggxDipoleKernel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMDecaygx2ggxDipoleKernel class.
//

#include "FIMDecaygx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMDecaygx2ggxDipoleKernel::FIMDecaygx2ggxDipoleKernel() 
  : DipoleSplittingKernel(){}

IBPtr FIMDecaygx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMDecaygx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMDecaygx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.incomingDecaySpectator() && !ind.incomingDecayEmitter() &&
    ind.emitterData()->id() == ParticleID::g &&
    !(ind.spectatorData()->mass() == ZERO) &&
    // Initial state here refers to the entire event
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool FIMDecaygx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    sk.emitter(b)->id() == ParticleID::g &&
    abs(sk.spectator(b)->mass()) == abs(spectator(a)->mass());

}


tcPDPtr FIMDecaygx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIMDecaygx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIMDecaygx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIMDecaygx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);

  // Sudakov parameterisation variables,
  // needed to calculate y.
  double z = split.lastZ();
  Energy pt = split.lastPt();

  // Construct mass squared variables
  // Note for q->qg can use the emitterMass
  // (i.e. mass of emitter before splitting = mass of emitter after)
  Energy2 Qijk = sqr(split.scale());
  Energy2 mk2 = sqr(split.recoilMass());
  Energy2 sbar = Qijk - mk2;

  // Note this should be the same as Qijk
  Energy2 ma2 = sqr(split.spectatorMass());


  // Calculate y
  double y = sqr(pt) / sbar / z / (1.-z);

  if( sqr(2.*mk2+sbar*(1.-y)) - 4.*mk2*Qijk < ZERO ){
    generator()->logWarning( Exception()
    << "error in FIMDecayqx2qgxDipoleKernel::evaluate -- " <<
    "mk2 " << mk2/GeV2 << "  y " << y << Exception::warning );
    return 0.0;
  }

  // zi, used in dipole splitting kernel
  double zi = split.lastSplittingParameters()[0];
  
  double vijk = sqrt( sqr(2.*mk2 + sbar*(1.-y)) - 4.*mk2*Qijk ) / sbar / (1.-y);
  double vtilde = 1.;
  double viji = 1.;

  double zip = 0.5*(1.+viji*vijk);
  double zim = 0.5*(1.-viji*vijk);

  // how to choose kappa?
  double kappa = 0.;

  double S1 = 0.5*3.*(2.*y + 1.)/((1.+y)-zi*(1.-y)) +
    (!strictLargeN() ? 4./3. : 3./2.)*
    y/(1.-zi*(1.-y)) * ( 2.*(2.*y + 1.)/((1.+y)-zi*(1.-y))
			- (vtilde/vijk)*(2. + 2.*ma2/((1.-zi*(1.-y))*sbar)) );
  double S2 = 0.5*3.*(2.*y + 1.)/((1.+y)-(1.-zi)*(1.-y)) +
    (!strictLargeN() ? 4./3. : 3./2.)*
    y/(1.-(1.-zi)*(1.-y)) * ( 2.*(2.*y + 1.)/((1.+y)-(1.-zi)*(1.-y))
			     - (vtilde/vijk)*(2. + 2.*ma2/((1.-(1.-zi)*(1.-y))*sbar)) );  
  double NS = 0.5*3.*(zi*(1.-zi)-(1.-kappa)*zip*zim - 2.)/vijk;

  if( theAsymmetryOption == 0 ){
    ret *= 2.*S1 + NS;
  }else if ( theAsymmetryOption == 1 ){
    ret *= 2.*zi*( S1 + S2 + NS );
  }else{
    ret *= S1 + S2 + NS;
  }
  
  return ret > 0. ? ret : 0.;
  
}

vector< pair<int, Complex> >
FIMDecaygx2ggxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo,
                                        const RhoDMatrix& rho) const {

  // Need variables for the AP kernels
  double z = dInfo.lastSplittingParameters()[0];

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_pmp = (1.-z)*sqrt( (1.-z)/z );

  //double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp)) - 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);

  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp))/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );
  
  return distPhiDep;
}


DecayMEPtr FIMDecaygx2ggxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  // Need variables for the AP kernels
  double z = dInfo.lastSplittingParameters()[0];  
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_pmp = (1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(2,2,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,2,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,2) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,2) = 0;
  (*kernelPhiDep)(2,0,0) = 0;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMDecaygx2ggxDipoleKernel::persistentOutput(PersistentOStream & os) const {
  os<<theAsymmetryOption;
}

void FIMDecaygx2ggxDipoleKernel::persistentInput(PersistentIStream & is, int) {
  is>>theAsymmetryOption;
}

ClassDescription<FIMDecaygx2ggxDipoleKernel> FIMDecaygx2ggxDipoleKernel::initFIMDecaygx2ggxDipoleKernel;
// Definition of the static class description member.

void FIMDecaygx2ggxDipoleKernel::Init() {

  static ClassDocumentation<FIMDecaygx2ggxDipoleKernel> documentation
    ("FIMDecaygx2ggxDipoleKernel");

  static Parameter<FIMDecaygx2ggxDipoleKernel,int> interfacetheAsymmetryOption
    ("AsymmetryOption",
     "The asymmetry option for final state gluon spliitings.",
     &FIMDecaygx2ggxDipoleKernel::theAsymmetryOption, 0, 0, 0,
     false, false, Interface::lowerlim);

}
