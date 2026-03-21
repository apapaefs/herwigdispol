#line 1 "./UEDBase.cc"
// -*- C++ -*-
//
// UEDBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDBase class.
//

#include "UEDBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Repository/Repository.h" 
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

UEDBase::UEDBase() : theRadCorr(true), theInvRadius(500.*GeV), 
		     theLambdaR(20.), theMbarH(), theSinThetaOne(0.),
		     theVeV(246.*GeV), includeSMMass_(true), fixedCouplings_(false), includeGaugeMixing_(true)
{}

void UEDBase::doinit() {
  readDecays(false);
  BSMModel::doinit();
  //level-1 masses and mixing angle
  calculateKKMasses(1);
  writeSpectrum();
  //add the level-1 vertices.
  addVertex(theF1F1Z0Vertex);
  addVertex(theF1F1G0Vertex);
  addVertex(theF1F0G1Vertex);
  addVertex(theG1G1G0Vertex);
  addVertex(theG0G0G1G1Vertex);
  addVertex(theF1F1P0Vertex);
  addVertex(theF1F1W0Vertex);
  addVertex(theF1F0W1Vertex);
  addVertex(theF1F0H1Vertex);
  addVertex(theP0H1H1Vertex);
  addVertex(theZ0H1H1Vertex);
  addVertex(theW0A1H1Vertex);
  addVertex(theZ0A1h1Vertex);
  addVertex(theW0W1W1Vertex);
  readDecays(true);
  if(decayFile()=="") return;
  decayRead();
}

void UEDBase::persistentOutput(PersistentOStream & os) const {
  os << theRadCorr << ounit(theInvRadius, GeV) << theLambdaR 
     << theF1F1Z0Vertex << theF1F1G0Vertex << theF1F0G1Vertex
     << theG1G1G0Vertex << theG0G0G1G1Vertex << theF1F1P0Vertex
     << theF1F1W0Vertex << theF1F0W1Vertex << theF1F0H1Vertex 
     << theP0H1H1Vertex << theZ0H1H1Vertex << theW0A1H1Vertex 
     << theZ0A1h1Vertex << theW0W1W1Vertex << ounit(theVeV,GeV) 
     << ounit(theMbarH, GeV) << theSinThetaOne << includeSMMass_
     << fixedCouplings_ << includeGaugeMixing_;
}

void UEDBase::persistentInput(PersistentIStream & is, int) {
  is >> theRadCorr >> iunit(theInvRadius, GeV) >> theLambdaR
     >> theF1F1Z0Vertex >> theF1F1G0Vertex >> theF1F0G1Vertex
     >> theG1G1G0Vertex >> theG0G0G1G1Vertex >> theF1F1P0Vertex
     >> theF1F1W0Vertex >> theF1F0W1Vertex >> theF1F0H1Vertex 
     >> theP0H1H1Vertex >> theZ0H1H1Vertex >> theW0A1H1Vertex 
     >> theZ0A1h1Vertex >> theW0W1W1Vertex >> iunit(theVeV,GeV) 
     >> iunit(theMbarH, GeV) >> theSinThetaOne >> includeSMMass_
     >> fixedCouplings_ >> includeGaugeMixing_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDBase,BSMModel>
describeHerwigUEDBase("Herwig::UEDBase", "HwUED.so");

void UEDBase::Init() {

  static ClassDocumentation<UEDBase> documentation
    ("This class implements/stores the necessary information for the simulation"
     " of a Universal Extra Dimensions model.",
     "Universal extra dimensions model based on \\cite{Cheng:2002iz,Appelquist:2000nn}.",
     "%\\cite{Cheng:2002iz}\n"
     "\\bibitem{Cheng:2002iz}\n"
     "  H.~C.~Cheng, K.~T.~Matchev and M.~Schmaltz,\n"
     "  ``Radiative corrections to Kaluza-Klein masses,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 66}, 036005 (2002)\n"
     "  [arXiv:hep-ph/0204342].\n"
     "  %%CITATION = PHRVA,D66,036005;%%\n"
     "%\\cite{Appelquist:2000nn}\n"
     "\\bibitem{Appelquist:2000nn}\n"
     "  T.~Appelquist, H.~C.~Cheng and B.~A.~Dobrescu,\n"
     "  ``Bounds on universal extra dimensions,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 64}, 035002 (2001)\n"
     "  [arXiv:hep-ph/0012100].\n"
     "  %%CITATION = PHRVA,D64,035002;%%\n"
   );

  static Switch<UEDBase,bool> interfaceRadiativeCorrections
    ("RadiativeCorrections",
     "Calculate the radiative corrections to the masses",
     &UEDBase::theRadCorr, true, false, false);
  static SwitchOption interfaceRadiativeCorrectionsYes
    (interfaceRadiativeCorrections,
     "Yes",
     "Calculate the radiative corrections to the masses",
     true);
  static SwitchOption interfaceRadiativeCorrectionsNo
    (interfaceRadiativeCorrections,
     "No",
     "Leave the masses of the KK particles as n/R",
     false);

  static Parameter<UEDBase,Energy> interfaceInverseRadius
    ("InverseRadius",
     "The inverse radius of the compactified dimension ",
     &UEDBase::theInvRadius, GeV, 500.*GeV, ZERO, ZERO,
     true, false, Interface::nolimits);

  static Parameter<UEDBase,double> interfaceLambdaR
    ("LambdaR",
     "The product of the cut-off scale  and the radius of compactification",
     &UEDBase::theLambdaR, 20.0, 0.0, 0,
     false, false, Interface::lowerlim);

    static Parameter<UEDBase,Energy> interfaceBoundaryMass
    ("HiggsBoundaryMass",
     "The boundary mass for the Higgs",
     &UEDBase::theMbarH, GeV, ZERO, ZERO, ZERO,
     false, false, Interface::lowerlim);

  static Parameter<UEDBase,Energy> interfaceVeV
    ("HiggsVEV",
     "The vacuum expectation value of the Higgs field",
     &UEDBase::theVeV, GeV, 246.*GeV, ZERO, ZERO,
     true, false, Interface::nolimits);
    
  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F1Z
    ("Vertex/F1F1Z",
     "The F1F1Z UED Vertex",
     &UEDBase::theF1F1Z0Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F1G0
    ("Vertex/F1F1G0",
     "The F1F1G UED Vertex",
     &UEDBase::theF1F1G0Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F0G1
    ("Vertex/F1F0G1",
     "The F1F0G0 UED Vertex",
     &UEDBase::theF1F0G1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVVVVertex> interfaceG1G1G0
    ("Vertex/G1G1G0",
     "The G1G1G0 UED Vertex",
     &UEDBase::theG1G1G0Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVVVVVertex> interfaceG0G0G1G1
    ("Vertex/G0G0G1G1",
     "The G0G0G1G1 UED Vertex",
     &UEDBase::theG0G0G1G1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F1P
    ("Vertex/F1F1P",
     "The F1F1P UED Vertex",
     &UEDBase::theF1F1P0Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F1W
    ("Vertex/F1F1W",
     "The F1F1W UED Vertex",
     &UEDBase::theF1F1W0Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F0W1
    ("Vertex/F1F0W1",
     "The F1F0W1 UED Vertex",
     &UEDBase::theF1F0W1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFSVertex> interfaceF1F0H1
    ("Vertex/F1F0H1",
     "The F1F0H1 UED Vertex",
     &UEDBase::theF1F0H1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVSSVertex> interfaceP0H1H1
    ("Vertex/P0H1H1",
     "The P0H1H1 UED Vertex",
     &UEDBase::theP0H1H1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVSSVertex> interfaceZ0H1H1
    ("Vertex/Z0H1H1",
     "The Z0H1H1 UED Vertex",
     &UEDBase::theZ0H1H1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVSSVertex> interfaceW0A1H1
    ("Vertex/W0A1H1",
     "The W0A1H1 UED Vertex",
     &UEDBase::theW0A1H1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVSSVertex> interfaceZ0A1h1
    ("Vertex/Z0A1h1",
     "The W0A1H1 UED Vertex",
     &UEDBase::theZ0A1h1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVVVVertex> interfaceW0W1W1
    ("Vertex/W0W1W1",
     "The W0W1W1 UED Vertex",
     &UEDBase::theW0W1W1Vertex, false, false, true, false, false);

  static Switch<UEDBase,bool> interfaceIncludeSMMass
    ("IncludeSMMass",
     "Whether or not to include the SM mass in the calculation of the masses of the KK states.",
     &UEDBase::includeSMMass_, true, false, false);
  static SwitchOption interfaceIncludeSMMassYes
    (interfaceIncludeSMMass,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncludeSMMassNo
    (interfaceIncludeSMMass,
     "No",
     "Don't include them",
     false);

  static Switch<UEDBase,bool> interfaceFixedCouplings
    ("FixedCouplings",
     "Use fixed or running couplings to calculate the masses.",
     &UEDBase::fixedCouplings_, false, false, false);
  static SwitchOption interfaceFixedCouplingsYes
    (interfaceFixedCouplings,
     "Yes",
     "Use fixed couplings",
     true);
  static SwitchOption interfaceFixedCouplingsNo
    (interfaceFixedCouplings,
     "No",
     "Use running couplings",
     false);

  static Switch<UEDBase,bool> interfaceIncludeGaugeMixing
    ("IncludeGaugeMixing",
     "Whether or not to include mixing between the KK photon"
     " and Z in the vertices, always included in the mass",
     &UEDBase::includeGaugeMixing_, true, false, false);
  static SwitchOption interfaceIncludeGaugeMixingYes
    (interfaceIncludeGaugeMixing,
     "Yes",
     "Include the mixing",
     true);
  static SwitchOption interfaceIncludeGaugeMixingNo
    (interfaceIncludeGaugeMixing,
     "No",
     "Don't include the mixing",
     false);

}

void UEDBase::calculateKKMasses(const unsigned int n) {
  useMe();
  if(n == 0)
    throw InitException() << "UEDBase::resetKKMasses - "
			  << "Trying to reset masses with KK number == 0!"
			  << Exception::warning;
  if(theRadCorr) {
    fermionMasses(n);
    bosonMasses(n);
  }
  else {
    cerr << 
      "Warning: Radiative corrections to particle masses have been "
      "turned off.\n  The masses will be set to (n/R + m_sm)^1/2 and "
      "the spectrum will be\n  highly degenerate so that no decays "
      "will occur.\n  This is only meant to be used for debugging "
      "purposes.\n";
    //set masses to tree level for each kk mode
    long level1 = 5000000 + n*100000;
    long level2 = 6000000 + n*100000;
    Energy2 ndmass2 = sqr(n*theInvRadius);
    for ( int i = 1; i < 38; ++i ) {
      if(i == 7 || i == 17) i += 4;
      if(i == 26) i += 10;
      Energy kkmass = sqrt( ndmass2 + sqr(getParticleData(i)->mass()) );
      resetMass(level1 + i, kkmass);
      if( i < 7 || i == 11 || i == 13 || i == 15 )
	resetMass(level2 + i, kkmass);
    }
  }
}

void UEDBase::bosonMasses(const unsigned int n) {
  using Constants::zeta3;
  using Constants::pi;
  // Common constants
  const Energy2 invRad2 = theInvRadius*theInvRadius;
  const double g_em2 = fixedCouplings_ ? 
    4.*Constants::pi*alphaEMMZ() : 4.*pi*alphaEM(invRad2);
  const double g_s2 = fixedCouplings_ ? 
    4.*Constants::pi*alphaS()    : 4.*pi*alphaS(invRad2);
  const double g_W2 = g_em2/sin2ThetaW();
  const double g_P2 = g_em2/(1-sin2ThetaW());
  //Should probably use a function  to calculate zeta.
  const Energy2 nmass2 = sqr(n*theInvRadius); 
  const double pi2 = sqr(pi);
  const double norm = 1./16./pi2;
  const double nnlogLR = n*n*log(theLambdaR);
  long level = 5000000 + n*100000;
  //gluon
  Energy2 deltaGB = g_s2*invRad2*norm*(23.*nnlogLR - 3.*zeta3/2./pi2 );
  resetMass(level + 21, sqrt(nmass2 + deltaGB));

  //W+/-
  Energy2 deltaGW = g_W2*invRad2*norm*( 15.*nnlogLR - 5.*zeta3/2./pi2 );

  Energy2 mw2 = sqr(getParticleData(24)->mass());
  resetMass(level + 24, sqrt(mw2 + nmass2 + deltaGW));

  //Z and gamma are a mixture of Bn and W3n
  deltaGB = -g_P2*invRad2*norm*( 39.*zeta3/2./pi2 + nnlogLR/3. );
  Energy2 mz2 = sqr(getParticleData(23)->mass());
  Energy2 fp = 0.5*(mz2 + deltaGB + deltaGW + 2.*nmass2);
  Energy2 sp = 0.5*sqrt( sqr(deltaGB - deltaGW - 2.*mw2 + mz2)
			 - 4.*mw2*(mw2 - mz2) );
  resetMass(level + 22, sqrt(fp - sp));
  resetMass(level + 23, sqrt(fp + sp));
  //mixing angle will now depend on both Z* and gamma* mass
  //Derived expression:
  // 
  // cos^2_theta_N = ( (n/R)^2 + delta_GW + mw^2 - m_gam*^2)/(m_z*^2 - m_gam*^2)
  //
  if(includeGaugeMixing_) {
    double cn2 = (nmass2 + deltaGW + mw2 - fp + sp)/2./sp;
    double sn = sqrt(1. - cn2);
    theMixingAngles.insert(make_pair(n, sn));
    if( n == 1 ) theSinThetaOne = sn;
  }
  else {
    theMixingAngles.insert(make_pair(n,0.));
    if( n == 1 ) theSinThetaOne = 0.;
  }
  //scalars
  Energy2 mh2 = sqr(getParticleData(25)->mass());
  double lambda_H = mh2/sqr(theVeV);
  deltaGB = nnlogLR*norm*invRad2*(3.*g_W2 + (3.*g_P2/2.) - 2.*lambda_H) 
    + sqr(theMbarH);
  //H0
  Energy2 new_m2 = nmass2 + deltaGB;
  resetMass(level + 25, sqrt( mh2 + new_m2 ));
  //A0
  resetMass(level + 36, sqrt( mz2 + new_m2 ));
  //H+
  resetMass(level + 37, sqrt( mw2 + new_m2 ));
}

void UEDBase::fermionMasses(const unsigned int n) {
  using Constants::pi;
  const Energy2 invRad2 = theInvRadius*theInvRadius;
  const double g_em2 = fixedCouplings_ ? 
    4.*pi*alphaEMMZ() : 4.*pi*alphaEM(invRad2);
  const double g_s2 = fixedCouplings_ ? 
    4.*pi*alphaS() : 4.*pi*alphaS(invRad2); 
  const double g_W2 = g_em2/sin2ThetaW();
  const double g_P2 = g_em2/(1-sin2ThetaW());
  const Energy nmass = n*theInvRadius;
  const Energy norm = 
    nmass*log(theLambdaR)/16./pi/pi;
  const Energy topMass = getParticleData(6)->mass();
  const double ht = sqrt(2)*topMass/theVeV;
  //doublets
  Energy deltaL = norm*(6.*g_s2 + (27.*g_W2/8.) + (g_P2/8.));
  Energy deltaQ = deltaL;
  Energy2 shift = sqr(nmass + deltaL);
  long level = 5000000 + n*100000;
  for(long i = 1; i < 17; ++i) {
    if(i == 5)  {
      i += 6;
      deltaL = norm*( (27.*g_W2/8.) + (9.*g_P2/8.) );
      shift = sqr(nmass + deltaL); 
    }
    Energy2 new_m2 = includeSMMass_ ? sqr(getParticleData(i)->mass()) + shift : shift;
    resetMass(level + i, sqrt(new_m2));
  }
  //singlet shifts
  const Energy  deltaU = norm*(6.*g_s2 + 2.*g_P2);
  const Energy  deltaD = norm*(6.*g_s2 + 0.5*g_P2);
  const Energy2 shiftU = sqr(nmass + deltaU);
  const Energy2 shiftD = sqr(nmass + deltaD);
  
  //Top quarks seperately as they have different corrections
  const Energy2 mt2 = sqr(topMass);
  const Energy delta_Q3 = -3.*ht*ht*norm/2.;
  const Energy deltaTD = deltaQ + delta_Q3;
  const Energy deltaTS = deltaU + 2.*delta_Q3;
  Energy second_term = 0.5*sqrt( sqr(2.*nmass + deltaTS + deltaTD) + 4.*mt2 );
  //doublet
  resetMass(level           + 6, abs(0.5*(deltaTD - deltaTS) - second_term) );
  //singlet
  resetMass(level + 1000000 + 6, 0.5*(deltaTD - deltaTS) + second_term);
  //Bottom quarks
  const Energy2 mb2 = sqr(getParticleData(5)->mass());
  const Energy deltaBS = deltaD;
  second_term = 0.5*sqrt( sqr(2.*nmass + deltaBS + deltaTD) + 4.*mb2 );
  //doublet
  resetMass(level + 1000000 + 5, abs(0.5*(deltaTD - deltaBS) - second_term) );
  //singlet
  resetMass(level           + 5, 0.5*(deltaTD - deltaBS) + second_term);
  // others
  //lepton 
  Energy delta = 9.*norm*g_P2/2.;
  shift = sqr(nmass + delta);

  level += 1000000;
  for(long i = 1; i < 17; ) {
    if(i == 5) i += 6;
    Energy2 smMass2(sqr(getParticleData(i)->mass()));
    if(i < 6) {
      Energy2 new_m2 = includeSMMass_ ? smMass2 : ZERO;
      if( i % 2 == 0) new_m2 = shiftU;
      else            new_m2 = shiftD;
      resetMass(level + i, sqrt(new_m2));
      ++i;
    }
    else {
      if(includeSMMass_)
	resetMass(level + i, sqrt(smMass2 + shift));
      else
	resetMass(level + i, sqrt(shift));
      i += 2;
    }
  }
}

void UEDBase::resetMass(long id, Energy mass) {
  theMasses.push_back(make_pair(id, mass));
  StandardModel::resetMass(id,mass);
}

void UEDBase::writeSpectrum() {
  sort(theMasses.begin(), theMasses.end(), lowerMass);
  ostream & ofs = CurrentGenerator::current().misc();
  ofs << "# MUED Model Particle Spectrum\n"
      << "# R^-1: " << theInvRadius/GeV << " GeV\n"
      << "# Lambda * R: " << theLambdaR << "\n"
      << "# Higgs Mass: " << getParticleData(25)->mass()/GeV << " GeV\n";
  ofs << "#\n# ID\t\t\tMass(GeV)\n";
  while (!theMasses.empty()) {
    IDMassPair tmp = theMasses.back();
    tcPDPtr data = getParticleData(tmp.first);
    ofs << tmp.first << "\t\t\t" << tmp.second/GeV << "\t\t" << (data? data->PDGName() : "") 
	<< endl;
    theMasses.pop_back();
  }
  ofs << "#\n";
}

double UEDBase::sinThetaN(const unsigned int n) const {
  WAMap::const_iterator pos = theMixingAngles.find(n);
  if(pos != theMixingAngles.end())
    return pos->second;
  else {
    throw Exception() << "UEDBase::sinThetaN() - A mixing angle has "
		      << "been requested for a level that does not "
		      << "exist. Check that the radiative corrections "
		      << "for the " << n << "th level have been "
		      << "calculated." << Exception::warning;
    return 0.0;
  }
}
#line 1 "./UEDF1F1Z0Vertex.cc"
// -*- C++ -*-
//
// UEDF1F1Z0Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F1Z0Vertex class.
//

#include "UEDF1F1Z0Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F1Z0Vertex::UEDF1F1Z0Vertex() : theSin2ThW(0.0), theCosThW(0.0), theRadius(),
				     theID1Last(0), theID2Last(0) ,
				     theq2Last(ZERO), theCoupLast(0.), 
				     theLeftLast(0.), theRightLast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

void UEDF1F1Z0Vertex::doinit() {
  long boson = 23;
  //QQ, uu, dd
  for(long i = 5100001; i < 6100007; ++i) {
    if(i == 5100007) i += 999994;
    addToList(-i, i, boson);
  }
  //top/bottom quark l/r mixing
  addToList(-5100006, 6100006, boson); 
  addToList(-6100006, 5100006, boson); 
  addToList(-5100005, 6100005, boson); 
  addToList(-6100005, 5100005, boson); 
  //leptons
  for(long i = 5100011; i < 5100017; ++i) {
    addToList(-i, i, boson);
  }
  for(long i = 6100011; i < 6100017; i +=2) {
    addToList(-i, i, boson);
  }
  FFVVertex::doinit();
  UEDBasePtr model = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "UEDF1F1Z0Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  
  theSin2ThW = sin2ThetaW();
  theCosThW = sqrt(1. - theSin2ThW); 
  theRadius = model->compactRadius();
}

void UEDF1F1Z0Vertex::persistentOutput(PersistentOStream & os) const {
  os << theSin2ThW << theCosThW << ounit(theRadius,1/GeV);
}

void UEDF1F1Z0Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theSin2ThW >> theCosThW >> iunit(theRadius,1/GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDF1F1Z0Vertex,FFVVertex>
describeHerwigUEDF1F1Z0Vertex("Herwig::UEDF1F1Z0Vertex", "HwUED.so");

void UEDF1F1Z0Vertex::Init() {

  static ClassDocumentation<UEDF1F1Z0Vertex> documentation
    ("This is the implementation of the level-1 fermion pair Z_0 boson "
     "coupling.");

}

void UEDF1F1Z0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  if( part3->id() != 23 ) {
    norm(0.0);
    left(0.0);
    right(0.0);  
    throw HelicityLogicalError()
      << "UEDF1F1Z0Vertex::setCoupling - The vector boson in this vertex "
      << "is not a Z^0 boson. ID: " << part3->id() << "\n"
      << Exception::warning;
    return;
  }
  long ianti(abs(part1->id())), iferm(abs(part2->id()));
  bool ferma = (iferm >= 5100001 && iferm <= 5100006) ||
    (iferm >= 6100001 && iferm <= 6100006) || 
    (iferm >= 5100011 && iferm <= 5100016) ||
    (iferm >= 6100011 && iferm <= 6100016); 
  bool fermb = (ianti >= 5100001 && ianti <= 5100006) ||
    (ianti >= 6100001 && ianti <= 6100006) || 
    (ianti >= 5100011 && ianti <= 5100016) ||
    (ianti >= 6100011 && ianti <= 6100016);
  if( ferma && fermb  ) {
    if(q2 != theq2Last || theCoupLast == 0. ) {
	theq2Last = q2;
	theCoupLast = 0.5*weakCoupling(q2)/theCosThW;
    }
    if( ianti != theID1Last || iferm != theID2Last) {
      theID1Last = ianti;
      theID2Last = iferm;
      int stateA = ianti/1000000;
      int stateB = iferm/1000000;
      long smID = (stateA == 6) ? ianti - 6100000 : ianti - 5100000;
      // L/R mixing
      double alpha = atan(getParticleData(smID)->mass()*theRadius)/2.;
      double sin2al = sqr(sin(alpha));
      double cos2al = 1. - sin2al;
      
      if(stateA == 5 && stateB == 5) {
	if(smID >= 11 && smID <= 16)
	  theLeftLast = -cos2al + 2.*theSin2ThW;
	else if(smID <= 6 && smID % 2 == 0)
	  theLeftLast = cos2al - 4.*theSin2ThW/3.;
	else
	  theLeftLast = -cos2al + 2.*theSin2ThW/3.;

	theRightLast = theLeftLast;
      }
      else if(stateA == 6 && stateB == 6) {
	if(smID >= 11 && smID <= 16)
	  theLeftLast = -sin2al + 2.*theSin2ThW;
	else if(smID <=6 && smID % 2 == 0)
	  theLeftLast = sin2al - 4.*theSin2ThW/3.;
	else
	  theLeftLast = -sin2al + 2.*theSin2ThW/3.;

	theRightLast = theLeftLast;
      }
      else {
	theLeftLast = sqrt(sin2al*cos2al);
	if(smID % 2 == 0) theLeftLast *= -1.;
	theRightLast = -theLeftLast;
      }
    }
    norm(theCoupLast);
    left(theLeftLast);
    right(theRightLast);
  }
  else {
    throw HelicityLogicalError() << "UEDF1F1Z0Vertex::setCoupling - "
				 << "There is an unknown particle(s) in the "
				 << "UED F^(1) F^(1) Z^(0) vertex. ID: " 
				 << ianti << " " << iferm 
				 << Exception::warning;      
    norm(0.0);
    left(0.0);
    right(0.0);  
  }
}
#line 1 "./UEDF1F1G0Vertex.cc"
// -*- C++ -*-
//
// UEDF1F1G0Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F1G0Vertex class.
//

#include "UEDF1F1G0Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F1G0Vertex::UEDF1F1G0Vertex() 
  : theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TFUND);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<UEDF1F1G0Vertex,FFVVertex>
describeHerwigUEDF1F1G0Vertex("Herwig::UEDF1F1G0Vertex", "HwUED.so");

void UEDF1F1G0Vertex::Init() {

  static ClassDocumentation<UEDF1F1G0Vertex> documentation
    ("This class implements the F^1 F^1 G^0 vertex.");

}

void UEDF1F1G0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long iferm;
  if(part1->id() == ParticleID::g)
    iferm = abs(part2->id());
  else if(part2->id() == ParticleID::g)
    iferm = abs(part1->id());
  else if(part3->id() == ParticleID::g)
    iferm = abs(part1->id());
  else
    throw HelicityLogicalError() << "UEDF1F1G0Vertex::setCoupling - "
				 << "There is no gluon in this vertex!"
				 << Exception::warning;
  if((iferm >= 5100001 && iferm <= 5100006) ||
     (iferm >= 6100001 && iferm <= 6100006)) {
    if(q2 != theq2Last || theCoupLast ==0. ) {
      theCoupLast = -strongCoupling(q2);
      theq2Last=q2;
    }
    norm(theCoupLast);
    left(1.);
    right(1.);
  }
  else
    throw HelicityLogicalError() << "UEDF1F1G0Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex! "
				 << iferm
				 << Exception::warning;
}
void UEDF1F1G0Vertex::doinit() {
  long boson = 21;
  //QQ
  for(long i = 5100001; i < 6100007; ++i) {
    if(i == 5100007) i += 999994;
    addToList(-i, i, boson);
  }
  FFVVertex::doinit();
}
#line 1 "./UEDF1F0G1Vertex.cc"
// -*- C++ -*-
//
// UEDF1F0G1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0G1Vertex class.
//

#include "UEDF1F0G1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F0G1Vertex::UEDF1F0G1Vertex() : theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TFUND);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<UEDF1F0G1Vertex,FFVVertex>
describeHerwigUEDF1F0G1Vertex("Herwig::UEDF1F0G1Vertex", "HwUED.so");

void UEDF1F0G1Vertex::Init() {

  static ClassDocumentation<UEDF1F0G1Vertex> documentation
    ("This class implements the F^1-F^0-G^1 vertex.");

}

void UEDF1F0G1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long ifermN;
  if(part1->id() == 5100021) {
    if(abs(part2->id()) > 5100000)
      ifermN = part2->id();
    else
      ifermN = part3->id();
  }
  else if(part2->id() == 5100021) {
    if(abs(part1->id()) > 5100000)
      ifermN = part1->id();
    else
      ifermN = part3->id();
  }
  else if(part3->id() == 5100021) {
    if(abs(part2->id()) > 5100000)
      ifermN = part2->id();
    else
      ifermN = part1->id();
  }
  else
    throw HelicityLogicalError() << "UEDF1F0G1Vertex::setCoupling - "
				 << "There is no KK gluon in this vertex!"
				 << Exception::warning;
  if((abs(ifermN) >= 5100001 && abs(ifermN) <= 5100006) ||
     (abs(ifermN) >= 6100001 && abs(ifermN) <= 6100006)) {
    if(q2 != theq2Last || theCoupLast == 0.) {
      theq2Last = q2;
      theCoupLast = -strongCoupling(q2);
    }
    norm(theCoupLast);
    int state = abs(ifermN)/1000000;
    if(state == 5) {
      left(1.);
      right(0.);
    }
    else {
      left(0.);
      right(1.);
    }
  }
  else
    throw HelicityLogicalError() << "UEDF1F0G1Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex! "
				 << ifermN
				 << Exception::warning;
}

void UEDF1F0G1Vertex::doinit() {
  long boson = 5100021;
  //QQ
  for(long i = 1; i < 7; ++i) {
    addToList(-i, i + 5100000, boson);
    addToList(-(i + 5100000), i, boson);

    addToList(-i, i + 6100000, boson);
    addToList(-(i + 6100000), i, boson);
  }
  FFVVertex::doinit();
}
#line 1 "./UEDG1G1G0Vertex.cc"
// -*- C++ -*-
//
// UEDG1G1G0Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDG1G1G0Vertex class.
//

#include "UEDG1G1G0Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDG1G1G0Vertex::UEDG1G1G0Vertex() 
  : theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3F);
}

void UEDG1G1G0Vertex::doinit() {
  long kkg1 = 5100021;
  addToList(kkg1, kkg1, 21);
  VVVVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<UEDG1G1G0Vertex,Helicity::VVVVertex>
describeUEDG1G1G0Vertex("Herwig::UEDG1G1G0Vertex", "HwUED.so");

void UEDG1G1G0Vertex::Init() {

  static ClassDocumentation<UEDG1G1G0Vertex> documentation
    ("The UEDG1G1G0Vertex class implements the coupling of the "
     "gluon to two KK excitations of the gluon in the UED model.");

}

void UEDG1G1G0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2, 
				  tcPDPtr part3) {
  long id1(part1->id()), id2(part2->id()), id3(part3->id());
  if( (id1 == ParticleID::g && id2 == 5100021 && id3 == 5100021) ||
      (id2 == ParticleID::g && id1 == 5100021 && id3 == 5100021) ||
      (id3 == ParticleID::g && id1 == 5100021 && id2 == 5100021) ) {
    if(q2 != theq2Last || theCoupLast == 0.) {
      theq2Last = q2;
      theCoupLast = strongCoupling(q2);
    }
    norm(theCoupLast);
  }
  else throw HelicityLogicalError() 
    << "UEDG1G1G0Vertex::setCoupling - "
    << "There is an unknown particle in this vertex "
    << id1 << " " << id2 << " " << id3 << Exception::runerror;
}
#line 1 "./UEDG0G0G1G1Vertex.cc"
// -*- C++ -*-
//
// UEDG0G0G1G1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDG0G0G1G1Vertex class.
//

#include "UEDG0G0G1G1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "UEDBase.h"

using namespace Herwig;

UEDG0G0G1G1Vertex::UEDG0G0G1G1Vertex() : 
  theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(2);
  orderInGem(0);
  colourStructure(ColourStructure::SU3FF);
}

void UEDG0G0G1G1Vertex::doinit() {
  long kk1g = 5100021, smgl = 21;
  addToList(smgl, smgl, kk1g, kk1g);
  VVVVVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<UEDG0G0G1G1Vertex,Helicity::VVVVVertex>
describeUEDG0G0G1G1Vertex("Herwig::UEDG0G0G1G1Vertex", "HwUED.so");

void UEDG0G0G1G1Vertex::Init() {

  static ClassDocumentation<UEDG0G0G1G1Vertex> documentation
    ("This class implements the coupling of a pair of SM gluons to"
     "a pair of UED level-1 KK gluons.");

}

void UEDG0G0G1G1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2, 
				    tcPDPtr part3, tcPDPtr part4) {
  int ismg(0), ikkg(0);
  const array<tcPDPtr, 4> particles{{ part1, part2, part3, part4 }};
  for(auto p : particles) {
    if(p->id() == ParticleID::g) ++ismg;
    if(p->id() == 5100021) ++ikkg;
  }
  assert(ismg == 2 && ikkg == 2);
  if(q2 != theq2Last || theCoupLast == 0. ) {
    theq2Last = q2;
    theCoupLast = sqr(strongCoupling(q2));
  }
  norm(theCoupLast);
  setType(1);
  setOrder(0,1,2,3);
}
#line 1 "./UEDF1F1P0Vertex.cc"
// -*- C++ -*-
//
// UEDF1F1P0Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F1P0Vertex class.
//

#include "UEDF1F1P0Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F1P0Vertex::UEDF1F1P0Vertex() : coupLast_(0.0), q2Last_(ZERO),
				     fermLast_(0), LRLast_(0.0), 
				     charges_(3) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

void UEDF1F1P0Vertex::persistentOutput(PersistentOStream & os) const {
  os << charges_;
}

void UEDF1F1P0Vertex::persistentInput(PersistentIStream & is, int) {
  is >> charges_;
}

void UEDF1F1P0Vertex::doinit() {
  long photon = 22;
  //quarks
  for(long i = 1; i < 7; ++i) {
    //left
    addToList(-5100000 - i, 5100000 + i, photon);
    //right
    addToList(-6100000 - i, 6100000 + i, photon);
  }
  //leptons
  for(long i = 11; i < 17; i += 2) {
    //left
    addToList(-5100000 - i, 5100000 + i, photon);
    //right
    addToList(-6100000 - i, 6100000 + i, photon);
  }
  FFVVertex::doinit();
  tUEDBasePtr UEDBase = 
    dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!UEDBase)
    throw InitException() << "UEDF1F1P0Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  charges_[0] = UEDBase->ee();
  charges_[1] = UEDBase->ed();
  charges_[2] = UEDBase->eu();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDF1F1P0Vertex,FFVVertex>
describeHerwigUEDF1F1P0Vertex("Herwig::UEDF1F1P0Vertex", "HwUED.so");

void UEDF1F1P0Vertex::Init() {

  static ClassDocumentation<UEDF1F1P0Vertex> documentation
    ("This class couples a pair of level-1 KK fermions to an SM "
     "photon.");

}


void UEDF1F1P0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr ,
#ifndef NDEBUG
				  tcPDPtr part3) {
#else
				  tcPDPtr ) {
#endif
  long iferm = abs(part1->id());
  assert(part3->id()==ParticleID::gamma);
  assert((iferm >= 5100001 && iferm <= 5100006) || 
	 (iferm >= 5100011 && iferm <= 5100016) ||
	 (iferm >= 6100001 && iferm <= 6100006) || 
	 (iferm >= 6100011 && iferm <= 6100016));
  if(q2 != q2Last_ || coupLast_ == 0. ) {
    q2Last_ = q2;
    coupLast_ = -electroMagneticCoupling(q2);
  }
  norm(coupLast_);
  if(iferm != fermLast_) {
    fermLast_ = iferm;
    int smtype = (iferm > 6000000) ? iferm - 6100000 : iferm - 5100000;
    if(smtype >= 11) 
      LRLast_ = charges_[0];
    else
      LRLast_ = (smtype % 2 == 0) ? charges_[2] : charges_[1];
  }
  left(LRLast_);
  right(LRLast_);
}
#line 1 "./UEDF1F1W0Vertex.cc"
// -*- C++ -*-
//
// UEDF1F1W0Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F1W0Vertex class.
//

#include "UEDF1F1W0Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F1W0Vertex::UEDF1F1W0Vertex(): includeMixing_(true),
				    theRadius(ZERO), theQ2Last(ZERO), 
				    theCoupLast(0.), 
				    thefermALast(0), thefermBLast(0) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

void UEDF1F1W0Vertex::doinit() {
  //outgoing W+
  for( long i = 2; i < 17; i += 2 ) {
    if( i == 7 ) i += 5;
    addToList(-5100000 - i, 5100000 + i - 1, 24);
    if( i < 7 ) {
      addToList(-6100000 - i, 6100000 + i - 1, 24);
    }
  }
  if(includeMixing_) {
    addToList(-6100006, 5100005, 24);
    addToList(-5100006, 6100005, 24);
  }
  //outgoing W-
  for( long i = 1; i < 16; i += 2 ) {
    if( i == 6 ) i += 5;
    addToList(-5100000 - i, 5100001 + i, -24);
    if( i < 6 ) {
      addToList(-6100000 - i, 6100001 + i, -24);
    }
  }
  if(includeMixing_) {
    addToList(-6100005, 5100006, -24);
    addToList(-5100005, 6100006, -24);
  }
  FFVVertex::doinit();
  tUEDBasePtr model = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "UEDF1F1W0Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theRadius = model->compactRadius();
}

void UEDF1F1W0Vertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theRadius,1/GeV) << includeMixing_;
}

void UEDF1F1W0Vertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theRadius,1/GeV) >> includeMixing_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDF1F1W0Vertex,FFVVertex>
describeHerwigUEDF1F1W0Vertex("Herwig::UEDF1F1W0Vertex", "HwUED.so");

void UEDF1F1W0Vertex::Init() {

  static ClassDocumentation<UEDF1F1W0Vertex> documentation
    ("This class implements the coupling of a pair of level-1 KK fermions"
     "to an SM W boson");

  static Switch<UEDF1F1W0Vertex,bool> interfaceIncludeMixing
    ("IncludeMixing",
     "Include the mixing",
     &UEDF1F1W0Vertex::includeMixing_, true, false, false);
  static SwitchOption interfaceIncludeMixingYes
    (interfaceIncludeMixing,
     "Yes",
     "Include mixing",
     true);
  static SwitchOption interfaceIncludeMixingNo
    (interfaceIncludeMixing,
     "No",
     "Don't include mixing",
     false);

}

void UEDF1F1W0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
#ifndef NDEBUG
				  tcPDPtr part3) {
#else
				  tcPDPtr) {
#endif
  long ianti(abs(part1->id())), iferm(abs(part2->id()));
  assert( abs(part3->id()) == 24 );
  bool ferma = (iferm >= 5100001 && iferm <= 5100006) ||
    (iferm >= 6100001 && iferm <= 6100006) || 
    (iferm >= 5100011 && iferm <= 5100016) ||
    (iferm >= 6100011 && iferm <= 6100016); 
  bool fermb = (ianti >= 5100001 && ianti <= 5100006) ||
    (ianti >= 6100001 && ianti <= 6100006) || 
    (ianti >= 5100011 && ianti <= 5100016) ||
    (ianti >= 6100011 && ianti <= 6100016);
  if( !ferma || !fermb ) 
    throw HelicityLogicalError() << "UEDF1F1W0Vertex::setCoupling - "
				 << "There is an unknown particle(s) in the "
				 << "UED F^(1) F^(1) W^(0) vertex. ID: " 
				 << ianti << " " << iferm 
				 << Exception::runerror;
  if(q2 != theQ2Last || theCoupLast == 0. ) {
    theQ2Last = q2;
    theCoupLast = sqrt(0.5)*weakCoupling(q2);
  }
  if(iferm != thefermALast || ianti != thefermBLast) {
    thefermALast = iferm;
    thefermBLast = ianti;
    int stateA(ianti/1000000), stateB(iferm/1000000);
    long sma = (stateA == 6) ? ianti - 6100000 : ianti - 5100000;
    long smb = (stateB == 6) ? iferm - 6100000 : iferm - 5100000;
    double afu(0.), afd(0.);
    if(includeMixing_) {
      if( sma % 2 == 0 ) {
	afu = atan(getParticleData(sma)->mass()*theRadius)/2.;
	afd = atan(getParticleData(smb)->mass()*theRadius)/2.;
      }
      else {
	afd = atan(getParticleData(sma)->mass()*theRadius)/2.;
	afu = atan(getParticleData(smb)->mass()*theRadius)/2.;
      }
    }
    else {
      afd = afu = 0.;
    }
    if( stateA == stateB ) {
      if( stateA == 5 ) {
	left(cos(afu)*cos(afd));
	right(cos(afu)*cos(afd));
      }
      else {
	left(sin(afu)*sin(afd));
	right(sin(afu)*sin(afd));
      }
    }
    else {
      if( sma % 2 == 0 ) {
	if( stateA == 5 ) {
	  left(cos(afu)*sin(afd));
	  right(-cos(afu)*sin(afd));
	}
	else {
	  left(sin(afu)*cos(afd));
	  right(-sin(afu)*cos(afd));
	}
      }
      else {
	if( stateA == 5 ) {
	  left(sin(afu)*cos(afd));
	  right(-sin(afu)*cos(afd));
	}
	else {
	  left(cos(afu)*sin(afd));
	  right(-cos(afu)*sin(afd));
	}
      }
    }
  }
  norm(theCoupLast);
}
#line 1 "./UEDF1F0W1Vertex.cc"
// -*- C++ -*-
//
// UEDF1F0W1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0W1Vertex class.
//

#include "UEDF1F0W1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F0W1Vertex::UEDF1F0W1Vertex() : theSinW(0.), theCosW(0.), theSinOne(0.),
				     theCosOne(0.), theSinWmO(0.), 
				     theCosWmO(0.), 
				     theCKM(0, vector<Complex>(0, 0.)),
				     theq2last(),
				     theCouplast(0.), theLlast(0.),
				     theRlast(0.), theGBlast(0), 
				     theKKlast(0), theSMlast(0) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

void UEDF1F0W1Vertex::doinit() {
  //outgoing W+
  for(long i = 2; i < 7; i += 2) {
    for(long j = 1; j < 6; j += 2) {
      addToList( -i, 5100000 + j, 5100024 );
      addToList( -(5100000 + i), j, 5100024 );
    }
  }
  for(long i = 11; i < 17; i += 2) {
    addToList( -i-1, 5100000 + i, 5100024 );
    addToList( -(5100001 + i), i, 5100024 );
  }
  //outgoing W-
  for(long i = 1; i < 6; i += 2) {
    for(long j = 2 ; j < 7; j += 2) {
      addToList( -i, 5100000 + j, -5100024 );
      addToList( -(5100000 + i), j, -5100024 );
    }
  }
  for(long i = 11; i < 17; i += 2) {
    addToList( -i, 5100001 + i, -5100024 );
    addToList(-(5100000 + i), i + 1, -5100024);
  }
  long boson[2] = {5100022,5100023}; 
  for(long b = 0; b < 2; ++b) { 
    //QQ
    for(int i = 1; i < 7; ++i) {
      addToList( -i, i + 5100000, boson[b]);
      addToList(-(i + 5100000), i, boson[b]);

      addToList(-i, i + 6100000, boson[b]);
      addToList(-(i + 6100000), i, boson[b]);
    }
    //LL
    for(int i = 11; i < 17; ++i) {
      addToList( -i, i + 5100000, boson[b]);
      addToList(-(i + 5100000), i, boson[b]);
      if( i % 2 != 0 ) {
	addToList(-i, i + 6100000, boson[b]);
	addToList(-(i + 6100000), i, boson[b]);
      }
    }
  }
  FFVVertex::doinit();
  tUEDBasePtr UEDBase = 
    dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!UEDBase)
    throw InitException() << "UEDF1F0W1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theSinW = sqrt(sin2ThetaW());
  theCosW = sqrt( 1. - sqr(theSinW));
  theSinOne = UEDBase->sinThetaOne();
  theCosOne = sqrt(1. - sqr(theSinOne)); 
  theSinWmO = theSinW*theCosOne - theSinOne*theCosW;
  theCosWmO = theCosW*theCosOne + theSinW*theSinOne;
  theCKM = dynamic_ptr_cast<Ptr<StandardCKM>::transient_pointer>
    (UEDBase->CKM())->getUnsquaredMatrix(UEDBase->families());
}


void UEDF1F0W1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theSinW << theCosW << theSinOne << theCosOne
     << theSinWmO << theCosWmO << theCKM;
}

void UEDF1F0W1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theSinW >> theCosW >> theSinOne >> theCosOne
     >> theSinWmO >> theCosWmO >> theCKM;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDF1F0W1Vertex,FFVVertex>
describeHerwigUEDF1F0W1Vertex("Herwig::UEDF1F0W1Vertex", "HwUED.so");

void UEDF1F0W1Vertex::Init() {

  static ClassDocumentation<UEDF1F0W1Vertex> documentation
    ("This is the coupling of a KK1 W boson to a KK1 fermion and "
     "a SM fermion.");

}

void UEDF1F0W1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long id1(abs(part1->id())), id2(abs(part2->id())),
    gboson(abs(part3->id())), kkparticle(0), smID(0);
  assert( gboson == 5100022  || gboson == 5100023 || gboson == 5100024 );
  if( id1 > 5000000 ) {
    kkparticle = id1;
    smID = id2;
  }
  else {
    kkparticle = id2;
    smID = id1;
  }
  if( (kkparticle >= 5100001 && kkparticle <= 5100006) ||
      (kkparticle >= 6100001 && kkparticle <= 6100006) ||
      (kkparticle >= 5100011 && kkparticle <= 5100016) ||
      (kkparticle >= 6100011 && kkparticle <= 6100016) ) {
    if(q2 != theq2last || theCouplast == 0.) {
      theq2last = q2;
      theCouplast = electroMagneticCoupling(q2);
    }
    if( gboson != theGBlast || kkparticle != theKKlast || smID != theSMlast ) {
      theGBlast = gboson;
      theKKlast = kkparticle;
      theSMlast = smID;
      if( gboson == 5100024 ) {
	Complex ckm(1.);
	if( smID >= 1 && smID <= 6 ) {
	  long smIDb(kkparticle - 5100000);
	  if( smID % 2 != 0 ) swap(smID, smIDb);
	  ckm = theCKM[smID/2 - 1][(smIDb - 1)/2];
	}
	theLlast = -ckm/sqrt(2)/theSinW;
	theRlast = 0.;
      }
      else if( gboson == 5100022 || gboson == 5100023 ) {
	double Qf = getParticleData(smID)->charge()/eplus;
	if( kkparticle/1000000 == 5 ) {
	  theRlast = 0.;
	  double I3f = (abs(smID) % 2 == 0) ? 0.5 : -0.5;
	  if( gboson == 5100023 )
	    theLlast = (Qf*theSinOne 
			- I3f*theCosWmO/theSinW)/theCosW;
	  else
	    theLlast = -(Qf*theCosOne 
			 - I3f*theSinWmO/theSinW)/theCosW;
	}
	else {
	  theLlast = 0.;
	  if( gboson == 5100023 )
	    theRlast = Qf*theSinOne/theCosW;
	  else
	    theRlast = -Qf*theCosOne/theCosW;
	}
      }
    }
    norm(theCouplast);
    left(theLlast);
    right(theRlast);
  }
  else
    throw HelicityLogicalError() << "UEDF1F0W1Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex! "
				 << kkparticle
				 << Exception::warning;
}

#line 1 "./UEDP0H1H1Vertex.cc"
// -*- C++ -*-
//
// UEDP0H1H1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDP0H1H1Vertex class.
//

#include "UEDP0H1H1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

UEDP0H1H1Vertex::UEDP0H1H1Vertex() : theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void UEDP0H1H1Vertex::doinit() {
  addToList(22, 5100037, -5100037);
  VSSVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<UEDP0H1H1Vertex,VSSVertex>
describeHerwigUEDP0H1H1Vertex("Herwig::UEDP0H1H1Vertex", "HwUED.so");

void UEDP0H1H1Vertex::Init() {

  static ClassDocumentation<UEDP0H1H1Vertex> documentation
    ("This is the coupling of the SM photon to the level-1 charged higgs.");

}
#ifndef NDEBUG
void UEDP0H1H1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr ) {
#else
void UEDP0H1H1Vertex::setCoupling(Energy2 q2, tcPDPtr , tcPDPtr part2,
				  tcPDPtr ) {
#endif

  assert(part1->id()==ParticleID::gamma);
  assert(abs(part2->id()) == 5100037);
  if(q2 != theq2Last || theCoupLast == 0.) {
    theq2Last = q2;
    theCoupLast = electroMagneticCoupling(q2);
  }
  if(part2->id()>0) 
    norm(-theCoupLast);
  else
    norm( theCoupLast);
}
#line 1 "./UEDZ0H1H1Vertex.cc"
// -*- C++ -*-
//
// UEDZ0H1H1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDZ0H1H1Vertex class.
//

#include "UEDZ0H1H1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDZ0H1H1Vertex::UEDZ0H1H1Vertex() : theCosThetaW(0.), theCosTheta2W(0.), theMw2(), 
				     theR2(), theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void UEDZ0H1H1Vertex::doinit() {
  addToList(23, 5100037, -5100037);
  VSSVertex::doinit();
  tUEDBasePtr theUEDBase = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!theUEDBase)
    throw InitException() << "UEDZ0H1H1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theCosThetaW = sqrt(1. - sin2ThetaW());
  theCosTheta2W = 1. - 2.*sin2ThetaW();
  theMw2 = sqr(getParticleData(24)->mass());
  theR2 = sqr(theUEDBase->compactRadius());
}

void UEDZ0H1H1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theCosThetaW << theCosTheta2W
     << ounit(theMw2,GeV2) << ounit(theR2,1/GeV2);
}

void UEDZ0H1H1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theCosThetaW >> theCosTheta2W
     >> iunit(theMw2,GeV2) >> iunit(theR2,1/GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDZ0H1H1Vertex,VSSVertex>
describeHerwigUEDZ0H1H1Vertex("Herwig::UEDZ0H1H1Vertex", "HwUED.so");

void UEDZ0H1H1Vertex::Init() {

  static ClassDocumentation<UEDZ0H1H1Vertex> documentation
    ("This is the coupling of the SM Z boson to the level-1 charged higgs.");

}

void UEDZ0H1H1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long kkhiggs(0);
  if(part1->id() == ParticleID::Z0)
    kkhiggs = abs(part2->id());
  else if(part2->id() == ParticleID::Z0)
    kkhiggs = abs(part1->id());
  else if(part3->id() == ParticleID::Z0)
    kkhiggs = abs(part1->id());
  else {
    throw HelicityLogicalError() << "UEDZ0H1H1Vertex::setCoupling - There is no "
				 << "SM photon in this vertex!." 
				 << Exception::warning;
    return;
  }
  if(kkhiggs == 5100037) {
    if(q2 != theq2Last || theCoupLast == 0.) {
      theq2Last = q2;
      theCoupLast = 
	Complex(0., 1.)*weakCoupling(q2);
      theCoupLast *= ( (theCosTheta2W/2./theCosThetaW) 
		       - sqr(theCosThetaW)*theMw2*theR2 )/(1. + theMw2*theR2);
    }
    norm(theCoupLast);
  }
  else
    throw HelicityLogicalError() << "UEDZ0H1H1Vertex::setCoupling - There is no "
				 << "level-1 higgs in this vertex! " << kkhiggs
				 << Exception::warning;
}
#line 1 "./UEDW0A1H1Vertex.cc"
// -*- C++ -*-
//
// UEDW0A1H1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDW0A1H1Vertex class.
//

#include "UEDW0A1H1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

void UEDW0A1H1Vertex::doinit() {
  addToList( 24, 5100036, -5100037);
  addToList(-24, 5100036,  5100037);
  VSSVertex::doinit();
  tUEDBasePtr UEDBase = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!UEDBase) throw InitException() 
    << "UEDW0A1H1Vertex::doinit() - The pointer to "
    << "the UEDBase object is null!" << Exception::runerror;
  theMw2 = sqr(getParticleData(24)->mass());
  theMz2 = sqr(getParticleData(23)->mass());
  theR2 = sqr(UEDBase->compactRadius());
}

UEDW0A1H1Vertex::UEDW0A1H1Vertex() : theMw2(), theMz2(), theR2(), 
				     theq2Last(), theCoupLast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void UEDW0A1H1Vertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMw2,GeV2) << ounit(theMz2,GeV2) << ounit(theR2,1/GeV2);  
}

void UEDW0A1H1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMw2,GeV2) >> iunit(theMz2,GeV2) >> iunit(theR2,1/GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDW0A1H1Vertex,VSSVertex>
describeHerwigUEDW0A1H1Vertex("Herwig::UEDW0A1H1Vertex", "HwUED.so");

void UEDW0A1H1Vertex::Init() {

  static ClassDocumentation<UEDW0A1H1Vertex> documentation
    ("The coupling of a SM W boson to a level-1 charged higgs and the "
     "level-1 heavy neutral higgs");

}

void UEDW0A1H1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long chiggs(0);
  if(abs(part1->id()) == ParticleID::Wplus)
    chiggs = (abs(part2->id()) == 5100037) ? part2->id() : part3->id();
  else if(abs(part2->id()) == ParticleID::Wplus)
    chiggs = (abs(part1->id()) == 5100037) ? part1->id() : part3->id();
  else if(abs(part3->id()) == ParticleID::Wplus)
    chiggs = (abs(part1->id()) == 5100037) ? part1->id() : part2->id();
  else {
    throw HelicityLogicalError() << "UEDW0A1H1Vertex::setCoupling - "
				 << "There is no SM W boson in this vertex"
				 << Exception::warning;
    return;
  }
  if(abs(chiggs) == 5100037) {
    if(q2 != theq2Last || theCoupLast == 0.) {
      theq2Last = q2;
      theCoupLast = weakCoupling(q2);
      double mwRs = theMw2*theR2;
      double denom = sqrt( (1 + mwRs)*(1. + theMw2*theR2) );
      theCoupLast *= ( 0.5 + mwRs )/denom;
    }
    if(chiggs > 0) theCoupLast *= -1.;
    norm(theCoupLast);
  }
  else
    throw HelicityLogicalError() << "UEDW0A1H1Vertex::setCoupling - "
				 << "There is an unknown particle in this " 
				 << "vertex " << chiggs
				 << Exception::runerror;
}

#line 1 "./UEDZ0A1h1Vertex.cc"
// -*- C++ -*-
//
// UEDZ0A1h1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDZ0A1h1Vertex class.
//

#include "UEDZ0A1h1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDZ0A1h1Vertex::UEDZ0A1h1Vertex() : theSin2ThetaW(0.), theKappa(0.),	    
				     theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void UEDZ0A1h1Vertex::doinit() {
  addToList(23, 5100036, 5100025);
  VSSVertex::doinit();
  tUEDBasePtr UEDBase = 
    dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!UEDBase)
    throw InitException() << "UEDZ0A1h1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  double sw2 = sin2ThetaW();
  theSin2ThetaW = 2.*sqrt(sw2*(1. - sw2));
  Energy2 mz2 = sqr(getParticleData(23)->mass());
  InvEnergy2 rad2 = sqr(UEDBase->compactRadius());
  theKappa = 1./sqrt(1. + mz2*rad2);
}

void UEDZ0A1h1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theSin2ThetaW << theKappa;
}

void UEDZ0A1h1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theSin2ThetaW >> theKappa;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDZ0A1h1Vertex,VSSVertex>
describeHerwigUEDZ0A1h1Vertex("Herwig::UEDZ0A1h1Vertex", "HwUED.so");

void UEDZ0A1h1Vertex::Init() {

  static ClassDocumentation<UEDZ0A1h1Vertex> documentation
    ("The coupling of an SM Z boson to a level-1 CP-Odd pseudo-scalar "
     "and level 1 higgs.");

}

void UEDZ0A1h1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2, 
				  tcPDPtr part3) {
  long scaA(0), scaB(0);
  if(part1->id() == ParticleID::Z0) {
    scaA = part2->id();
    scaB = part3->id();
  }
  else if(part2->id() == ParticleID::Z0) {
    scaA = part1->id();
    scaB = part3->id();
  }
  else if(part3->id() == ParticleID::Z0) {
    scaA =  part1->id();
    scaB = part2->id();
  }
  else {
    throw HelicityLogicalError() << "UEDZ0A1h1Vertex::setCoupling - "
				 << "There is no SM Z boson in this vertex"
				 << Exception::warning;
  }
  if( (scaA == 5100036 && scaB == 5100025) ||
      (scaB == 5100036 && scaA == 5100025) ) {
    if(q2 != theq2Last || theCoupLast == 0.) {
      theq2Last = q2;
      theCoupLast = theKappa*electroMagneticCoupling(q2)/theSin2ThetaW;
    }
    norm(theCoupLast); 
  }
  else
    throw HelicityLogicalError() << "UEDZ0A1h1Vertex::setCoupling - "
				 << "There is an unknown particle in this "
				 << "vertex. " << scaA << " " << scaB
				 << Exception::warning;
}
#line 1 "./UEDW0W1W1Vertex.cc"
// -*- C++ -*-
//
// UEDW0W1W1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDW0W1W1Vertex class.
//

#include "UEDW0W1W1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDW0W1W1Vertex::UEDW0W1W1Vertex() : theSinW(0.), theCosW(0.),
				     theSinThetaOne(0.), theCosThetaOne(0.),
				     theq2last(), theElast(0.), theCouplast(0.),
				     theSMlast(0), theKKlast(0) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void UEDW0W1W1Vertex::doinit() {
  addToList( 22, -5100024, 5100024);
  addToList( 23, -5100024, 5100024);

  addToList( 24, -5100024, 5100022);
  addToList( 24, -5100024, 5100023);

  addToList(-24,  5100024, 5100022);
  addToList(-24,  5100024, 5100023);
  VVVVertex::doinit();
  tUEDBasePtr model = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "UEDW0W1W1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theSinW = sqrt(sin2ThetaW());
  theCosW = sqrt( 1. - sqr(theSinW) );
  theSinThetaOne = model->sinThetaOne();
  theCosThetaOne = sqrt( 1. - sqr(theSinThetaOne));
}

void UEDW0W1W1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theSinW << theCosW << theSinThetaOne << theCosThetaOne;
}

void UEDW0W1W1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theSinW >> theCosW >> theSinThetaOne >> theCosThetaOne;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDW0W1W1Vertex,VVVVertex>
describeHerwigUEDW0W1W1Vertex("Herwig::UEDW0W1W1Vertex", "HwUED.so");

void UEDW0W1W1Vertex::Init() {

  static ClassDocumentation<UEDW0W1W1Vertex> documentation
    ("The coupling of an SM W boson to a level 1 KK W and KK Z and KK photon");

}

/// \todo look again
void UEDW0W1W1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long id1(abs(part1->id())), id2(abs(part2->id())), id3(abs(part3->id())), 
    smID(0), kkparticle(0);
  double perm(1.);
  if( id1 == 22 || id1 == 23) {
    smID = id1;
    kkparticle = id2;
    if(part2->id()>0) perm=-1.;
  }
  else if(id2 == 22 || id2 == 23) {
    smID = id2;
    kkparticle = id1;
    if(part1->id()<0) perm=-1.;
  }
  else if(id3 == 22 || id3 == 23) {
    smID = id3;
    kkparticle = id1;
    if(part1->id()>0) perm=-1.;
  }
  else if(id1 == 24 ) {
    if( part1->id() == 24 ) perm = -1.;
    smID = id1;
    kkparticle = (id2 == 5100024) ? id3 : id2;
    if(id3 == 5100024) perm *=-1.;
  }
  else if( id2 == 24 ) {
    if( part2->id() == 24 ) perm = -1.;
    smID = id2;
    kkparticle = (id1 == 5100024) ? id3 : id1;
    if(id1 == 5100024) perm *=-1.;
  }
  else if( id3 == 24 ) {
    if( part3->id() == 24 ) perm = -1.;
    smID = id3;
    kkparticle = (id1 == 5100024) ? id2 : id1;
    if(id2 == 5100024) perm *=-1.;
  }
  else {
    throw HelicityLogicalError()
      << "UEDW0W1W1Vertex::setCoupling() - There is no SM gauge boson in "
      << "this vertex. " << id1 << " " << id2 << " " << id3 
      << Exception::warning; 
    norm(0.);
    return;
  }
  if( q2 != theq2last || theElast == 0.) {
    theq2last = q2;
    theElast = electroMagneticCoupling(q2);
  }
  
  if( smID != theSMlast || kkparticle != theKKlast ) { 
    theSMlast = smID;
    theKKlast = kkparticle;
    if( smID == 22 )
      theCouplast = 1.;
    else if(smID == 23) 
      theCouplast = theCosW/theSinW;
    else {
      if( kkparticle == 5100023 )
	theCouplast = theCosThetaOne/theSinW;
      else
	theCouplast = theSinThetaOne/theSinW;
    }
  }
  norm(perm*theElast*theCouplast);
}
#line 1 "./UEDF1F0H1Vertex.cc"
// -*- C++ -*-
//
// UEDF1F0H1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0H1Vertex class.
//

#include "UEDF1F0H1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F0H1Vertex::UEDF1F0H1Vertex() : theRadius(ZERO), theMw(ZERO), 
				     theSinThetaW(0.), theq2Last(ZERO),
				     theCoupLast(0.), theLeftLast(0.),
				     theRightLast(0.), theAntiLast(0),
				     theFermLast(0), theHLast(0) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

void UEDF1F0H1Vertex::doinit() {
  long heavy[3] = {5, 6, 15};
  //h0
  for( int i = 0; i < 3; ++i ) {
    addToList(-5100000 - i, 5100000 + i, 25);
    addToList(-6100000 - i, 6100000 + i, 25);
    addToList(-5100000 - i, 6100000 + i, 25);
    addToList(-6100000 - i, 5100000 + i, 25);
  }
  // Neutral KK-Higgs
  long higgs[2] = {5100025, 5100036};
  for( unsigned int h = 0; h < 2; ++h ) {
    for( int i = 0; i < 3; ++i ) {
      addToList(-heavy[i], 5100000 + heavy[i], higgs[h]);
      addToList(-5100000 - heavy[i], heavy[i], higgs[h]);
      addToList(-heavy[i], 6100000 + heavy[i], higgs[h]);
      addToList(-6100000 - heavy[i], heavy[i], higgs[h]);
    }
  }

  //KK-charged higgs
  //outgoing H+
  addToList(-5100006, 5, 5100037);
  addToList(-6100006, 5, 5100037);

  addToList(-6, 5100005, 5100037);
  addToList(-6, 6100005, 5100037);

  addToList(-5100016, 15, 5100037);
  addToList(-6100016, 15, 5100037);

  addToList(-16, 5100015, 5100037);
  addToList(-16, 6100015, 5100037);

  //outgoing H-
  addToList(-5100005, 6,-5100037);
  addToList(-6100005, 6,-5100037);

  addToList(-5, 5100006,-5100037);
  addToList(-5, 6100006,-5100037);

  addToList(-5100015, 16,-5100037);
  addToList(-6100015, 16,-5100037);

  addToList(-15, 5100016,-5100037);
  addToList(-15, 6100016,-5100037);
  FFSVertex::doinit();
  tUEDBasePtr UEDBase = 
    dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!UEDBase)
    throw InitException() << "UEDF1F0H1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theRadius = UEDBase->compactRadius();
  theSinThetaW = sqrt(sin2ThetaW());
  theCosThetaW = sqrt(1. - sin2ThetaW());
  theMw = getParticleData(24)->mass();
  theMz = getParticleData(23)->mass();
}


void UEDF1F0H1Vertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theRadius,1/GeV) << ounit(theMw,GeV) << theSinThetaW 
     << ounit(theMz, GeV) << theCosThetaW;
}

void UEDF1F0H1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theRadius,1/GeV) >> iunit(theMw,GeV) >> theSinThetaW
     >> iunit(theMz, GeV) >> theCosThetaW;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDF1F0H1Vertex,FFSVertex>
describeHerwigUEDF1F0H1Vertex("Herwig::UEDF1F0H1Vertex", "HwUED.so");

void UEDF1F0H1Vertex::Init() {

  static ClassDocumentation<UEDF1F0H1Vertex> documentation
    ("The coupling involving a KK-Higgs and a pair of fermions.");

}

void UEDF1F0H1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long anti(abs(part1->id())), ferm(abs(part2->id())), higgs(part3->id());
  if( ferm > 17 ) swap( ferm, anti);

  if( anti != theAntiLast || ferm != theFermLast || higgs != theHLast ) { 
    theAntiLast = anti;
    theFermLast = ferm;
    theHLast = higgs;

    tcPDPtr pd;
    if( higgs != 25 ) {
      pd = getParticleData(ferm);
    }
    else {
      long smid = ( ferm/1000000 == 5 ) ? ferm - 5100000 : ferm - 6100000;
      pd = getParticleData(smid);
    }
    Energy mf = pd->mass();
    double alpha = mf*theRadius/2.;
    double salpha = sin(alpha);
    double calpha = cos(alpha);
    double fact(0.);
    if( abs(higgs) == 5100037 ) {
      fact = theRadius/2./sqrt(1. + sqr(theMw*theRadius)) * UnitRemoval::E;
      theRightLast = theRadius*theMw;
      if(anti/1000000 == 5) {
	Energy mfk  = getParticleData(anti - 5100000)->mass();
	theLeftLast = (theMw*calpha - (mfk*salpha/theRadius/theMw)) 
	  * UnitRemoval::InvE;
      
	theRightLast *= calpha*mfk* UnitRemoval::InvE;
      }
      else {
	Energy mfk = getParticleData(anti - 6100000)->mass();
	theLeftLast = (theMw*salpha +(mfk*calpha/theRadius/theMw))
	  * UnitRemoval::InvE;
	theRightLast *= -salpha*mfk*UnitRemoval::InvE;
      }
      theLeftLast *= fact;
      theRightLast *= fact;
      if( higgs < 0 ) swap( theLeftLast, theRightLast );
    }
    else if( higgs == 5100025 ) {
      fact = mf/theMw/2.;
      if( anti/1000000 == 5 )
	theLeftLast = salpha + calpha;
      else 
	theLeftLast = salpha - calpha;
      theRightLast = theLeftLast;

      theLeftLast *= fact;
      theRightLast *= fact;
    }
    else if( higgs == 5100036 ) {
      fact = theRadius/theCosThetaW/sqrt(1.+sqr(theMz*theRadius))*UnitRemoval::E;
    
      double i3f = ( ferm % 2 == 0 ) ? 0.5 : -0.5;
      double qf = pd->charge()/eplus;
      if( anti/1000000 == 5 ) {
	theLeftLast = (theMz*calpha*(i3f - qf*sqr(theSinThetaW))
		       - mf*salpha/theRadius/theMw) * UnitRemoval::InvE;
	theRightLast = (-theMz*salpha*qf*sqr(theSinThetaW) 
			+ mf*calpha/theRadius/theMw) * UnitRemoval::InvE; 
      }
      else {
	theLeftLast = (theMz*salpha*(i3f - qf*sqr(theSinThetaW))
		       - mf*calpha/theRadius/theMw) * UnitRemoval::InvE;
	theRightLast = (-theMz*calpha*qf*sqr(theSinThetaW) 
			+ mf*salpha/theRadius/theMw)*UnitRemoval::InvE; 
      }
      theLeftLast *= fact;
      theRightLast *= fact;
    }
    else {
      theLeftLast = mf*calpha*salpha/2./theMw;
      if( ferm/1000000 == 5 ) theLeftLast *= -1.;  
      theRightLast = theLeftLast;
    }
  }


  if(q2 != theq2Last || theCoupLast == 0.) {
    theq2Last = q2;
    theCoupLast = weakCoupling(q2);
  }

  norm(theCoupLast);
  left(theLeftLast);
  right(theRightLast);
}
