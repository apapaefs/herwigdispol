#line 1 "./RPV.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPV class.
//

#include "RPV.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RPV::persistentOutput(PersistentOStream & os ) const {
  os << lambdaLLE_ << lambdaLQD_ << lambdaUDD_ << ounit(vnu_,GeV)
     << upSquarkMix_ << downSquarkMix_ << triLinearOnly_
     << LLEVertex_ << LQDVertex_ << UDDVertex_ 
     << ounit(epsilon_,GeV) << ounit(epsB_,GeV);
}

void RPV::persistentInput(PersistentIStream & is, int) {
  is >> lambdaLLE_ >> lambdaLQD_ >> lambdaUDD_ >> iunit(vnu_,GeV)
     >> upSquarkMix_ >> downSquarkMix_ >> triLinearOnly_
     >> LLEVertex_ >> LQDVertex_ >> UDDVertex_ 
     >> iunit(epsilon_,GeV) >> iunit(epsB_,GeV);
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<RPV,MSSM>
describeHerwigRPV("Herwig::RPV", "HwSusy.so HwRPV.so");

void RPV::Init() {

  static ClassDocumentation<RPV> documentation
    ("The RPV class is the base class for the implementation of the"
     " R-parity violating MSSM.");

  static Reference<RPV,AbstractFFSVertex> interfaceLLEVertex
    ("Vertex/LLE",
     "The vertex for the trillinear LLE interaction",
     &RPV::LLEVertex_, false, false, true, false, false);

  static Reference<RPV,AbstractFFSVertex> interfaceLQDVertex
    ("Vertex/LQD",
     "The vertex for the trillinear LQD interaction",
     &RPV::LQDVertex_, false, false, true, false, false);

  static Reference<RPV,AbstractFFSVertex> interfaceUDDVertex
    ("Vertex/UDD",
     "The vertex for the trillinear UDD interaction",
     &RPV::UDDVertex_, false, false, true, false, false);

  static Switch<RPV,bool> interfaceTriLinearOnly
    ("TriLinearOnly",
     "Only include trilinears and take rest of model to be MSSM",
     &RPV::triLinearOnly_, false, false, false);
  static SwitchOption interfaceTriLinearOnlyYes
    (interfaceTriLinearOnly,
     "Yes",
     "Trilinears + MSSM",
     true);
  static SwitchOption interfaceTriLinearOnlyNo
    (interfaceTriLinearOnly,
     "No",
     "All RPV couplings and mixings",
     false);

}

void RPV::extractParameters(bool checkmodel) {
  MSSM::extractParameters(false);
  if(checkmodel) {
    map<string,ParamMap>::const_iterator pit;
    pit = parameters().find("modsel");
    if(pit == parameters().end()) return;
    ParamMap::const_iterator it;
    // nmssm or mssm
    it = pit->second.find(3);
    int inmssm = (it != pit->second.end()) ? int(it->second) : 0;
    if(inmssm != 0) 
      throw Exception() << "R-parity violating MSSM model"
			<< " used but NMSSM read in." << Exception::runerror; 
    // RPV
    it = pit->second.find(4);
    int irpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(irpv != 1) throw Exception() << "RPV model used but no RPV in input file"
				    << Exception::runerror; 
    // CPV
    it = pit->second.find(5);
    int icpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(icpv != 0) throw Exception() << "RPV model does not support CPV" 
				  << Exception::runerror; 
    // flavour violation
    it = pit->second.find(6);
    int ifv = (it != pit->second.end()) ? int(it->second) : 0;
    if(ifv != 0) throw Exception() << "RPV model does not support "
				 << "flavour violation"
				 << Exception::runerror;
  }
  // get the RPV parameters
  // lambda
  map<string,ParamMap>::const_iterator pit;
  pit=parameters().find("rvlamlle");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaLLE_[i][j][k] = it->second;
    }
  }
  // lambda'
  pit=parameters().find("rvlamlqd");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaLQD_[i][j][k] = it->second;
    }
  }
  // lambda''
  pit=parameters().find("rvlamudd");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaUDD_[i][j][k] = it->second;
    }
  }
  // sneutrino vevs
  pit=parameters().find("rvsnvev");
  vnu_.resize(3);
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first>0) {
	assert(it->first>=1&&it->first<=3);
	vnu_[it->first-1] = it->second*GeV;
      }
    }
  }
  // bilinears
  pit=parameters().find("rvkappa");
  epsilon_.resize(3);
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first>0) {
	assert(it->first>=1&&it->first<=3);
	epsilon_[it->first-1] = it->second*GeV;
      }
    }
  }
  // blinear soft terms
  pit=parameters().find("rvd");
  epsB_.resize(3);
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first>0) {
	assert(it->first>=1&&it->first<=3);
	epsB_[it->first-1] = it->second*GeV;
      }
    }
  }
}

void RPV::createMixingMatrices() {
  map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    // pseudo-scalar higgs mixing
    if (name == "rvamix") {
      MixingMatrixPtr temp;
      createMixingMatrix(temp,name,it->second.second,it->second.first);
      CPoddHiggsMix(temp);
    }
    else if (name == "rvlmix") {
      MixingMatrixPtr temp;
      createMixingMatrix(temp,name,it->second.second,it->second.first);
      ChargedHiggsMix(temp);
    }
    else if (name == "dsqmix" ) {
      createMixingMatrix(downSquarkMix_,name,it->second.second,it->second.first);
    }
    else if (name == "usqmix" ) {
      createMixingMatrix(upSquarkMix_,name,it->second.second,it->second.first);
    }
  }
  // base class for neutralinos and charginos
  MSSM::createMixingMatrices();
  // now adjust the mixing matrices to have our structure
  // first bloodly SPHENO as it doesn't obey the SLHA  
  map<string,StringMap>::const_iterator sit = info().find("spinfo");
  string program;
  if(sit!=info().end()) {
    StringMap::const_iterator pit = sit->second.find(1);
    if(pit!=sit->second.end()) program = pit->second;
  }
  if(program=="SPheno") {
    map<string,ParamMap>::const_iterator fit=parameters().find("mass");
    if(fit==parameters().end()) 
      throw Exception() << "BLOCK MASS not found in input file"
			<< " can't set masses of SUSY particles"
			<< Exception::runerror;
    // adjust the charged scalars
    map<double,long> massMap;
    massMap[findValue(fit,     37,"mass",     "37")] =      37;
    massMap[findValue(fit,1000011,"mass","1000011")] = 1000011;
    massMap[findValue(fit,1000013,"mass","1000013")] = 1000013;
    massMap[findValue(fit,1000015,"mass","1000015")] = 1000015;
    massMap[findValue(fit,2000011,"mass","2000011")] = 2000011;
    massMap[findValue(fit,2000013,"mass","2000013")] = 2000013;
    massMap[findValue(fit,2000015,"mass","2000015")] = 2000015;
    massMap[getParticleData(24)->mass()/GeV        ] =      24;
    vector<int> move;
    for(map<double,long>::iterator mit=massMap.begin();mit!=massMap.end();++mit) {
      if     (mit->second==     37) move.push_back(0);
      else if(mit->second==1000011) move.push_back(1);
      else if(mit->second==1000013) move.push_back(2);
      else if(mit->second==1000015) move.push_back(3);
      else if(mit->second==2000011) move.push_back(4);
      else if(mit->second==2000013) move.push_back(5);
      else if(mit->second==2000015) move.push_back(6);
      else if(mit->second==     24) move.push_back(7);
    }
    CMatrix oldMat = ChargedHiggsMix()->getMatrix();
    CMatrix newMat(8,vector<Complex>(8,0.));
    for(unsigned int ix=0;ix<8;++ix) {
      for(unsigned int iy=0;iy<8;++iy)
	newMat[move[ix]][iy] = oldMat[ix][iy];
    }
    ChargedHiggsMix(new_ptr(MixingMatrix(newMat,ChargedHiggsMix()->getIds())));
    // adjust the pseudoscalars
    massMap.clear();
    massMap[findValue(fit,     36,"mass",     "36")] =      36;
    // extract the pseudoscalar masses and change the ids for the pseudoscalars
    // to those from the SLHA if needed
    if(fit->second.find(2000012)!=fit->second.end()) {
      massMap[findValue(fit,2000012,"mass","2000012")] = 1000017;
      idMap().insert(make_pair(2000012,1000017));
    }
    else {
      massMap[findValue(fit,1000017,"mass","1000017")] = 1000017;
    }
    if(fit->second.find(2000014)!=fit->second.end()) {
      massMap[findValue(fit,2000014,"mass","2000014")] = 1000018;
      idMap().insert(make_pair(2000014,1000018));
    }
    else {
      massMap[findValue(fit,1000018,"mass","1000018")] = 1000018;
    }
    if(fit->second.find(2000016)!=fit->second.end()) {
      massMap[findValue(fit,2000016,"mass","2000016")] = 1000019;
      idMap().insert(make_pair(2000016,1000019));
    }
    else {
      massMap[findValue(fit,1000019,"mass","1000019")] = 1000019;
    }
    move.clear(); move.push_back(4);
    for(map<double,long>::iterator mit=massMap.begin();mit!=massMap.end();++mit) {
      if     (mit->second==     36) move.push_back(0);
      else if(mit->second==1000017) move.push_back(1);
      else if(mit->second==1000018) move.push_back(2);
      else if(mit->second==1000019) move.push_back(3);
    }
    oldMat = CPoddHiggsMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[move[ix]][iy] = oldMat[ix][iy];
    }
    CPoddHiggsMix(new_ptr(MixingMatrix(newMat,CPoddHiggsMix()->getIds())));
    // adjust the neutral scalars
    massMap.clear();
    massMap[findValue(fit,     25,"mass",     "25")] =      25;
    massMap[findValue(fit,     35,"mass",     "35")] =      35;
    massMap[findValue(fit,1000012,"mass","1000012")] = 1000012;
    massMap[findValue(fit,1000014,"mass","1000014")] = 1000014;
    massMap[findValue(fit,1000016,"mass","1000016")] = 1000016;
    move.clear();
    for(map<double,long>::iterator mit=massMap.begin();mit!=massMap.end();++mit) {
      if     (mit->second==     25) move.push_back(0);
      else if(mit->second==     35) move.push_back(1);
      else if(mit->second==1000012) move.push_back(2);
      else if(mit->second==1000014) move.push_back(3);
      else if(mit->second==1000016) move.push_back(4);
    }
    oldMat = CPevenHiggsMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[move[ix]][iy] = oldMat[ix][iy];
    }
    CPevenHiggsMix(new_ptr(MixingMatrix(newMat,CPevenHiggsMix()->getIds())));
    // neutralino mixing
    move.resize(7);
    move[0] = 3; move[1] = 4; move[2] = 5; move[3] = 6;
    move[4] = 0; move[5] = 1; move[6] = 2;
    oldMat = neutralinoMix()->getMatrix();
    newMat = CMatrix(7,vector<Complex>(7,0.));
    for(unsigned int ix=0;ix<7;++ix) {
      for(unsigned int iy=0;iy<7;++iy)
	newMat[ix][move[iy]] = oldMat[ix][iy];
    }
    neutralinoMix(new_ptr(MixingMatrix(newMat,neutralinoMix()->getIds())));
    // chargino mixing
    move.resize(5);
    move[0] = 3; move[1] = 4;
    move[2] = 0; move[3] = 1; move[4] = 2;
    oldMat = charginoUMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[ix][move[iy]] = oldMat[ix][iy];
    }
    charginoUMix(new_ptr(MixingMatrix(newMat,charginoUMix()->getIds())));
    oldMat = charginoVMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[ix][move[iy]] = oldMat[ix][iy];
    }
    charginoVMix(new_ptr(MixingMatrix(newMat,charginoVMix()->getIds())));
  }
  // we don't want neutrinos first then neutralinos so swap them
  // neutralinos first then neutrinos
  if ( neutralinoMix()->size().first == 7 ) {
    vector<int> move(7);
    move[0] = 4; move[1] = 5; move[2] = 6;
    move[3] = 0; move[4] = 1; move[5] = 2; move[6] = 3;
    CMatrix oldMat = neutralinoMix()->getMatrix();
    CMatrix newMat(7,vector<Complex>(7,0.));
    for(unsigned int ix=0;ix<7;++ix) {
      for(unsigned int iy=0;iy<7;++iy)
	newMat[move[ix]][move[iy]] = oldMat[ix][iy];
    }
    neutralinoMix(new_ptr(MixingMatrix(newMat,neutralinoMix()->getIds())));
  }
  // charginos the same, i.e. charginos first then charged leptons
  if(charginoUMix()->size().first  != charginoVMix()->size().first || 
     charginoUMix()->size().second != charginoVMix()->size().second )
    throw Exception() << "Chargino U and V mixing matrices must have the same size.\n"
		      << "Check your SLHA file!" << Exception::runerror;
  if ( charginoUMix()->size().first == 5 ) {
    vector<int> move(5);
    move[0] = 2; move[1] = 3; move[2] = 4;
    move[3] = 0; move[4] = 1;
    CMatrix oldMat = charginoUMix()->getMatrix();
    CMatrix newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[move[ix]][move[iy]] = oldMat[ix][iy];
    }
    charginoUMix(new_ptr(MixingMatrix(newMat,charginoUMix()->getIds())));
    oldMat = charginoVMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[move[ix]][move[iy]] = oldMat[ix][iy];
    }
    charginoVMix(new_ptr(MixingMatrix(newMat,charginoVMix()->getIds())));
  }

  const MatrixSize & n = neutralinoMix()->size();
  const MatrixSize & u = charginoUMix()->size();
  const MatrixSize & h = CPevenHiggsMix()->size();
  const MatrixSize & a = CPoddHiggsMix()->size();
  const MatrixSize & l = ChargedHiggsMix()->size();

  bool nBig   = n.first == 7 && n.second == 7;
  bool nSmall = n.first == 4 && n.second == 4;
  if ( ! (nBig || nSmall) )
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "(RV)Nmix " << n.first << ',' << n.second << '\n'
      << Exception::runerror; 

  bool uBig   = u.first == 5 && u.second == 5;
  bool uSmall = u.first == 2 && u.second == 2;
  if ( ! (uBig || uSmall) )
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "(RV)Umix " << u.first << ',' << u.second << '\n'
      << Exception::runerror; 

  bool hBig   = h.first == 5 && h.second == 5;
  bool hSmall = h.first == 2 && h.second == 2;
  if ( ! (hBig || hSmall) )
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "(RV)Hmix " << h.first << ',' << h.second << '\n'
      << Exception::runerror; 

  bool aBig = (a.first == 4 || a.first == 5) && a.second == 5;
  bool aSmall = a.first == 1 && a.second == 2;
  if ( ! (aBig || aSmall) )
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "RVAmix " << a.first << ',' << a.second << '\n'
      << Exception::runerror; 

  bool lBig = (l.first == 7 || l.first == 8) && l.second == 8;
  bool lSmall = l.first == 1 && l.second == 2;
  if ( ! (lBig || lSmall) )
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "RVLmix " << l.first << ',' << l.second << '\n'
      << Exception::runerror; 


  bool allBig   = nBig && uBig && hBig && aBig && lBig;
  bool allSmall = nSmall && uSmall && hSmall && aSmall && lSmall;

  bool allSmallExceptN = nBig && uSmall && hSmall && aSmall && lSmall;
  
  if ( allSmallExceptN ) {
    cerr << "Warning: Truncating Nmix to 4,4 for consistency "
	 << "with other mixing matrices.\n";
    CMatrix oldMat = neutralinoMix()->getMatrix();
    CMatrix newMat(4,vector<Complex>(4,0.));
    for(unsigned int ix=0;ix<4;++ix)
      for(unsigned int iy=0;iy<4;++iy)
	newMat[ix][iy] = oldMat[ix][iy];
    assert( neutralinoMix()->getIds().size() >= 4 );
    vector<long>::const_iterator beg = neutralinoMix()->getIds().begin();
    vector<long>::const_iterator end = beg + 4;
    neutralinoMix(new_ptr(MixingMatrix(newMat,vector<long>(beg,end))));
    return;
  }
  else if ( ! (allBig || allSmall) ) {
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "(RV)Nmix " << n.first << ',' << n.second << '\n'
      << "(RV)Umix " << u.first << ',' << u.second << '\n'
      << "(RV)Hmix " << h.first << ',' << h.second << '\n'
      << "RVAmix   " << a.first << ',' << a.second << '\n'
      << "RVLmix   " << l.first << ',' << l.second << '\n'
      << Exception::runerror;
  }
  // reduce to MSSM + trilinear if requested
  if(triLinearOnly_) {
    // reduce size of n
    if(neutralinoMix()->size().first ==7) {
      CMatrix mix(4,vector<Complex>(4,0.));
      unsigned int irow=0;
      int imax[7]={-1,-1,-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<7;++ix) {
      	double maxComp(0.);
      	for(unsigned int iy=0;iy<7;++iy) {
       	  double value = abs((*neutralinoMix())(ix,iy));
       	  if(value>maxComp) {
       	    maxComp = value;
       	    imax[ix] = iy;
       	  }
      	}
	// neutralino
	if(imax[ix]<=3) {
	  for(unsigned int iy=0;iy<4;++iy) mix[irow][iy] = (*neutralinoMix())(ix,iy);
	  ++irow;
	  assert(irow<=4);
	}
	// neutrino
	else {
	  idMap()[neutralinoMix()->getIds()[ix]] = neutralinoMix()->getIds()[imax[ix]];
	}
      }
      vector<long> ids = neutralinoMix()->getIds();
      ids.resize(4);
      neutralinoMix(new_ptr(MixingMatrix(mix,ids)));
    }
    // reduce size of u
    if(charginoUMix()->size().first ==5) {
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int irow=0;
      int imax[5]={-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<5;++ix) {
	double maxComp(0.);
	for(unsigned int iy=0;iy<5;++iy) {
       	  double value = abs((*charginoUMix())(ix,iy));
       	  if(value>maxComp) {
       	    maxComp = value;
       	    imax[ix] = iy;
       	  }
      	}
	// chargino
       	if(imax[ix]<=1) {
     	  for(unsigned int iy=0;iy<2;++iy) mix[irow][iy] = (*charginoUMix())(ix,iy);
     	  ++irow;
     	  assert(irow<=2);
      	}
	// charged lepton
      	else {
      	  idMap()[abs(charginoUMix()->getIds()[ix])] = abs(charginoUMix()->getIds()[imax[ix]]);
     	}
      }
      vector<long> ids = charginoUMix()->getIds();
      ids.resize(2);
      charginoUMix(new_ptr(MixingMatrix(mix,ids)));
    }
    // reduce size of v
    if(charginoVMix()->size().first ==5) {
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int irow=0;
      int imax[5]={-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<5;++ix) {
	double maxComp(0.);
	for(unsigned int iy=0;iy<5;++iy) {
       	  double value = abs((*charginoVMix())(ix,iy));
       	  if(value>maxComp) {
       	    maxComp = value;
       	    imax[ix] = iy;
       	  }
      	}
	// chargino
       	if(imax[ix]<=1) {
     	  for(unsigned int iy=0;iy<2;++iy) mix[irow][iy] = (*charginoVMix())(ix,iy);
     	  ++irow;
     	  assert(irow<=2);
      	}
      }
      vector<long> ids = charginoVMix()->getIds();
      ids.resize(2);
      charginoVMix(new_ptr(MixingMatrix(mix,ids)));
    }
    // reduce size of pseudo scalar mixing
    if(CPoddHiggsMix()) {
      MixingVector hmix;
      double beta = atan(tanBeta());
      hmix.push_back(MixingElement(1,1,sin(beta)));
      hmix.push_back(MixingElement(1,2,cos(beta)));
      vector<long> ids(1,36);
      MixingMatrixPtr newMix = new_ptr(MixingMatrix(1,2));
      (*newMix).setIds(ids);
      for(unsigned int ix=0; ix < hmix.size(); ++ix)
	(*newMix)(hmix[ix].row-1,hmix[ix].col-1) = hmix[ix].value;
      CPoddHiggsMix(newMix);
    }
    // reduce size of true scalar mixing
    if(CPevenHiggsMix()->size().first==5) {
      double alpha(0.);
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int irow=0;
      int imax[5]={-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<5;++ix) {
      	double maxComp(0.);
      	for(unsigned int iy=0;iy<5;++iy) {
       	  double value = abs((*CPevenHiggsMix())(ix,iy));
       	  if(value>maxComp) {
       	    maxComp = value;
       	    imax[ix] = iy;
       	  }
      	}
       	// neutral Higgs
       	if(imax[ix]<=1) {
      	  for(unsigned int iy=0;iy<2;++iy) mix[irow][iy] = (*CPevenHiggsMix())(ix,iy);
	  if(irow==0)
	    alpha = atan2(-(*CPevenHiggsMix())(ix,0).real(),(*CPevenHiggsMix())(ix,1).real());
      	  ++irow;
      	  assert(irow<=2);
       	}
      }
      vector<long> ids = CPevenHiggsMix()->getIds();
      ids.resize(2);
      CPevenHiggsMix(new_ptr(MixingMatrix(mix,ids)));
      higgsMixingAngle(alpha);
    }
    else {
      bool readAlpha = false;
      map<string,ParamMap>::const_iterator pit=parameters().find("alpha");
      if(pit!=parameters().end()) {
	ParamMap::const_iterator it = pit->second.find(1);
	if(it!=pit->second.end()) {
	  readAlpha = true;
	  higgsMixingAngle(it->second);
	}
      }
      if(!readAlpha) 
	throw Exception() << "In the RPV model BLOCK ALPHA which must be"
			  << " present in the SLHA file is missing"
			  << Exception::runerror;
    }
    // reduce size of charged scalar mixing
    if(ChargedHiggsMix()->size().first>=7) {
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int istau(0);
      int imax[8]={-1,-1,-1,-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<ChargedHiggsMix()->size().first;++ix) {
	double maxComp(0.);
	for(unsigned int iy=0;iy<8;++iy) {
	  double value = abs((*ChargedHiggsMix())(ix,iy));
	  if(value>maxComp) {
	    maxComp = value;
	    imax[ix] = iy;
	  }
	}
	if(imax[ix]<=1) imax[ix]=0;
	else --imax[ix];
	idMap()[abs(ChargedHiggsMix()->getIds()[ix])] = abs(ChargedHiggsMix()->getIds()[imax[ix]]);
      	if(abs(ChargedHiggsMix()->getIds()[imax[ix]])%10==5) {
      	  mix[istau][0] = (*ChargedHiggsMix())(ix,4);
      	  mix[istau][1] = (*ChargedHiggsMix())(ix,7);
	  if(istau==0) idMap()[abs(ChargedHiggsMix()->getIds()[ix])] = 1000015;
	  else         idMap()[abs(ChargedHiggsMix()->getIds()[ix])] = 2000015;
      	  istau+=1;
      	  assert(istau<=2);
      	}
      }
      // set up the stau mixing matrix
      vector<long> ids(2);
      ids[0] = 1000015;
      ids[1] = 2000015;
      stauMix(new_ptr(MixingMatrix(mix,ids)));
      // delete 7x7
      MixingVector hmix;
      double beta = atan(tanBeta());
      hmix.push_back(MixingElement(1,1,sin(beta)));
      hmix.push_back(MixingElement(1,2,cos(beta)));
      ids.clear();
      ids.resize(1,37);
      MixingMatrixPtr newMix = new_ptr(MixingMatrix(1,2));
      (*newMix).setIds(ids);
      for(unsigned int ix=0; ix < hmix.size(); ++ix)
      	(*newMix)(hmix[ix].row-1,hmix[ix].col-1) = hmix[ix].value;
      ChargedHiggsMix(newMix);
    }
    // reduce size of up squark mixing
    if( upSquarkMix_ ) {
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int istop(0);
      int imax[6]={-1,-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<6;++ix) {
	double maxComp(0.);
	for(unsigned int iy=0;iy<6;++iy) {
	  double value = abs((*upSquarkMix_)(ix,iy));
	  if(value>maxComp) {
	    maxComp = value;
	    imax[ix] = iy;
	  }
	}
	idMap()[upSquarkMix_->getIds()[ix]] = upSquarkMix_->getIds()[imax[ix]];
	if(upSquarkMix_->getIds()[imax[ix]]%10==6) {
	  mix[istop][0] = (*upSquarkMix_)(ix,2);
	  mix[istop][1] = (*upSquarkMix_)(ix,5);
	  if(istop==0) idMap()[upSquarkMix_->getIds()[ix]] = 1000006;
	  else         idMap()[upSquarkMix_->getIds()[ix]] = 2000006;
	  istop+=1;
	  assert(istop<=2);
	}
      }
      // set up the stop mixing matrix
      vector<long> ids(2);
      ids[0] = 1000006;
      ids[1] = 2000006;
      stopMix(new_ptr(MixingMatrix(mix,ids)));
      // delete 6x6
      upSquarkMix_ = MixingMatrixPtr();
    }
    // reduce size of down squark mixing
    if( downSquarkMix_ ) {
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int isbot(0);
      int imax[6]={-1,-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<6;++ix) {
      	double maxComp(0.);
      	for(unsigned int iy=0;iy<6;++iy) {
      	  double value = abs((*downSquarkMix_)(ix,iy));
       	  if(value>maxComp) {
       	    maxComp = value;
       	    imax[ix] = iy;
       	  }
       	}
      	idMap()[downSquarkMix_->getIds()[ix]] = downSquarkMix_->getIds()[imax[ix]];
      	if(downSquarkMix_->getIds()[imax[ix]]%10==5) {
      	  mix[isbot][0] = (*downSquarkMix_)(ix,2);
      	  mix[isbot][1] = (*downSquarkMix_)(ix,5);
	  if(isbot==0) idMap()[downSquarkMix_->getIds()[ix]] = 1000005;
	  else         idMap()[downSquarkMix_->getIds()[ix]] = 2000005;
      	  isbot+=1;
      	  assert(isbot<=2);
      	}
      }
      // set up the sbottom mixing matrix
      vector<long> ids(2);
      ids[0] = 1000005;
      ids[1] = 2000005;
      sbottomMix(new_ptr(MixingMatrix(mix,ids)));
      // delete 6x6
      downSquarkMix_ = MixingMatrixPtr();
    }
    // get the masses as we will have to mess with them
    map<string,ParamMap>::const_iterator fit=parameters().find("mass");
    ParamMap theMasses = fit->second;
    for(long id=1000017;id<=1000019;++id) {
      if(theMasses.find(id)==theMasses.end()) continue;
      double mass = abs(theMasses.find(id)->second);
      // find scalar partner
      double mdiff=1e30;
      long new_id = 0;
      for(ParamMap::const_iterator it=theMasses.begin();it!=theMasses.end();++it) {
	double diff = abs(abs(it->second)-mass);
	if(diff<mdiff) {
	  mdiff = diff;
	  new_id = it->first;
	}
      }
      if(idMap().find(new_id)!=idMap().end()) {
	idMap()[id] = idMap()[new_id];
      }
      else {
	idMap()[id] = new_id;
      }
    }
  }
}

void RPV::doinit() {
  MSSM::doinit();
  addVertex(LLEVertex_);
  addVertex(LQDVertex_);
  addVertex(UDDVertex_);
}
#line 1 "./RPVLLEVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVLLEVertex class.
//

#include "RPVLLEVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

RPVLLEVertex::RPVLLEVertex() {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr RPVLLEVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVLLEVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVLLEVertex::persistentOutput(PersistentOStream & os) const {
  os << lambda_ << stau_;
}

void RPVLLEVertex::persistentInput(PersistentIStream & is, int) {
  is >> lambda_ >> stau_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<RPVLLEVertex,FFSVertex>
describeHerwigRPVLLEVertex("Herwig::RPVLLEVertex", "HwSusy.so HwRPV.so");

void RPVLLEVertex::Init() {

  static ClassDocumentation<RPVLLEVertex> documentation
    ("The RPVLLEVertex class implements the trilinear LLE"
     " coupling in R-parity violating models.");

}

void RPVLLEVertex::doinit() {
  RPVPtr rpv = dynamic_ptr_cast<RPVPtr>(generator()->standardModel());
  if(!rpv)
    throw InitException() << "Must have the RPV model in"
			  << " RPVLLEVertex::doinit()"
			  << Exception::abortnow;
  // get the coupling
  lambda_ = rpv->lambdaLLE();
  stau_   = rpv->stauMix();
  // set the particles in the vertex if coupling non-zero
  for( int i=0;i<3;++i) {
    for( int j=0; j<3; ++j) {
      for( int k=0; k<3; ++k) {
	// continue if zero
	if(lambda_[i][j][k]==0.) continue;
	// particles in the vertex
	// sneutrino
	addToList( -2*i-11,  2*k+11, -1000012-2*j );
	addToList(  2*i+11, -2*k-11,  1000012+2*j );
	// left slepton
	addToList( -2*i-12,  2*k+11, -1000011-2*j );
	addToList(  2*i+12, -2*k-11, +1000011+2*j );
	if(j==2) {
	  addToList( -2*i-12,  2*k+11, -2000011-2*j );
	  addToList(  2*i+12, -2*k-11, +2000011+2*j );
	}
	// right slepton
	addToList(  2*i+12,  2*j+11, -2000011-2*k );
	addToList( -2*i-12, -2*j-11, +2000011+2*k );
	if(k==2) {
	  addToList(  2*i+12,  2*j+11, -1000011-2*k );
	  addToList( -2*i-12, -2*j-11, +1000011+2*k );
	}
      }
    }
  }
  FFSVertex::doinit();
}

void RPVLLEVertex::setCoupling(Energy2, tcPDPtr part1,
			       tcPDPtr part2, tcPDPtr part3) {
  int islep = part3->id();
  int i(-1),j(-1),k(-1);
  Complex mix=1.;
  // sneutrino case
  if( abs(islep) == ParticleID::SUSY_nu_eL  || 
      abs(islep) == ParticleID::SUSY_nu_muL || 
      abs(islep) == ParticleID::SUSY_nu_tauL) {
    j = (abs(islep) - 1000012)/2;
    i = (abs(part1->id())-11)/2;
    k = (abs(part2->id())-11)/2;
    if(part1->id()*islep<0) swap(i,k);
    if(islep<0) {
      left (1.);
      right(0.);
    }
    else {
      left (0.);
      right(1.);
    }
  }
  // charged slepton case
  else if( abs(islep) == ParticleID::SUSY_e_Lminus  || 
	   abs(islep) == ParticleID::SUSY_mu_Lminus || 
	   abs(islep) == ParticleID::SUSY_e_Rminus  || 
	   abs(islep) == ParticleID::SUSY_mu_Rminus || 
	   abs(islep) == ParticleID::SUSY_tau_1minus|| 
	   abs(islep) == ParticleID::SUSY_tau_2minus) {
    // right charged slepton
    if(part1->id()*part2->id()>0) {
      if(abs(part1->id())%2==0) {
	i = (abs(part1->id())-12)/2;
	j = (abs(part2->id())-11)/2;
      }
      else {
	i = (abs(part2->id())-12)/2;
	j = (abs(part1->id())-11)/2;
      }
      if(abs(islep)>2000000) {
	k = (abs(islep)-2000011)/2;
	if(k==2) mix = (*stau_)(1,1);
      }
      else {
	assert(abs(islep)==ParticleID::SUSY_tau_1minus);
	k = 2;
	mix = (*stau_)(0,1);
      }
      if(islep>0) {
	left (-1.);
	right(0.);
      }
      else {
	left (0.);
	right(-1.);
      }
    }
    // left charged
    else {
      if(abs(part1->id())%2==0) {
	i = (abs(part1->id())-12)/2;
	k = (abs(part2->id())-11)/2;
      }
      else {
	i = (abs(part2->id())-12)/2;
	k = (abs(part1->id())-11)/2;
      }
      if(abs(islep)<2000000) {
	j = (abs(islep)-1000011)/2;
	if(j==2) mix = (*stau_)(0,0);
      }
      else {
	assert(abs(islep)==ParticleID::SUSY_tau_2minus);
	j = 2;
	mix = (*stau_)(1,0);
      }
      if(islep<0) {
	left (-1.);
	right(0.);
      }
      else {
	left (0.);
	right(-1.);
      }
    }
  }
  else
    assert(false);
  assert( i>=0 && i<=2 &&  j>=0 && j<=2 &&  k>=0 && k<=2 );
  norm(mix*lambda_[i][j][k]);
}
  
#line 1 "./RPVLQDVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVLQDVertex class.
//

#include "RPVLQDVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

RPVLQDVertex::RPVLQDVertex() {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

IBPtr RPVLQDVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVLQDVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVLQDVertex::persistentOutput(PersistentOStream & os) const {
  os << lambda_ << stop_ << sbot_ << stau_;
}

void RPVLQDVertex::persistentInput(PersistentIStream & is, int) {
  is >> lambda_ >> stop_ >> sbot_ >> stau_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<RPVLQDVertex,FFSVertex>
describeHerwigRPVLQDVertex("Herwig::RPVLQDVertex", "HwSusy.so HwRPV.so");

void RPVLQDVertex::Init() {

  static ClassDocumentation<RPVLQDVertex> documentation
    ("The RPVLQDVertex class implements the trilinear LQD"
     " coupling in R-parity violating models.");

}

void RPVLQDVertex::doinit() {
  RPVPtr rpv = dynamic_ptr_cast<RPVPtr>(generator()->standardModel());
  if(!rpv)
    throw InitException() << "Must have the RPV model in"
			  << " RPVLQDVertex::doinit()"
			  << Exception::abortnow;
  // get the coupling
  lambda_ = rpv->lambdaLQD();
  stop_   = rpv->stopMix();
  sbot_   = rpv->sbottomMix();
  stau_   = rpv->stauMix();
  // set the particles in the vertex if coupling non-zero
  for( int i=0;i<3;++i) {
    for( int j=0; j<3; ++j) {
      for( int k=0; k<3; ++k) {
	// continue if zero
	if(lambda_[i][j][k]==0.) continue;
	// particles in the vertex
	// sneutrino
	addToList( -2*j-1 ,  2*k+1 , -1000012-2*i );
	addToList(  2*j+1 , -2*k-1 ,  1000012+2*i );
	// left slepton
	addToList( -2*j-2 ,  2*k+1 , -1000011-2*i );
	addToList(  2*j+2 , -2*k-1 , +1000011+2*i );
	if(i==2) {
	  addToList( -2*j-2 ,  2*k+1 , -2000011-2*i );
	  addToList(  2*j+2 , -2*k-1 , +2000011+2*i );
	}
	// up squark
	addToList( -2*i-11,  2*k+1 , -1000002-2*j );
	addToList(  2*i+11, -2*k-1 ,  1000002+2*j );
	if(j==2) {
	  addToList( -2*i-11,  2*k+1 , -2000002-2*j );
	  addToList(  2*i+11, -2*k-1 ,  2000002+2*j );
	}
	// left down squark
	addToList( -2*i-12,  2*k+1 , -1000001-2*j );
	addToList(  2*i+12, -2*k-1 ,  1000001+2*j );
	if(j==2) {
	  addToList( -2*i-12,  2*k+1 , -2000001-2*j );
	  addToList(  2*i+12, -2*k-1 ,  2000001+2*j );
	}
	// right down squark
	addToList(  2*i+12,  2*j+1 , -2000001-2*k );
	addToList( -2*i-12, -2*j-1 , +2000001+2*k );
	if(k==2) {
	  addToList(  2*i+12,  2*j+1 , -1000001-2*k );
	  addToList( -2*i-12, -2*j-1 , +1000001+2*k );
	}
	addToList(  2*i+11,  2*j+2 , -2000001-2*k );
	addToList( -2*i-11, -2*j-2 , +2000001+2*k );
	if(k==2) {
	  addToList(  2*i+11,  2*j+2 , -1000001-2*k );
	  addToList( -2*i-11, -2*j-2 , +1000001+2*k );
	}
      }
    }
  }
  FFSVertex::doinit();
}

void RPVLQDVertex::setCoupling(Energy2, tcPDPtr part1,
			       tcPDPtr part2, tcPDPtr part3) {
  int islep = part3->id();
  int i(-1),j(-1),k(-1);
  Complex mix(1.);
  // sneutrino case
  if( abs(islep) == ParticleID::SUSY_nu_eL  || 
      abs(islep) == ParticleID::SUSY_nu_muL || 
      abs(islep) == ParticleID::SUSY_nu_tauL) {
    i = (abs(islep)-1000012)/2;
    j = (abs(part1->id())-1)/2;
    k = (abs(part2->id())-1)/2;
    if( islep*part1->id() <0 ) swap(j,k);
    if(islep>0) {
      left ( 0.);
      right(-1.);
    }
    else {
      left (-1.);
      right( 0.);
    }
  }
  else if(abs(islep) == ParticleID::SUSY_e_Lminus|| 
	  abs(islep) == ParticleID::SUSY_mu_Lminus|| 
	  abs(islep) == ParticleID::SUSY_tau_1minus|| 
	  abs(islep) == ParticleID::SUSY_tau_2minus) { 
    if(abs(islep)<2000000) {
      i = (abs(islep)-1000011)/2;
      if(i==2) mix = (*stau_)(0,0);
    }
    else {
      i = 2;
      mix = (*stau_)(1,0);
      assert(abs(islep)==ParticleID::SUSY_tau_2minus);
    }
    if(abs(part1->id())%2==1) {
      j = (abs(part2->id())-2)/2;
      k = (abs(part1->id())-1)/2;
    }
    else {
      j = (abs(part1->id())-2)/2;
      k = (abs(part2->id())-1)/2;
    }
    if(islep<0) {
      left ( 1.);
      right( 0.);
    }
    else {
      left ( 0.);
      right( 1.);
    }
  }
  // left up squark
  else if(abs(islep) == ParticleID::SUSY_u_L  || 
	  abs(islep) == ParticleID::SUSY_c_L || 
	  abs(islep) == ParticleID::SUSY_t_1 || 
	  abs(islep) == ParticleID::SUSY_t_2) {
    if(abs(islep)<2000000) {
      j = (abs(islep)-1000002)/2;
      if(j==2) mix = (*stop_)(0,0); 
    }
    else {
      j = 2;
      mix = (*stop_)(1,0);
      assert(abs(islep)==ParticleID::SUSY_t_2);
    }
    if(part1->coloured()) {
      i = (abs(part2->id())-11)/2;
      k = (abs(part1->id())- 1)/2;
    }
    else {
      i = (abs(part1->id())-11)/2;
      k = (abs(part2->id())- 1)/2;
    }
    if(islep<0) {
      left ( 1.);
      right( 0.);
    }
    else {
      left ( 0.);
      right( 1.);
    }
  }
  else if(abs(islep) == ParticleID::SUSY_d_L || 
	  abs(islep) == ParticleID::SUSY_s_L || 
	  abs(islep) == ParticleID::SUSY_b_1 || 
	  abs(islep) == ParticleID::SUSY_b_2 ||
	  abs(islep) == ParticleID::SUSY_d_R || 
	  abs(islep) == ParticleID::SUSY_s_R) {
    // right down squark
    if(part1->id()*part2->id()>0) {
      if(part1->coloured()) {
	if(abs(part1->id())%2==1) {
	  i = (abs(part2->id())-12)/2;
	  j = (abs(part1->id())-1 )/2;
	  mix *= -1.;
	  assert(i>=0);
	}
	else {
	  i = (abs(part2->id())-11)/2;
	  j = (abs(part1->id())-2 )/2;
	  assert(i>=0);
	}
      }
      else {
	if(abs(part2->id())%2==1) {
	  i = (abs(part1->id())-12)/2;
	  j = (abs(part2->id())-1 )/2;
	  mix *= -1.;
	  assert(i>=0);
	}
	else {
	  i = (abs(part1->id())-11)/2;
	  j = (abs(part2->id())-2 )/2;
	  assert(i>=0);
	}
      }
      if(abs(islep)>2000000) {
	k = (abs(islep)-2000001)/2;
	if(k==2) mix *= (*sbot_)(1,1);
      }
      else {
	k = 2;
	mix *= (*sbot_)(0,1);
	assert(abs(islep)==ParticleID::SUSY_b_1);
      }
      if(islep<0) {
	left ( 0.);
	right( 1.);
      }
      else {
	left ( 1.);
	right( 0.);
      }
    }
    // left   down  squark
    else {
      if(abs(islep)<2000000) {
	j = (abs(islep)-1000001)/2;
	if(j==2) mix = (*sbot_)(0,0); 
      }
      else {
	j = 2;
	mix = (*sbot_)(1,0);
	assert(abs(islep)==ParticleID::SUSY_b_2);
      }
      if(part1->coloured()) {
	i = (abs(part2->id())-12)/2;
	k = (abs(part1->id())- 1)/2;
      }
      else {
	i = (abs(part1->id())-12)/2;
	k = (abs(part2->id())- 1)/2;
      }
      if(islep<0) {
	left (-1.);
	right( 0.);
      }
      else {
	left ( 0.);
	right(-1.);
      }
      assert(i>=0);
    }
  }
  else
    assert(false);
  assert( i>=0 && i<=2 &&  j>=0 && j<=2 &&  k>=0 && k<=2 );
  norm(mix*lambda_[i][j][k]);
}
#line 1 "./RPVUDDVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UDDVertex class.
//

#include "RPVUDDVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

RPVUDDVertex::RPVUDDVertex() {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::EPS);
}

IBPtr RPVUDDVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVUDDVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVUDDVertex::persistentOutput(PersistentOStream & os) const {
  os << lambda_ << stop_ << sbot_;
}

void RPVUDDVertex::persistentInput(PersistentIStream & is, int) {
  is >> lambda_ >> stop_ >> sbot_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<RPVUDDVertex,FFSVertex>
describeHerwigRPVUDDVertex("Herwig::RPVUDDVertex", "HwSusy.so HwRPV.so");

void RPVUDDVertex::Init() {

  static ClassDocumentation<RPVUDDVertex> documentation
    ("The RPVUDDVertex class implements the trilinear UDD"
     " coupling in R-parity violating models.");

}

void RPVUDDVertex::doinit() {
  RPVPtr rpv = dynamic_ptr_cast<RPVPtr>(generator()->standardModel());
  if(!rpv)
    throw InitException() << "Must have the RPV model in"
			  << " RPVUDDVertex::doinit()"
			  << Exception::abortnow;
  // get the coupling
  lambda_ = rpv->lambdaUDD();
  stop_   = rpv->stopMix();
  sbot_   = rpv->sbottomMix();
  // set the particles in the vertex if coupling non-zero
  for( int i=0;i<3;++i) {
    for( int j=0; j<3; ++j) {
      for( int k=0; k<3; ++k) {
 	// continue if zero
 	if(lambda_[i][j][k]==0.) continue;
 	// particles in the vertex
 	// right up squark
	if(j<k) {
	  addToList( -2*j-1 , -2*k-1 , -2000002-2*i );
	  addToList(  2*j+1 ,  2*k+1 , +2000002+2*i );
	}
	if(i==2) {
	  addToList( -2*j-1 , -2*k-1 , -1000002-2*i );
	  addToList(  2*j+1 ,  2*k+1 , +1000002+2*i );
	}
 	// right down squark
	addToList( -2*i-2 , -2*j-1 , -2000001-2*k );
	addToList(  2*i+2 ,  2*j+1 , +2000001+2*k );
	if(k==2) {
	  addToList( -2*i-2 , -2*j-1 , -1000001-2*k );
	  addToList(  2*i+2 ,  2*j+1 , +1000001+2*k );
 	}
      }
    }
  }
  FFSVertex::doinit();
}

void RPVUDDVertex::setCoupling(Energy2, tcPDPtr part1,
			       tcPDPtr part2, tcPDPtr part3) {
  int islep = part3->id();
  int i(-1),j(-1),k(-1);
  Complex mix(1.);
  // left up squark
  if(abs(islep) == ParticleID::SUSY_u_R  || 
     abs(islep) == ParticleID::SUSY_c_R || 
     abs(islep) == ParticleID::SUSY_t_1 || 
     abs(islep) == ParticleID::SUSY_t_2) {
    if(abs(islep)>2000000) {
      i = (abs(islep)-2000002)/2;
      if(i==2) mix = (*stop_)(1,1); 
    }
    else {
      i = 2;
      mix = (*stop_)(0,1);
      assert(abs(islep)==ParticleID::SUSY_t_1);
    }
    j = (abs(part1->id())- 1)/2;
    k = (abs(part2->id())- 1)/2;
  }
  else if(abs(islep) == ParticleID::SUSY_d_R || 
	  abs(islep) == ParticleID::SUSY_s_R || 
	  abs(islep) == ParticleID::SUSY_b_1 || 
	  abs(islep) == ParticleID::SUSY_b_2) {
    if(abs(islep)>2000000) {
      k = (abs(islep)-2000001)/2;
      if(k==2) mix = (*sbot_)(1,1);
    }
    else {
      k = 2;
      mix = (*sbot_)(0,1);
      assert(abs(islep)==ParticleID::SUSY_b_1);
    }
    if(abs(part1->id())%2==0) {
      i = (abs(part1->id())- 2)/2;
      j = (abs(part2->id())- 1)/2;
    }
    else {
      i = (abs(part2->id())- 2)/2;
      j = (abs(part1->id())- 1)/2;
      // mix *= -1.;
    }
  }
  else 
    assert(false);
  assert( i>=0 && i<=2 &&  j>=0 && j<=2 &&  k>=0 && k<=2 );
  if(islep>0) {
    left (1.);
    right(0.);
  }
  else {
    left (0.);
    right(1.);
  }
  norm(mix*lambda_[i][j][k]);
}
#line 1 "./RPVFFZVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVFFZVertex class.
//

#include "RPVFFZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "RPVhelper.h"

using namespace Herwig;

RPVFFZVertex::RPVFFZVertex()  : _sw(0.), _cw(0.), _id1last(0), 
				_id2last(0), _q2last(), _couplast(0.),
				_leftlast(0.), _rightlast(0.),
				_gl(17,0.0), _gr(17,0.0), _gblast(0),
				_interactions(0) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

IBPtr RPVFFZVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVFFZVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVFFZVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _cw << _theN << _theU << _theV 
     << _gl << _gr << _interactions;
}

void RPVFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _theN >> _theU >> _theV 
     >> _gl >> _gr >> _interactions;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RPVFFZVertex,Helicity::FFVVertex>
describeHerwigRPVFFZVertex("Herwig::RPVFFZVertex", "HwSusy.so HwRPV.so");

void RPVFFZVertex::Init() {

  static ClassDocumentation<RPVFFZVertex> documentation
    ("The RPVFFZVertex class implements trhe coupling of the Z to all"
     " fermion-antifermion pairs in models with bilinear RPV.");

  static Switch<RPVFFZVertex,unsigned int> interfaceInteractions
    ("Interactions",
     "Which interactions to include",
     &RPVFFZVertex::_interactions, 0, false, false);
  static SwitchOption interfaceInteractionsAll
    (interfaceInteractions,
     "All",
     "Include all the interactions",
     0);
  static SwitchOption interfaceInteractionsSM
    (interfaceInteractions,
     "SM",
     "Only include what would have been the interactions with the SM"
     " fermions in the absence of mixing",
     1);
  static SwitchOption interfaceInteractionsNeutralino
    (interfaceInteractions,
     "Neutralino",
     "Only include what would have been the interactions with the "
     "neutralinos in the absence of mixing",
     2);
  static SwitchOption interfaceInteractionsChargino
    (interfaceInteractions,
     "Chargino",
     "Only include what would have been the interactions with the "
     "charginos in the absence of mixing",
     3);
}

void RPVFFZVertex::doinit() {
  // extract the mixing matrices
  tSusyBasePtr model = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!model) throw InitException() << "RPVFFZVertex::doinit() - "
				   << "The model pointer is null."
				   << Exception::abortnow;
  _theN = model->neutralinoMix();
  _theU = model->charginoUMix();
  _theV = model->charginoVMix();
  if( !_theN || !_theU || !_theV )
    throw InitException() << "RPVFFZVertex::doinit - "
			  << "A mixing matrix pointer is null.  U: " 
			  << _theU << "  V: " << _theV << "  N: " << _theN
			  << Exception::abortnow;
  // Standard Model fermions
  if(_interactions==0||_interactions==1) {
    // PDG codes for the particles
    // the quarks
    for(int ix=1;ix<7;++ix) {
      addToList(-ix, ix, 23);
    }
    // the leptons
    for(int ix=11;ix<17;ix+=2) {
      addToList(-ix, ix, 23);
    }
    for(int ix=12;ix<17;ix+=2) {
      if(_theN->size().first==7) {
	long inu = (ix-12)/2+17;
	addToList( inu, inu, 23);
      }
      else
	addToList(-ix, ix, 23);
    }
  }
  // neutralinos
  if(_interactions==0||_interactions==2) {
    vector<long> neu(4);
    neu[0] =  1000022; neu[1] = 1000023;
    neu[2] =  1000025; neu[3] = 1000035;
    if(_theN->size().first==7) {
      if(model->majoranaNeutrinos()) {
	neu.push_back(17);
	neu.push_back(18);
	neu.push_back(19);
      }
      else {
	neu.push_back(12);
	neu.push_back(14);
	neu.push_back(16);
      }
    }
    for(unsigned int i = 0; i < neu.size(); ++i) {
      for(unsigned int j = 0; j < neu.size(); ++j) {
	if(!(i>3&&i==j)) addToList(neu[i], neu[j], 23);
      }
    }
  }
  // charginos
  if(_interactions==0||_interactions==3) {
    addToList(-1000024, 1000024, 22);
    addToList(-1000037, 1000037, 22);
    vector<long> cha(2);
    cha[0] = 1000024; cha[1] = 1000037;
    if(_theV->size().first==5) {
      cha.push_back(-11);
      cha.push_back(-13);
      cha.push_back(-15);
    }
    for(unsigned int i = 0; i < cha.size(); ++i) {
      for(unsigned int j = 0; j < cha.size(); ++j) {
	if(!(i>1&&i==j)) addToList(-cha[i], cha[j], 23);
      }
    }
  }
  Helicity::FFVVertex::doinit();
  // weak mixing
  double sw2 = sin2ThetaW();
  _cw  = sqrt(1. - sw2);
  _sw  = sqrt(   sw2  );
  // Standard Model couplings
  for(int ix=1;ix<4;++ix) {
    _gl[2*ix-1]  = -0.25*(model->vd()  + model->ad() );
    _gl[2*ix ]   = -0.25*(model->vu()  + model->au() );
    _gl[2*ix+9 ] = -0.25*(model->ve()  + model->ae() );
    _gl[2*ix+10] = -0.25*(model->vnu() + model->anu());
    _gr[2*ix-1]  = -0.25*(model->vd()  - model->ad() );
    _gr[2*ix ]   = -0.25*(model->vu()  - model->au() );
    _gr[2*ix+9 ] = -0.25*(model->ve()  - model->ae() );
    _gr[2*ix+10] = -0.25*(model->vnu() - model->anu());
  }
}

void RPVFFZVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			       tcPDPtr part2,tcPDPtr part3) {
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = electroMagneticCoupling(q2);
    _q2last=q2;
  }
  long iferm1(part1->id()), iferm2(part2->id()), boson(part3->id());
  long iferm = abs(iferm1);
  // chargino coupling to the photon
  if(part3->id()==ParticleID::gamma) {
    assert(iferm == abs(iferm2));
    _gblast = boson;
    _id1last = iferm1;
    _id2last = iferm2;
    _leftlast  = -1.;
    _rightlast = -1.;
    if(iferm1>0) {
      Complex temp = _leftlast;
      _leftlast  = -_rightlast;
      _rightlast = -temp;
    }
  }
  // coupling to the Z
  else {
    assert(part3->id()==ParticleID::Z0);
    // quarks
    if(iferm<=6) {
      _leftlast  = _gl[iferm]/(_sw*_cw);
      _rightlast = _gr[iferm]/(_sw*_cw);
    }
    // charged leptons and charginos
    else if(part1->iCharge()!=0) {
      if(boson != _gblast || iferm1 != _id1last || iferm2 != _id2last) {
	_gblast = boson;
	_id1last = iferm1;
	_id2last = iferm2;
	unsigned int ic1(0);
	if(_theV->size().first==2&&iferm<=16) {
	  _leftlast  = -_gr[iferm];
	  _rightlast = -_gl[iferm];
	}
	else {
	  ic1 = RPV_helper::charginoIndex(iferm1);
	  unsigned int ic2 = RPV_helper::charginoIndex(iferm2);
	  _leftlast = -(*_theV)(ic1, 0)*conj((*_theV)(ic2, 0)) - 
	    0.5*(*_theV)(ic1, 1)*conj((*_theV)(ic2, 1));
	  _rightlast = -conj((*_theU)(ic1, 0))*(*_theU)(ic2, 0) - 
	    0.5*conj((*_theU)(ic1, 1))*(*_theU)(ic2, 1);
	  if(abs(iferm1) == abs(iferm2)) {
	    _leftlast  += sqr(_sw);
	    _rightlast += sqr(_sw);
	  }
	  if(_theV->size().first==5) {
	    for(unsigned int ix=0;ix<3;++ix) {
	      _rightlast += -0.5*(*_theU)(ic1, 2+ix)*conj((*_theU)(ic2, 2+ix));
	    }
	  }
	}
	if((ic1<2&&iferm1>0)||(ic1>=2&&iferm1<0)) {
	  Complex temp = _leftlast;
	  _leftlast  = -_rightlast;
	  _rightlast = -temp;
	}
	Complex temp = _leftlast;
	_leftlast  = -_rightlast;
	_rightlast = -temp;
	_leftlast  /= _sw*_cw;
	_rightlast /= _sw*_cw;
      }
    }
    // neutrinos and neutralinos
    else {
      // case where only 4x4 matrix and neutrino
      if(_theN->size().first==4&&iferm<=16) {
	assert(iferm==12||iferm==14||iferm==16);
	_leftlast  = _gl[iferm]/(_sw*_cw);
	_rightlast = _gr[iferm]/(_sw*_cw);
      }
      // neutralino
      else {
	long ic1 = part2->id();
	long ic2 = part1->id();
	assert(ic1 == ParticleID::SUSY_chi_10 || ic1 == ParticleID::SUSY_chi_20 ||
	       ic1 == ParticleID::SUSY_chi_30 || ic1 == ParticleID::SUSY_chi_40 ||
	       abs(ic1) == 12 || abs(ic1) == 14 || abs(ic1) == 16 ||
	       abs(ic1) == 17 || abs(ic1) == 18 || abs(ic1) == 19 );
	assert(ic2 == ParticleID::SUSY_chi_10 || ic2 == ParticleID::SUSY_chi_20 ||
	       ic2 == ParticleID::SUSY_chi_30 || ic2 == ParticleID::SUSY_chi_40 ||
	       abs(ic2) == 12 || abs(ic2) == 14 || abs(ic2) == 16 ||
	       abs(ic2) == 17 || abs(ic2) == 18 || abs(ic2) == 19 );
	if(ic1 != _id1last || ic2 != _id2last) {
	  _id1last = ic1;
	  _id2last = ic2;
	  unsigned int neu1 = RPV_helper::neutralinoIndex(ic1);
	  unsigned int neu2 = RPV_helper::neutralinoIndex(ic2);
	  _leftlast = 0.5*( (*_theN)(neu1, 3)*conj((*_theN)(neu2, 3)) -
			    (*_theN)(neu1, 2)*conj((*_theN)(neu2, 2)) );
	  if(_theN->size().first>4) {
	    for(unsigned int k=0;k<3;++k)
	      _leftlast -= 0.5*(*_theN)(neu1, 4+k)*conj((*_theN)(neu2, 4+k));
	  }
	  _rightlast = -conj(_leftlast);
	  _leftlast  /= _sw*_cw;
	  _rightlast /= _sw*_cw;
	}
      }
    }
  }
  norm ( _couplast);
  left ( _leftlast);
  right(_rightlast);
}
#line 1 "./RPVFFWVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVFFWVertex class.
//

#include "RPVFFWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "RPVhelper.h"

using namespace Herwig;

RPVFFWVertex::RPVFFWVertex() : _diagonal(false), _ckm(3,vector<Complex>(3,0.0)),
			       _sw(0.), _couplast(0.), _q2last(ZERO), 
			       _id1last(0), _id2last(0), _leftlast(0.),
			       _rightlast(0.), _interactions(0) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

IBPtr RPVFFWVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVFFWVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVFFWVertex::doinit() {
  // SUSY mixing matrices
  tSusyBasePtr model = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "RPVFFWVertex::doinit() - The model pointer is null!"
			  << Exception::abortnow;
  _theN = model->neutralinoMix();
  _theU = model->charginoUMix();
  _theV = model->charginoVMix();
  if(!_theN || !_theU || ! _theV)
    throw InitException() << "RPVFFWVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " N: " << _theN << " U: " << _theU << " V: "
			  << _theV << Exception::abortnow;
  // SM interactions
  if(_interactions==0 || _interactions==1) {
    // particles for outgoing W-
    // quarks
    for(int ix=1;ix<6;ix+=2) {
      for(int iy=2;iy<7;iy+=2) {
	bool isOff = iy/2 != (ix+1)/2;
	if ( isOff && _diagonal )
	  continue;
	addToList(-ix, iy, -24);
      }
    }
    // leptons
    for(int ix=11;ix<17;ix+=2) {
      int inu = model->majoranaNeutrinos() ? (ix+23)/2 : ix+1;
      addToList(-ix, inu, -24);
    }
    // particles for outgoing W+
    // quarks
    for(int ix=2;ix<7;ix+=2) {
      for(int iy=1;iy<6;iy+=2) {
	bool isOff = ix/2 != (iy+1)/2;
	if ( isOff && _diagonal )
	  continue;
	addToList(-ix, iy, 24);
      }
    }
    // leptons
    for(int ix=11;ix<17;ix+=2) {
      int inu = model->majoranaNeutrinos() ? (ix+23)/2 : -ix-1;
      addToList(inu, ix, 24);
    }
  }
  // neutralino and chargino
  if(_interactions==0 || _interactions==2) {
    vector<long> neu(4);
    neu[0] = 1000022; neu[1] = 1000023;
    neu[2] = 1000025; neu[3] = 1000035;
    if(_theN->size().first==7) {
      if(model->majoranaNeutrinos()) {
	neu.push_back(17);
	neu.push_back(18);
	neu.push_back(19);
      }
      else {
	neu.push_back(12);
	neu.push_back(14);
	neu.push_back(16);
      }
    }
    vector<long> cha(2);
    cha[0] = 1000024; cha[1] = 1000037;
    if(_theV->size().first==5) {
      cha.push_back(-11);
      cha.push_back(-13);
      cha.push_back(-15);
    }
    // sign == -1 outgoing W-, sign == +1 outgoing W+
    for(int sign = -1; sign < 2; sign += 2) {
      for(unsigned int ine = 0; ine < neu.size(); ++ine) {
	for(unsigned int ic = 0; ic < cha.size(); ++ic ) {
	  if(ic>1&&ine>3&&ic==ine-2) continue;
	  addToList(-sign*cha[ic], neu[ine], sign*24);
	}
      }
    }
  }
  Helicity::FFVVertex::doinit();
  // CKM matric
  if ( !_diagonal ) {
    Ptr<CKMBase>::transient_pointer CKM = model->CKM();
    // cast the CKM object to the HERWIG one
    ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer 
      hwCKM = ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardCKM>::
					transient_const_pointer>(CKM);
    if(hwCKM) {
      vector< vector<Complex > > CKM;
      CKM = hwCKM->getUnsquaredMatrix(generator()->standardModel()->families());
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  _ckm[ix][iy]=CKM[ix][iy];
	}
      }
    }
    else {
      throw Exception() << "Must have access to the Herwig::StandardCKM object"
			<< "for the CKM matrix in RPVFFWVertex::doinit()"
			<< Exception::runerror;
    }
  }
  _sw = sqrt(sin2ThetaW());
}

void RPVFFWVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _theN << _theU << _theV  << _diagonal << _ckm << _interactions;
}

void RPVFFWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _theN >> _theU >> _theV  >> _diagonal >> _ckm >> _interactions;
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<RPVFFWVertex,Helicity::FFVVertex>
  describeHerwigRPVFFWVertex("Herwig::RPVFFWVertex", "HwSusy.so HwRPV.so");

void RPVFFWVertex::Init() {

  static ClassDocumentation<RPVFFWVertex> documentation
    ("The couplings of the fermions to the W boson in the RPV model"
     " with bilinear R-parity violation");

  static Switch<RPVFFWVertex,unsigned int> interfaceInteractions
    ("Interactions",
     "Which interactions to include",
     &RPVFFWVertex::_interactions, 0, false, false);
  static SwitchOption interfaceInteractionsAll
    (interfaceInteractions,
     "All",
     "Include all the interactions",
     0);
  static SwitchOption interfaceInteractionsSM
    (interfaceInteractions,
     "SM",
     "Only include the MS terms",
     1);
  static SwitchOption interfaceInteractionsSUSY
    (interfaceInteractions,
     "SUSY",
     "Include the neutralino/chargino terms",
     2);

  static Switch<RPVFFWVertex,bool> interfaceDiagonal
    ("Diagonal",
     "Use a diagonal CKM matrix (ignoring the CKM object of the StandardModel).",
     &RPVFFWVertex::_diagonal, false, false, false);
  static SwitchOption interfaceDiagonalYes
    (interfaceDiagonal,
     "Yes",
     "Use a diagonal CKM matrix.",
     true);
  static SwitchOption interfaceDiagonalNo
    (interfaceDiagonal,
     "No",
     "Use the CKM object as used by the StandardModel.",
     false);

}

void RPVFFWVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			       tcPDPtr part2,
#ifndef NDEBUG
  tcPDPtr part3) {
#else
  tcPDPtr) {
#endif
  assert(abs(part3->id()) == ParticleID::Wplus);
  // normalization
  // first the overall normalisation
  if(q2 != _q2last||_couplast==0.) {
    _couplast = weakCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // left and right couplings for quarks
  if(abs(part1->id()) <= 6) {
    int iferm=abs(part1->id());
    int ianti=abs(part2->id());
    if(iferm%2!=0) swap(iferm,ianti);
    iferm = iferm/2;
    ianti = (ianti+1)/2;
    assert( iferm>=1 && iferm<=3 && ianti>=1 && ianti<=3);
    left(-sqrt(0.5)*_ckm[iferm-1][ianti-1]);
    right(0.);
  }
  else {
    long neu, cha;
    if(part1->charged()) {
      cha = part1->id();
      neu = part2->id();
    }
    else {
      cha = part2->id();
      neu = part1->id();
    }
    if(_theV->size().first==2&&abs(neu)<=16) {
      left(-sqrt(0.5));
      right(0.);
    }
    else {
      if(cha != _id1last || neu != _id2last) {
        _id1last = cha;
        _id2last = neu;
	unsigned int eigc = RPV_helper::charginoIndex(cha);
	unsigned int eign = RPV_helper::neutralinoIndex(neu);
	_leftlast = (*_theN)(eign, 1)*conj((*_theV)(eigc, 0)) - 
	  ( (*_theN)(eign, 3)*conj((*_theV)(eigc, 1))/sqrt(2));
	_rightlast = conj((*_theN)(eign, 1))*(*_theU)(eigc, 0) +
	  ( conj((*_theN)(eign, 2))*(*_theU)(eigc, 1)/sqrt(2));
	if(_theV->size().first==5) {
	  for(unsigned int k=0;k<3;++k)
	    _rightlast += ( conj((*_theN)(eign, 4+k))*(*_theU)(eigc, 2+k)/sqrt(2));
	}
      }
      Complex ltemp = _leftlast;
      Complex rtemp = _rightlast;
      bool chapart = abs(cha)>1000000 ? cha>0 : cha<0;
      // conjugate if +ve chargino
      if(chapart) {
	ltemp = conj(ltemp);
	rtemp = conj(rtemp);
      }
      if((part1->id()==cha&&chapart)||(part2->id()==cha&&!chapart)) {
	Complex temp = ltemp;
	ltemp  = -rtemp;
	rtemp = -temp;
      }
      left (ltemp);
      right(rtemp);
    }
  }
}
#line 1 "./RPVWSSVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVWSSVertex class.
//

#include "RPVWSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "RPV.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

RPVWSSVertex::RPVWSSVertex() :_sw(0.), _cw(0.), _s2w(0.), _c2w(0.),
			      _interactions(0), _q2last(), 
			      _ulast(0), _dlast(0), _gblast(0),
			      _factlast(0.), _couplast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}				 

IBPtr RPVWSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVWSSVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVWSSVertex::doinit() {
  // extract the model 
  tRPVPtr model = dynamic_ptr_cast<tRPVPtr>(generator()->standardModel());
  if( !model ) throw InitException() << "RPVWSSVertex::doinit() - The"
				     << " pointer to the RPV object is null!"
				     << Exception::abortnow;
  // get the mixing matrices
  MixingMatrixPtr mixH = model->CPevenHiggsMix() ;
  MixingMatrixPtr mixP = model->CPoddHiggsMix()  ;
  MixingMatrixPtr mixC = model->ChargedHiggsMix();
  // find the codes for the Higgs bosons
  vector<long> pseudo(1,36);
  vector<long> scalar(2);
  vector<long> charged(1,37);
  scalar[0] = 25;
  scalar[1] = 35;
  if(mixH&&mixH->size().first>2) {
    scalar.push_back(1000012);
    scalar.push_back(1000014);
    scalar.push_back(1000016);
  }
  if(mixP&&mixP->size().first>1) {
    pseudo.push_back(1000017);
    pseudo.push_back(1000018);
    pseudo.push_back(1000019);
  }
  if(mixC&&mixC->size().first>2) {
    charged.push_back(-1000011);
    charged.push_back(-1000013);
    charged.push_back(-1000015);
    charged.push_back(-2000011);
    charged.push_back(-2000013);
    charged.push_back(-2000015);
  }
  // sfermion interactions
  if(_interactions==0||_interactions==1) {
    // squarks
    //W-
    //LL-squarks
    for(long ix=1000001;ix<1000006;ix+=2) {
      addToList(-24,ix+1,-ix);
    }
    //1-2 stop sbottom
    addToList(-24,1000006,-2000005);
    //2-1 stop sbottom
    addToList(-24,2000006,-1000005);
    //2-2 stop sbottom
    addToList(-24,2000006,-2000005);
    //W+
    for(long ix=1000001;ix<1000006;ix+=2) {
      addToList(24,-(ix+1),ix);
    }
    //1-2 stop sbottom
    addToList(24,-1000006,2000005);
    //2-1 stop sbottom
    addToList(24,-2000006,1000005);
    //2-2 stop sbottom
    addToList(24,-2000006,2000005);
    //LL squarks
    for(long ix=1000001;ix<1000007;++ix) {
      addToList(23,ix,-ix);
    }
    //RR squarks
    for(long ix=2000001;ix<2000007;++ix) {
      addToList(23,ix,-ix);
    }
    //L-Rbar stop
    addToList(23,1000006,-2000006);
    //Lbar-R stop
    addToList(23,-1000006,2000006);
    //L-Rbar sbottom
    addToList(23,1000005,-2000005);
    //Lbar-R sbottom
    addToList(23,-1000005,2000005);
    // gamma
    //squarks
    for(long ix=1000001;ix<1000007;++ix) {
      addToList(22,ix,-ix);
    }
    for(long ix=2000001;ix<2000007;++ix) {
      addToList(22,ix,-ix);
    }
    //sleptons
    // gamma
    for(long ix=1000011;ix<1000016;ix+=2) {
      addToList(22,ix,-ix);
    }
    for(long ix=2000011;ix<2000016;ix+=2) {
      addToList(22,ix,-ix);
    }
    // Z
    //LL-sleptons
    for(long ix=1000011;ix<1000017;ix+=2) {
      addToList(23,ix,-ix);
    }
    // RR-sleptons
    for(long ix=2000011;ix<2000016;ix+=2) {
      addToList(23,ix,-ix);
    }
    //L-Rbar stau
    addToList(23,1000015,-2000015);
    //Lbar-R stau
    addToList(23,-1000015,2000015);
    if(!(mixH&&mixH->size().first>2)) {
      for(long ix=1000012;ix<1000017;ix+=2) {
	// sneutrinos
	addToList(23,ix,-ix);
      }
      //LL-sleptons
      for(long ix=1000011;ix<1000016;ix+=2) {
	addToList(-24,-ix,ix+1);
      }
      //2-L stau
      addToList(-24,-2000015,1000016);
      //LL-sleptons
      for(long ix=1000011;ix<1000016;ix+=2) {
	addToList(24,ix,-ix-1);
      }
      //2-L stau
      addToList(24,2000015,-1000016);
    }
  }
  if(_interactions==0||_interactions==2) {
    // charged Higgs and photon
    addToList(22,37,-37);
    // charged Higgs and Z
    for(unsigned int ix=0;ix<charged.size();++ix) {
      for(unsigned int iy=0;iy<charged.size();++iy) {
	if(abs(charged[ix])>1000000&&abs(charged[iy])>100000 && 
	   ( charged[ix]==charged[iy] ||
	     (abs(charged[ix])%1000000==15&&abs(charged[iy])%1000000==15)))
	  continue;
	addToList(23,charged[ix],-charged[iy]);
      }
    }
    // neutral Higgs and Z
    for(unsigned int ix=0;ix<scalar.size();++ix) {
      for(unsigned int iy=0;iy<pseudo.size();++iy) {
	addToList(23,scalar[ix],pseudo[iy]);
      }
    }
    // charged higss, scalar higgs and W
    for(unsigned int ix=0;ix<charged.size();++ix) {
      for(unsigned int iy=0;iy<scalar.size();++iy) {
	addToList( 24,-charged[ix],scalar[iy]);
	addToList(-24, charged[ix],scalar[iy]);
      }
    }
    // charged higss, pseudoscalar higgs and W
    for(unsigned int ix=0;ix<charged.size();++ix) {
      for(unsigned int iy=0;iy<pseudo.size();++iy) {
	addToList( 24,-charged[ix],pseudo[iy]);
	addToList(-24, charged[ix],pseudo[iy]);
      }
      for(unsigned int iy=0;iy<scalar.size();++iy) {
	addToList( 24,-charged[ix],scalar[iy]);
	addToList(-24, charged[ix],scalar[iy]);
      }
    }
  }
  VSSVertex::doinit();
  // sfermion mixing
  _stop = model->stopMix();
  _sbottom = model->sbottomMix();
  if(!_stop || !_sbottom)
    throw InitException() << "RPVWSSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << _stop << " sbottom: " << _sbottom
			  << Exception::abortnow;
  _stau = model->stauMix();
  if(!_stau && (!mixC || mixC->size().first<2))
    throw InitException() << "RPVWSSVertex::doinit() either the stau"
			  << " mixing matrix must be set or the stau"
			  << " included in mixing with the"
			  << " charged Higgs bosons" << Exception::abortnow;
  // weak mixing
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt( 1. - sin2ThetaW() );
  _s2w = 2.*_cw*_sw;
  _c2w = 1. - 2.*sin2ThetaW();

  // Coupling of Z to scalar and pseudoscalar
  Cijeo_.resize(pseudo.size(),vector<Complex>(scalar.size(),0.));
  for(unsigned int i=0;i<pseudo.size();++i) {
    for(unsigned int j=0;j<scalar.size();++j) {
      for(unsigned int ix=0;ix<mixH->size().second;++ix) {
	double sign = ix!=1 ? 1. : -1.;
	Cijeo_[i][j] += sign*(*mixP)(i,ix)*(*mixH)(j,ix);
      }
    }
  }
  // Coupling of W to scalar charged
  Cijec_.resize(scalar.size(),vector<Complex>(charged.size(),0.));
  for(unsigned int i=0;i<scalar.size();++i) {
    for(unsigned int j=0;j<charged.size();++j) {
      for(unsigned int ix=0;ix<mixH->size().second;++ix) {
	double sign = ix!=1 ? 1. : -1.;
	Cijec_[i][j] += sign*(*mixH)(i,ix)*(*mixC)(j,ix);
      }
    }
  }
  // Coupling of W to pseudopseudo charged
  Cijco_.resize(pseudo.size(),vector<Complex>(charged.size(),0.));
  for(unsigned int i=0;i<pseudo.size();++i) {
    for(unsigned int j=0;j<charged.size();++j) {
      for(unsigned int ix=0;ix<mixP->size().second;++ix) {
	// not sure about this need it to get agreement with SPheno
	//double sign = ix!=1 ? 1. : -1.;
	double sign = 1.;
	Cijco_[i][j] += sign*(*mixP)(i,ix)*(*mixC)(j,ix);
      }
    }
  }
  // Coupling of Z to charged Higgs
  Cijc_.resize(charged.size(),vector<Complex>(charged.size(),0.));
  for(unsigned int i=0;i<charged.size();++i) {
    for(unsigned int j=0;j<charged.size();++j) {
      for(unsigned int ix=5;ix<mixC->size().second;++ix) {
	Cijc_[i][j] += (*mixC)(i,ix)*(*mixC)(j,ix);
      }
    }
  }
}

void RPVWSSVertex::persistentOutput(PersistentOStream & os) const {
  os << _interactions << _sw  << _cw
     << _stau << _stop << _sbottom << _s2w << _c2w
     << Cijeo_ << Cijec_ << Cijco_ << Cijc_;
}

void RPVWSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _interactions >> _sw >> _cw
     >> _stau >> _stop >> _sbottom >> _s2w >> _c2w
     >> Cijeo_ >> Cijec_ >> Cijco_ >> Cijc_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RPVWSSVertex,Helicity::VSSVertex>
describeHerwigRPVWSSVertex("Herwig::RPVWSSVertex", "HwSusy.so HwRPV.so");

void RPVWSSVertex::Init() {

  static ClassDocumentation<RPVWSSVertex> documentation
    ("There is no documentation for the RPVWSSVertex class");

  static Switch<RPVWSSVertex,unsigned int> interfaceInteractions
    ("Interactions",
     "Which interactions to include",
     &RPVWSSVertex::_interactions, 0, false, false);
  static SwitchOption interfaceInteractionsAll
    (interfaceInteractions,
     "All",
     "Include both the interactions which would have been sfermion"
     " and Higgs bosons with the gauge bosons in the MSSM",
     0);
  static SwitchOption interfaceInteractionsSfermions
    (interfaceInteractions,
     "Sfermions",
     "Include the sfermion interactions",
     1);
  static SwitchOption interfaceInteractionsHiggs
    (interfaceInteractions,
     "Higgs",
     "Include the Higgs boson interactions",
     2);

}


void RPVWSSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  long gboson = part1->id();
  assert(     gboson  == ParticleID::Z0    ||
	      gboson  == ParticleID::gamma || 
	  abs(gboson) == ParticleID::Wplus );
  long h1ID = part2->id();
  long h2ID = part3->id();
  // squarks and sleptons
  if(((abs(h1ID)>=1000001&&abs(h1ID)<=1000006)||
      (abs(h1ID)>=2000001&&abs(h1ID)<=2000006)) ||
     (((abs(h1ID)>=1000011&&abs(h1ID)<=1000016)||
       (abs(h1ID)>=2000011&&abs(h1ID)<=2000016))&&
      Cijeo_.size()==1)) {
    long sf1(abs(part2->id())),sf2(abs(part3->id()));
    if( sf1 % 2 != 0 ) swap(sf1, sf2);
    if( sf1 != _ulast || sf2 != _dlast || gboson != _gblast) {
      _gblast = gboson;
      _ulast = sf1;
      _dlast = sf2;
      //photon is simplest
      if( gboson == ParticleID::gamma )
	_factlast = getParticleData(sf1)->charge()/eplus;
      else {
	//determine which helicity state
	unsigned int alpha(sf1/1000000 - 1), beta(sf2/1000000 - 1);
	//mixing factors
	Complex m1a(0.), m1b(0.);
	if( sf1 == ParticleID::SUSY_t_1 || sf1 == ParticleID::SUSY_t_2 )
	  m1a = (*_stop)(alpha, 0);
	else if( sf1 == ParticleID::SUSY_b_1 || sf1 == ParticleID::SUSY_b_2 )
	  m1a = (*_sbottom)(alpha, 0);
	else if( sf1 == ParticleID::SUSY_tau_1minus || 
		 sf1 == ParticleID::SUSY_tau_2minus )
	  m1a = (*_stau)(alpha, 0);
	else
	  m1a = (alpha == 0) ? Complex(1.) : Complex(0.);
	if( sf2 == ParticleID::SUSY_t_1 || sf2 == ParticleID::SUSY_t_2 )
	  m1b = (*_stop)(beta, 0);
	else if( sf2 == ParticleID::SUSY_b_1 || sf2 == ParticleID::SUSY_b_2 )
	  m1b = (*_sbottom)(beta, 0);
	else if( sf2 == ParticleID::SUSY_tau_1minus || 
		 sf2 == ParticleID::SUSY_tau_2minus )
	  m1b = (*_stau)(beta, 0);
	else
	  m1b = (beta == 0) ? Complex(1.) : Complex(0.);
	//W boson
	if( abs(gboson) == ParticleID::Wplus ) {
	  _factlast = m1a*m1b/sqrt(2)/_sw;
	}
	//Z boson
	else {
	  if( sf1 == ParticleID::SUSY_nu_eL || sf1 == ParticleID::SUSY_nu_muL ||
	      sf1 == ParticleID::SUSY_nu_tauL ) {
	    _factlast = 1./_cw/2./_sw;
	  }
	  else {
	    double lmda(1.);
	    if( sf2 % 2 == 0 ) lmda = -1.;
	    _factlast = lmda*m1a*m1b;
	    if( alpha == beta) {
	      double ef = getParticleData(sf1)->charge()/eplus;
	      _factlast += 2.*ef*sqr(_sw);
	    }
	    _factlast *= -1./2./_cw/_sw; 
	  }
	}
      }
    }
  }
  // Higgs bosons
  else {
    _gblast = gboson;
    _ulast = h1ID;
    _dlast = h2ID;
    _factlast = 0.;
    if( gboson == ParticleID::Z0 ) {
      if( part2->charged() ) {
	unsigned int c1 = abs(h1ID) < 1000000 ? 0 : 
	  (abs(h1ID) < 2000000 ? (abs(h1ID)-1000009)/2 : (abs(h1ID)-2000003)/2);
	unsigned int c2 = abs(h2ID) < 1000000 ? 0 : 
	  (abs(h2ID) < 2000000 ? (abs(h2ID)-1000009)/2 : (abs(h2ID)-2000003)/2);
	if(c1==c2) _factlast = (_c2w-Cijc_[c1][c2])/_s2w;
	else       _factlast = -Cijc_[c1][c2]/_s2w;
	if(part2->iCharge()<0) _factlast *= -1.;
      }
      else {
	if(h1ID == ParticleID::h0         || h1ID  == ParticleID::H0         ||
	   h1ID == ParticleID::SUSY_nu_eL || h1ID == ParticleID::SUSY_nu_muL ||
	   h1ID == ParticleID::SUSY_nu_tauL ) {
	  unsigned int is = h1ID < 1000000 ? (h1ID-25)/10 : (h1ID-1000008)/2;
	  unsigned int ip = h2ID < 1000000 ? 0 : (h2ID-1000016);
	  _factlast =  Complex(0.,1.)*Cijeo_[ip][is]/_s2w;
	}
	else {
	  unsigned int is = h2ID < 1000000 ? (h2ID-25)/10 : (h2ID-1000008)/2;
	  unsigned int ip = h1ID < 1000000 ? 0 : (h1ID-1000016);
	  _factlast = -Complex(0.,1.)*Cijeo_[ip][is]/_s2w;
	}

      }
    }
    else if( gboson == ParticleID::gamma ) {
      _factlast = part2->iCharge()/3;
    }
    else {
      long scalar  = part2->charged() ? h2ID : h1ID;
      long charged = part2->charged() ? h1ID : h2ID;
      unsigned int ic = abs(charged) < 1000000 ? 0 : 
	(abs(charged) < 2000000 ? (abs(charged)-1000009)/2 : (abs(charged)-2000003)/2);
      if(scalar == ParticleID::h0         || scalar  == ParticleID::H0         ||
	 scalar == ParticleID::SUSY_nu_eL || scalar == ParticleID::SUSY_nu_muL ||
	 scalar == ParticleID::SUSY_nu_tauL ) {
	unsigned int ih = scalar < 1000000 ? (scalar-25)/10 : (scalar-1000008)/2;
	_factlast = -0.5*Cijec_[ih][ic]/_sw;
	if(gboson<0) _factlast *= -1.;
      }
      else {
	unsigned int ih = scalar < 1000000 ? 0 : (scalar-1000016);
	_factlast =  Complex(0., 0.5)*Cijco_[ih][ic]/_sw;
      } 
      if(part3->charged()) _factlast *= -1.;
    }
  }
  if( q2 != _q2last || _couplast==0. ) {
    _q2last = q2;
    _couplast = electroMagneticCoupling(q2);
  }
  if(part2->id()>0)
    norm(-_couplast*_factlast);
  else
    norm(+_couplast*_factlast);
}
#line 1 "./RPVFFSVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVFFSVertex class.
//

#include "RPVFFSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

namespace {

  unsigned int neutralinoIndex(long id) {
    if(id> 1000000)
      return id<1000025 ? id-1000022 : (id-1000005)/10;
    else if(abs(id)<=16) 
      return (abs(id)-4)/2;
    else
      return id-13;
  }
  
  unsigned int charginoIndex(long id) {
    return abs(id)>1000000 ? (abs(id)-1000024)/13 : (abs(id)-7)/2;
  }

}

RPVFFSVertex::RPVFFSVertex() : interactions_(0), mw_(ZERO),
			       _q2last(ZERO), _couplast(0.),
			       _leftlast(0.),_rightlast(0.),
			       _id1last(0), _id2last(0), _id3last(0),
			       yukawa_(true),
			       _sw(0.), _cw(0.),_sb(0.), _cb(0.),
			       vd_(ZERO), vu_(ZERO),
			       _massLast(make_pair(ZERO,ZERO)) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}  

IBPtr RPVFFSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVFFSVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVFFSVertex::doinit() {
  // cast the model to the RPV one
  model_ = dynamic_ptr_cast<tRPVPtr>(generator()->standardModel());
  if( !model_ ) throw InitException() << "RPVFFSVertex::doinit() - "
				      << "The pointer to the MSSM object is null!"
				      << Exception::abortnow;
  // get the various mixing matrices
  _stop = model_->stopMix();
  _sbot = model_->sbottomMix();
  _stau = model_->stauMix();
  _nmix = model_->neutralinoMix();
  _umix = model_->charginoUMix();
  _vmix = model_->charginoVMix();
  _mixH = model_->CPevenHiggsMix();
  _mixP = model_->CPoddHiggsMix();
  _mixC = model_->ChargedHiggsMix();
  if(!_stop || !_sbot ) throw InitException() << "RPVFFSVertex::doinit() - "
					      << "A squark mixing matrix pointer is null."
					      << " stop: " << _stop << " sbottom: "
					      << _sbot << Exception::abortnow;
  if(!_nmix) throw InitException() << "RPVFFSVertex::doinit() - "
				   << "A gaugino mixing matrix pointer is null."
				   << " N: " << _nmix << " U: " << _umix 
				   << " V: " << _vmix << Exception::abortnow;
  if(! ( _stau || (!_stau&& _mixC->size().first>1))) 
    throw InitException() << "RPVFFSVertex::doinit() - "
			  << "Must have either the stau mixing matrix or it"
			  << " must be part of the charged Higgs boson mixing."
			  << Exception::abortnow;
  // various interactions
  // scalar Higgs bosons
  vector<long> h0(2);
  h0[0] = 25; h0[1] = 35; 
  if(_mixH&&_mixH->size().first>2) {
    h0.push_back(1000012); h0.push_back(1000014); h0.push_back(1000016);
  }
  // pseudoscalar Higgs bosons
  vector<long> a0(1,36);
  if(_mixP&&_mixP->size().first>1) {
    a0.push_back(1000017); a0.push_back(1000018); a0.push_back(1000019);
  }
  // charged Higgs bosons
  vector<long> hp(1,37);
  if(_mixC->size().first>1) {
    hp.push_back(-1000011); hp.push_back(-1000013); hp.push_back(-1000015);
    hp.push_back(-2000011); hp.push_back(-2000013); hp.push_back(-2000015);
  }
  // neutralinos
  vector<long> neu(4);
  neu[0] = 1000022; neu[1] = 1000023;
  neu[2] = 1000025; neu[3] = 1000035;
  if(_nmix->size().first>4) {
    if(model_->majoranaNeutrinos()) {
      neu.push_back(17);
      neu.push_back(18);
      neu.push_back(19);
    }
    else {
      neu.push_back(12);
      neu.push_back(14);
      neu.push_back(16);
    }
  }
  // charginos
  vector<long> chg(2);
  chg[0] = 1000024; chg[1] = 1000037;
  if(_umix->size().first>2) {
    chg.push_back(-11); chg.push_back(-13); chg.push_back(-15);
  }
  // FFH
  if(interactions_==0||interactions_==1) {
    // quarks neutral scalar
    for ( unsigned int h = 0; h < h0.size(); ++h ) {
      for(long ix=1;ix<7;++ix) addToList(-ix,ix,h0[h]);
    } 
    // quarks neutral pseudoscalar
    for ( unsigned int h = 0; h < a0.size(); ++h ) {
      for(long ix=1;ix<7;++ix) addToList(-ix,ix,a0[h]);
    }
    // quarks charged higgs
    for(unsigned int h=0;h<hp.size();++h) {
      for(long ix=1;ix<6;ix+=2) {
	//outgoing H+
	addToList(-ix-1,  ix, hp[h]);
	//outgoing H-
	addToList(-ix  ,ix+1,-hp[h]);
      }
    }
    // charged leptons neutral scalar
    for ( unsigned int h = 0; h < h0.size(); ++h ) {
      for(long ix=11;ix<16;ix+=2) addToList(-ix,ix,h0[h]);
    } 
    // charged leptons neutral pseudoscalar
    for ( unsigned int h = 0; h < a0.size(); ++h ) {
      for(long ix=11;ix<16;ix+=2) addToList(-ix,ix,a0[h]);
    }
    // charged higgs to leptons, no mixing
    if(_nmix->size().first<=4) {
      for(unsigned int h=0;h<hp.size();++h) {
	for(long ix=11;ix<16;ix+=2) {
	  //outgoing H+
	  addToList(-ix-1,  ix, hp[h]);
	  //outgoing H-
	  addToList(-ix  ,ix+1,-hp[h]);
	}
      }
    }
  }
  // GOGOH
  if(interactions_==0||interactions_==2) {
    // neutral scalar Higgs neutralino
    for(unsigned int i=0;i<h0.size();++i) {
      for(unsigned int j = 0; j < neu.size(); ++j) {
	for(unsigned int k = j; k < neu.size(); ++k) {
	  addToList(neu[j], neu[k], h0[i]);
	}
      }
    }
    // neutral pseudoscalar Higgs neutralino
    for(unsigned int i=0;i<a0.size();++i) {
      for(unsigned int j = 0; j < neu.size(); ++j) {
	for(unsigned int k = j; k < neu.size(); ++k) {
	  addToList(neu[j], neu[k], a0[i]);
	}
      }
    }
    // neutral scalar Higgs chargino
    for(unsigned int i=0;i<h0.size();++i) {
      for(unsigned int j = 0; j < chg.size(); ++j) {
	for(unsigned int k = 0; k < chg.size(); ++k) {
	  if(j==k&&abs(chg[j])<16) continue;
	  addToList(-chg[j], chg[k], h0[i]);
	}
      }
    }
    // neutral scalar Higgs chargino
    for(unsigned int i=0;i<a0.size();++i) {
      for(unsigned int j = 0; j < chg.size(); ++j) {
	for(unsigned int k = 0; k < chg.size(); ++k) {
	  if(j==k&&abs(chg[j])<16) continue;
	  addToList(-chg[j], chg[k], a0[i]);
	}
      }
    }
    // charged higgs
    for(unsigned int i=0;i<hp.size();++i) {
      for(unsigned int j = 0; j < neu.size(); ++j) {
	for(unsigned int k = 0; k < chg.size(); ++k) {
 	  addToList(-chg[k], neu[j], hp[i]);
 	  addToList( chg[k], neu[j],-hp[i]);
	}
      }
    }
  }
  // neutralino sfermion
  if(interactions_==0||interactions_==3) {
    // quarks
    for(unsigned int nl = 0; nl < neu.size(); ++nl) {
      for(long ix=1;ix<7;++ix){
	addToList( neu[nl],  ix, -(1000000+ix) );
	addToList( neu[nl],  ix, -(2000000+ix) );
	addToList( neu[nl], -ix,  (1000000+ix) );
	addToList( neu[nl], -ix,  (2000000+ix) );
      }
    }
    // neutrino
    if(_nmix->size().first<=4) {
      for(unsigned int nl = 0; nl < neu.size(); ++nl) {
	for(long ix=12;ix<17;ix+=2) {
	  addToList( neu[nl],  ix, -(1000000+ix) );
	  addToList( neu[nl], -ix,  (1000000+ix) );
	}
      }
    }
    // charged leptons no mixing
    if(!_mixC || _mixC->size().first==1) {
      for(unsigned int nl = 0; nl < neu.size(); ++nl) {
	for(long ix=11;ix<17;ix+=2) {
	  addToList( neu[nl],  ix, -(1000000+ix) );
	  addToList( neu[nl], -ix,  (1000000+ix) );
	  addToList( neu[nl],  ix, -(2000000+ix) );
	  addToList( neu[nl], -ix,  (2000000+ix) );
	}
      }
    }
  }
  // chargino sfermion
  if(interactions_==0||interactions_==4) {
    //quarks 
    for(unsigned int ic = 0; ic < chg.size(); ++ic) {
      for(long ix = 1; ix < 7; ++ix) {
	if( ix % 2 == 0 ) {
	  addToList(-chg[ic], ix,-( 999999+ix));
	  addToList(-chg[ic], ix,-(1999999+ix));
	  addToList( chg[ic],-ix,   999999+ix );
	  addToList( chg[ic],-ix,  1999999+ix );
	}
	else {
	  addToList(-chg[ic],-ix,  1000001+ix );
	  addToList(-chg[ic],-ix,  2000001+ix );
	  addToList( chg[ic], ix,-(1000001+ix));
	  addToList( chg[ic], ix,-(2000001+ix));
	}
      }
    }
    // sneutrinos
    if(!_mixH || _mixH->size().first<=2) {
      for(unsigned int ic = 0; ic < chg.size(); ++ic) {
	for(long ix = 11; ix < 17; ix+=2) {
	  addToList(-chg[ic],-ix,1000001+ix);
	  addToList(chg[ic],ix,-(1000001+ix));
	}
      }
    }
    // charged leptons
    if(!_mixC || _mixC->size().first==1) {
      for(unsigned int ic = 0; ic < chg.size(); ++ic) {
	for(long ix = 12; ix < 17; ix+=2) {
	  addToList(-chg[ic], ix,-( 999999+ix));
	  addToList(-chg[ic], ix,-(1999999+ix));
	  addToList( chg[ic],-ix, ( 999999+ix));
	  addToList( chg[ic],-ix, (1999999+ix));
	}
      }
    }
  }
  FFSVertex::doinit();
  // various couplings and parameters
  mw_ = getParticleData(ParticleID::Wplus)->mass();
  _sw = sqrt(sin2ThetaW());
  double tb = model_->tanBeta();
  _cw = sqrt(1. - sqr(_sw));
  _sb = tb/sqrt(1 + sqr(tb));
  _cb = sqrt(1 - sqr(_sb));
  vector<Energy> vnu = model_->sneutrinoVEVs();
  double g = electroMagneticCoupling(sqr(mw_))/_sw;
  Energy v = 2.*mw_/g;
  vd_ = sqrt((sqr(v)-sqr(vnu[0])-sqr(vnu[1])-sqr(vnu[2]))/
		   (1.+sqr(tb)));
  vu_ = vd_*tb;
  // couplings of the neutral scalar Higgs to charginos
  Energy me   = model_->mass(sqr(mw_),getParticleData(ParticleID::eminus  ));
  Energy mmu  = model_->mass(sqr(mw_),getParticleData(ParticleID::muminus ));
  Energy mtau = model_->mass(sqr(mw_),getParticleData(ParticleID::tauminus));
  double h_E[3] = {sqrt(2.)*me/vd_/g,sqrt(2.)*mmu /vd_/g,sqrt(2.)*mtau/vd_/g};
  for(unsigned int ih=0;ih<_mixH->size().first;++ih ) {
    OCCHL_.push_back(vector<vector<Complex> >
		     (_umix->size().first,vector<Complex>(_umix->size().first,0.)));
    for(unsigned int i=0;i<_umix->size().first;++i) {
      for(unsigned int j=0;j<_umix->size().first;++j) {
	OCCHL_[ih][i][j]    = -sqrt(0.5)*
	  ( (*_mixH)(ih,0)*(*_vmix)(i,0)*(*_umix)(j,1) +
	    (*_mixH)(ih,1)*(*_vmix)(i,1)*(*_umix)(j,0));
	for(unsigned int k=2;k<_umix->size().first;++k) {
	  OCCHL_[ih][i][j] -= sqrt(0.5)*(+(*_mixH)(ih,k)*(*_vmix)(i,0)*(*_umix)(j,k)
					 +h_E[k-2]*(*_umix)(j,k)*(*_vmix)(i,k)*(*_mixH)(ih,0)
					 -h_E[k-2]*(*_umix)(j,1)*(*_vmix)(i,k)*(*_mixH)(ih,k));
	}
      }
    }
  }
  // couplings of the neutral scalar Higgs to neutralinos
  double tw = _sw/_cw;
  for(unsigned int ih=0;ih<_mixH->size().first;++ih) {
    ONNHL_.push_back(vector<vector<Complex> >
		     (_nmix->size().first,vector<Complex>(_nmix->size().first,0.)));
    for(unsigned int i=0;i<_nmix->size().first;++i) {
      for(unsigned int j=0;j<_nmix->size().first;++j) {
	ONNHL_[ih][i][j] =0.;
	for(unsigned int in=2;in<_nmix->size().first;++in) {
	  double sign = in!=3 ? 1. : -1.;
	  ONNHL_[ih][i][j] += 0.5*sign*(*_mixH)(ih,in-2)*
	    ( + (tw*(*_nmix)(i,0) - (*_nmix)(i,1) )*(*_nmix)(j,in)
	      + (tw*(*_nmix)(j,0) - (*_nmix)(j,1) )*(*_nmix)(i,in) );
	} 
      }
    }
  }
  // couplings of the neutral pseudoscalar Higgs to neutralinos
  for(unsigned int ih=0;ih<_mixP->size().first;++ih) {
    ONNAL_.push_back(vector<vector<Complex> >
		     (_nmix->size().first,vector<Complex>(_nmix->size().first,0.)));
    for(unsigned int i=0;i<_nmix->size().first;++i) {
      for(unsigned int j=0;j<_nmix->size().first;++j) {
	ONNAL_[ih][i][j] =0.;
	for(unsigned int in=2;in<_nmix->size().first;++in) {
	  double sign = in!=3 ? 1. : -1.;
	  ONNAL_[ih][i][j] += -0.5*sign*(*_mixP)(ih,in-2)*
	    ( + (tw*(*_nmix)(i,0) - (*_nmix)(i,1) )*(*_nmix)(j,in)
	      + (tw*(*_nmix)(j,0) - (*_nmix)(j,1) )*(*_nmix)(i,in) );
	}
      }
    }
  }
  // couplings of the neutral pseudoscalar higgs to charginos
  for(unsigned int ih=0;ih<_mixP->size().first;++ih) {
    OCCAL_.push_back(vector<vector<Complex> >
		     (_umix->size().first,vector<Complex>(_umix->size().first,0.)));
    for(unsigned int i=0;i<_umix->size().first;++i) {
      for(unsigned int j=0;j<_umix->size().first;++j) {
	OCCAL_[ih][i][j] = 
	  (*_mixP)(ih,0)*(*_vmix)(i,0)*(*_umix)(j,1)+
	  (*_mixP)(ih,1)*(*_vmix)(i,1)*(*_umix)(j,0);
	for(unsigned int k=2;k<_umix->size().first;++k) {
	  OCCAL_[ih][i][j] += (*_mixP)(ih,k)*(*_vmix)(i,0)*(*_umix)(j,k)
	    -h_E[k-2]*(*_umix)(j,k)*(*_vmix)(i,k)*(*_mixP)(ih,0)
	    +h_E[k-2]*(*_umix)(j,1)*(*_vmix)(i,k)*(*_mixP)(ih,k);
	}
	OCCAL_[ih][i][j] *= sqrt(0.5);
      }
    }
  }
  // couplings for the charged higgs
  for(unsigned int ih=0;ih<_mixC->size().first;++ih) {
    OCNSL_.push_back(vector<vector<Complex> >
		     (_nmix->size().first,vector<Complex>(_umix->size().first,0.)));
    OCNSR_.push_back(vector<vector<Complex> >
		     (_nmix->size().first,vector<Complex>(_umix->size().first,0.)));
    for(unsigned int i = 0; i < _nmix->size().first; ++i) {
      for(unsigned int j=0;j<_umix->size().first;++j) {
	OCNSL_[ih][i][j] = (*_mixC)(ih,1)*conj((*_nmix)(i, 3)*(*_vmix)(j,0)
					       +((*_nmix)(i,1) + (*_nmix)(i,0)*tw)*(*_vmix)(j,1)/sqrt(2));


	OCNSR_[ih][i][j] = (*_mixC)(ih,0)*    ((*_nmix)(i, 2)*(*_umix)(j,0)
					       -((*_nmix)(i,1) + (*_nmix)(i,0)*tw)*(*_umix)(j,1)/sqrt(2));
	for(unsigned int k=2;k<_umix->size().first;++k) {
	  OCNSL_[ih][i][j] += -h_E[k-2]*(*_mixC)(ih,0)*conj((*_nmix)(i,2+k)*(*_vmix)(j,k))
	    +(*_mixC)(ih,k  )*h_E[k-2]*conj((*_nmix)(i, 2)*(*_vmix)(j,k))
	    +(*_mixC)(ih,k+3)*tw*sqrt(2.)*conj((*_nmix)(i,0)*(*_vmix)(j,k));
	  OCNSR_[ih][i][j] += (*_mixC)(ih,k)*((*_nmix)(i,2+k)*(*_umix)(j,0)
					      -((*_nmix)(i,1) + (*_nmix)(i,0)*tw)*(*_umix)(j,k)/sqrt(2))
	    -(*_mixC)(ih,k+3)*h_E[k-2]*((*_nmix)(i,k+2)*(*_umix)(j,1)-(*_nmix)(i,2)*(*_umix)(j,k));
	}
      }
    }
  }
}

void RPVFFSVertex::persistentOutput(PersistentOStream & os) const {
  os << interactions_ << _stop << _sbot << _stau << _umix << _vmix
     << _nmix << _mixH << _mixP << _mixC << ounit(mw_,GeV) << yukawa_
     << model_  << _sw << _cw << _sb << _cb << ounit(vd_,GeV) << ounit(vu_,GeV)
     << OCCHL_ << ONNHL_ << ONNAL_ << OCCAL_ << OCNSL_ << OCNSR_;
}

void RPVFFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> interactions_ >> _stop >> _sbot >> _stau >> _umix >> _vmix
     >> _nmix >> _mixH >> _mixP >> _mixC >> iunit(mw_,GeV) >> yukawa_
     >> model_ >> _sw >> _cw  >> _sb >> _cb >> iunit(vd_,GeV) >> iunit(vu_,GeV)
     >> OCCHL_ >> ONNHL_ >> ONNAL_ >> OCCAL_ >> OCNSL_ >> OCNSR_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RPVFFSVertex,Helicity::FFSVertex>
describeHerwigRPVFFSVertex("Herwig::RPVFFSVertex", "HwSusy.so HwRPV.so");

void RPVFFSVertex::Init() {

  static ClassDocumentation<RPVFFSVertex> documentation
    ("The RPVFFSVertex class implements all the couplings of fermion-antiferion to scalars"
     " in R-Parity violating models, including sfermion fermion gaugino, SM ferimon antiferimon Higgs "
     "and gaugino-gaugino Higgs.");

  static Switch<RPVFFSVertex,unsigned int> interfaceInteractions
    ("Interactions",
     "Whice interactions to include",
     &RPVFFSVertex::interactions_, 0, false, false);
  static SwitchOption interfaceInteractionsAll
    (interfaceInteractions,
     "All",
     "Include all the interactions",
     0);
  static SwitchOption interfaceInteractionsHiggsSMFermions
    (interfaceInteractions,
     "HiggsSMFermions",
     "Interactions of Higgs with SM fermions",
     1);
  static SwitchOption interfaceInteractionsHiggsGaugino
    (interfaceInteractions,
     "HiggsGaugino",
     "Interactions of the Higgs with the gauginos",
     2);
  static SwitchOption interfaceInteractionsNeutralinoSfermion
    (interfaceInteractions,
     "NeutralinoSfermion",
     "Include the neutralino sfermion interactions",
     3);
  static SwitchOption interfaceInteractionsCharginoSfermions
    (interfaceInteractions,
     "CharginoSfermions",
     "Include the chargino sfermion interactions",
     4);

  static Switch<RPVFFSVertex,bool> interfaceYukawa
    ("Yukawa",
     "Whether or not to include the Yukawa type couplings in neutralino/chargino interactions",
     &RPVFFSVertex::yukawa_, true, false, false);
  static SwitchOption interfaceYukawaYes
    (interfaceYukawa,
     "Yes",
     "Include the terms",
     true);
  static SwitchOption interfaceYukawaNo
    (interfaceYukawa,
     "No",
     "Don't include them",
     false);
}

void RPVFFSVertex::neutralinoSfermionCoupling(Energy2 q2, tcPDPtr fermion,
					      tcPDPtr gaugino, tcPDPtr sfermion) {
  long ism(abs(fermion->id())),ig(abs(gaugino->id())),isc(sfermion->id());
  if( ig != _id1last || ism != _id2last || isc != _id3last ) {
    _id1last = ig;
    _id2last = ism;
    _id3last = isc;
    // sfermion mass eigenstate
    unsigned int alpha(abs(isc)/1000000 - 1);
    // neutralino state
    unsigned int nl = neutralinoIndex(ig);
    assert(nl<=6);
    // common primed neutralino matrices
    Complex n2prime = (*_nmix)(nl,1)*_cw - (*_nmix)(nl,0)*_sw;
    //handle neutrinos first
    if( ism == 12 || ism == 14 || ism == 16 ) {
      _leftlast = Complex(0., 0.);
      _rightlast = -sqrt(0.5)*n2prime/_cw;
    }
    else {
      Complex n1prime = (*_nmix)(nl,0)*_cw + (*_nmix)(nl,1)*_sw;
      tcPDPtr smf = getParticleData(ism);
      double qf = smf->charge()/eplus;
      Complex bracketl = qf*_sw*( conj(n1prime) - _sw*conj(n2prime)/_cw );
      double y = yukawa_ ? double(model_->mass(q2, smf)/2./mw_) : 0.;
      double lambda(0.);
      //neutralino mixing element
      Complex nlf(0.);
      if( ism % 2 == 0 ) {
	y /= _sb;
	lambda = -0.5 + qf*sqr(_sw);
	nlf = (*_nmix)(nl,3);
      }
      else { 
	y /= _cb;
	lambda = 0.5 + qf*sqr(_sw);
	nlf = (*_nmix)(nl,2);
      }
      Complex bracketr = _sw*qf*n1prime - n2prime*lambda/_cw;
      //heavy quarks/sleptons
      if( ism == 5 || ism == 6 || ism == 15 ) {
	Complex ma1(0.), ma2(0.);
	if( ism == 5 ) {
	  ma1 = (*_sbot)(alpha,0);
	  ma2 = (*_sbot)(alpha,1);
	} 
	else if( ism == 6 ) {
	  ma1 = (*_stop)(alpha,0);
	  ma2 = (*_stop)(alpha,1);
	} 
	else {
	  ma1 = (*_stau)(alpha,0);
	  ma2 = (*_stau)(alpha,1);
	}
	_leftlast = y*conj(nlf)*ma1 - ma2*bracketl;
	_rightlast = y*nlf*ma2 + ma1*bracketr;
      }
      else {
	if( alpha == 0 ) {
	  _leftlast = y*conj(nlf);
	  _rightlast = bracketr;
	} 
	else {
	  _leftlast = -bracketl;
	  _rightlast = y*nlf;
	}
      }
      _leftlast  *= -sqrt(2.);
      _rightlast *= -sqrt(2.);
    }
  }
  //determine the helicity order of the vertex
  if( fermion->id() < 0 ) {
    left (conj(_rightlast));
    right(conj( _leftlast));
  }
  else {
    left ( _leftlast);
    right(_rightlast);
  }
  norm(_couplast);
}

void RPVFFSVertex::charginoSfermionCoupling(Energy2 q2, tcPDPtr fermion,
					    tcPDPtr gaugino, tcPDPtr sfermion) {
  long ism(abs(fermion->id())),ig(abs(gaugino->id())),isc(sfermion->id());
  if( ig != _id1last || ism != _id2last || isc != _id3last ) {
    _id1last = ig;
    _id2last = ism;
    _id3last = isc;
    // sfermion mass eigenstate
    unsigned int alpha(abs(isc)/1000000 - 1);
    // get the type of chargino
    unsigned int ch = charginoIndex(ig);
    assert(ch<=4);
    // various mixing matrix elements
    Complex ul1 = (*_umix)(ch,0), ul2 = (*_umix)(ch,1);
    Complex vl1 = (*_vmix)(ch,0), vl2 = (*_vmix)(ch,1);
    // lepton/slepton
    if( ism >= 11 && ism <= 16 ) {
      long lept = ( ism % 2 == 0 ) ? ism - 1 : ism;
      double y = yukawa_ ? 
	double(model_->mass(q2, getParticleData(lept))/mw_/sqrt(2)/_cb) : 0.;
      if( ism == 12 || ism == 14 ) {
	_leftlast = 0.;
	_rightlast = alpha == 0 ? - ul1 : y*ul2;
      }
      else if( ism == 16 ) {
	_leftlast = 0.;
	_rightlast = -ul1*(*_stau)(alpha, 0) + y*(*_stau)(alpha,1)*ul2;
      }
      else if( ism == 11 || ism == 13 || ism == 15 ) {
	_leftlast = y*conj(ul2);
	_rightlast = -vl1;
      }
    }
    // squark/quark
    else {
      double yd(0.), yu(0.);
      if(yukawa_) {
	if( ism % 2 == 0) {
	  yu = model_->mass(q2, getParticleData(ism))/mw_/sqrt(2)/_sb;
	  yd = model_->mass(q2, getParticleData(ism - 1))/mw_/sqrt(2)/_cb;
	}
	else {
	  yu = model_->mass(q2, getParticleData(ism + 1))/mw_/sqrt(2)/_sb;
	  yd = model_->mass(q2, getParticleData(ism))/mw_/sqrt(2)/_cb;
	}
      }
      //heavy quarks
      if( ism == 5 ) {
	_leftlast =  yd*conj(ul2)*(*_stop)(alpha,0);
	_rightlast = -vl1*(*_stop)(alpha, 0) + yu*vl2*(*_stop)(alpha,1);
      }
      else if( ism == 6 ) {
	_leftlast =  yu*conj(vl2)*(*_sbot)(alpha,0);
	_rightlast = -ul1*(*_sbot)(alpha, 0) + yd*ul2*(*_sbot)(alpha,1);
      }
      else {
	if( alpha == 0 ) {
	  _leftlast  = (ism % 2 == 0) ? yu*conj(vl2) : yd*conj(ul2);
	  _rightlast = (ism % 2 == 0) ? -ul1 : -vl1;
	}
	else {
	  _leftlast = 0.;
	  _rightlast = (ism % 2 == 0) ? yd*ul2 : yu*vl2;
	}
      }
    }
  }
  //determine the helicity order of the vertex
  if( fermion->id() < 0 ) {
    left (conj(_rightlast));
    right(conj( _leftlast));
  }
  else {
    left ( _leftlast);
    right(_rightlast);
  }
  norm(_couplast);
}

void RPVFFSVertex::higgsFermionCoupling(Energy2 q2, tcPDPtr f1,
					tcPDPtr f2, tcPDPtr higgs) {
  long f1ID(f1->id()), f2ID(f2->id()), isc(higgs->id());
  // running fermion masses
  if( q2 != _q2last || _id1last  != f1ID) {
    _massLast.first  = model_->mass(q2,f1);
    _id1last  = f1ID;
  }
  if( q2 != _q2last || _id2last != f2ID) {
    _massLast.second = model_->mass(q2,f2);
    _id2last = f2ID;
  }
  if( q2 != _q2last) _id3last = isc;
  Complex fact(0.);
  // scalar neutral Higgs
  if(isc == ParticleID::h0         || isc  == ParticleID::H0         ||
     isc == ParticleID::SUSY_nu_eL || isc == ParticleID::SUSY_nu_muL ||
     isc == ParticleID::SUSY_nu_tauL ) {
    unsigned int ih = isc < 1000000 ? (isc-25)/10 : (isc-1000008)/2;
    unsigned int id = abs(f1ID);
    fact = -Complex(_massLast.first*
		    ((id%2==0) ? (*_mixH)(ih,1)/vu_ : (*_mixH)(ih,0)/vd_));
    left (1.);
    right(1.);
  }
  // pseudoscalar neutral Higgs
  else if(isc == ParticleID::A0 || isc == 1000017 || isc == 1000018 ||
	  isc == 1000019 ) {
    unsigned int ih = isc < 1000000 ? 0 : (isc-1000016);
    unsigned int id = abs(f1ID);
    if(_mixP) {
      fact = Complex(-Complex(0., 1.)*_massLast.first*
		     ( (id%2==0) ?  (*_mixP)(ih,1)/vu_ : (*_mixP)(ih,0)/vd_));
    }
    else {
      fact = Complex(-Complex(0., 1.)*_massLast.first*
		     ( (id%2==0) ?  _cb/vu_ : _sb/vd_));
    }
    left ( 1.);
    right(-1.);
  }
  // charged higgs
  else {
    if(!_mixC) {
      if( abs(f1ID) % 2 == 0 ) {
	_leftlast  =  _massLast.first /vu_*_cb;
	_rightlast =  _massLast.second/vd_*_sb;
      }
      else {
	_leftlast  =  _massLast.second/vu_*_cb;
	_rightlast =  _massLast.first /vd_*_sb;
      }
    }
    else {
      unsigned int ih;
      if(abs(isc)==ParticleID::Hplus) {
	ih = 0;
      }
      else {
	isc *= -1;
	ih = abs(isc)<2000000 ? (abs(isc)-1000009)/2 : (abs(isc)-2000003)/2;
      }
      if( abs(f1ID) % 2 == 0 ) {
	_leftlast  =  Complex(_massLast.first /vu_*(*_mixC)(ih,1));
	_rightlast =  Complex(_massLast.second/vd_*(*_mixC)(ih,0));
      }
      else {
	_leftlast  =  Complex(_massLast.second/vu_*(*_mixC)(ih,1));
	_rightlast =  Complex(_massLast.first /vd_*(*_mixC)(ih,0));
      }
    }
    if( isc > 0 ) swap(_leftlast,_rightlast);
    fact = sqrt(2.);
    left ( _leftlast);
    right(_rightlast);
  }
  norm(fact);
}

void RPVFFSVertex::higgsGauginoCoupling(Energy2, tcPDPtr f1,
					tcPDPtr f2, tcPDPtr higgs) {
  long f1ID(f1->id()), f2ID(f2->id()), isc(higgs->id());
  if( isc == _id3last && f1ID == _id1last && f2ID == _id2last ) {
    left ( _leftlast);
    right(_rightlast);
  }
  else {
    _id1last = f1ID;
    _id2last = f2ID;
    _id3last = isc;
    // scalar neutral Higgs
    if(isc == ParticleID::h0         || isc  == ParticleID::H0         ||
       isc == ParticleID::SUSY_nu_eL || isc == ParticleID::SUSY_nu_muL ||
       isc == ParticleID::SUSY_nu_tauL ) {
      unsigned int ih = isc < 1000000 ? (isc-25)/10 : (isc-1000008)/2;
      // charginos
      if(f1->charged()) {
	unsigned int ei = charginoIndex(f1ID);
	unsigned int ej = charginoIndex(f2ID);
	if     (ei< 2&&f1ID>0) swap(ei,ej);
	else if(ei>=2&&f1ID<0) swap(ei,ej);
	_rightlast  = conj(OCCHL_[ih][ej][ei]);
	_leftlast   =      OCCHL_[ih][ei][ej] ;
      }
      // neutralinos
      else {
	unsigned int ei = neutralinoIndex(f1ID);
	unsigned int ej = neutralinoIndex(f2ID);
	_leftlast  = conj(ONNHL_[ih][ej][ei]);
	_rightlast =      ONNHL_[ih][ei][ej] ;
      }
    }
    // pseudoscalar neutral Higgs
    else if(isc == ParticleID::A0 || isc == 1000017 || isc == 1000018 ||
	    isc == 1000019 ) {
      unsigned int ih = isc < 1000000 ? 0 : (isc-1000016);
      // charginos
      if(f1->charged()) {
	unsigned int ei = charginoIndex(f1ID);
	unsigned int ej = charginoIndex(f2ID);
	if     (ei< 2&&f1ID>0) swap(ei,ej);
	else if(ei>=2&&f1ID<0) swap(ei,ej);
	_rightlast = -Complex(0.,1.)*conj(OCCAL_[ih][ej][ei]);
	_leftlast  =  Complex(0.,1.)*     OCCAL_[ih][ei][ej] ;
      }
      // neutralinos
      else {
	unsigned int ei = neutralinoIndex(f1ID);
	unsigned int ej = neutralinoIndex(f2ID);
	_leftlast  =  Complex(0.,1.)*conj(ONNAL_[ih][ej][ei]);
	_rightlast = -Complex(0.,1.)*     ONNAL_[ih][ei][ej] ;
      }
    }
    // charged higgs
    else {
      unsigned int ih = abs(isc) < 1000000 ? 0 : 
	(abs(isc) < 2000000 ? (abs(isc)-1000009)/2 : (abs(isc)-2000003)/2);
      long chg(f2ID), neu(f1ID);
      if(f1->charged()) swap(chg, neu);
      unsigned int ei = neutralinoIndex(neu);
      unsigned int ej = charginoIndex(chg);
      _leftlast  = -OCNSL_[ih][ei][ej];
      _rightlast = -OCNSR_[ih][ei][ej];
      bool chargedSwap = abs(isc)<1000000 ? isc<0 : isc>0;
      if( chargedSwap ) {
	Complex tmp = _leftlast;
	_leftlast  = conj(_rightlast);
	_rightlast = conj(tmp);
      }
    }
    left ( _leftlast);
    right(_rightlast);
  }
  norm(_couplast);
}

void RPVFFSVertex::setCoupling(Energy2 q2, tcPDPtr part1, 
			       tcPDPtr part2,tcPDPtr part3) {
  // overall normalisation
  if(q2!=_q2last || _couplast==0.) {
    _couplast = weakCoupling(q2);
    _q2last=q2;
  }
  long f1ID(part1->id()), f2ID(part2->id()), isc(abs(part3->id()));
  // squark quark
  if(part3->coloured()) {
    tcPDPtr smfermion = part1, gaugino = part2;
    if(gaugino->coloured()) swap(smfermion,gaugino);
    if(gaugino->charged())
      charginoSfermionCoupling(q2,smfermion,gaugino,part3);
    else
      neutralinoSfermionCoupling(q2,smfermion,gaugino,part3);
  }
  // slepton/lepton without mixing
  else if((( isc >= 1000011 && isc <= 1000016) ||
	   ( isc >= 2000011 && isc <= 2000016)) && 
	  (!_mixC || _mixC->size().first<=1 || 
	   !_mixP || _mixP->size().first<=1 )) {
    tcPDPtr smfermion = part1, gaugino = part2;
    if(abs(gaugino->id())<1000000) swap(smfermion,gaugino);
    if(gaugino->charged())
      charginoSfermionCoupling(q2,smfermion,gaugino,part3);
    else
      neutralinoSfermionCoupling(q2,smfermion,gaugino,part3);
  }
  // SM quarks and Higgs
  else if((abs(f1ID) <=  6 && abs(f2ID) <=  6) ||
	  ((abs(f1ID) >= 11 && abs(f1ID) <= 16) &&
	   (abs(f2ID) >= 11 && abs(f2ID) <= 16) && 
	   _umix->size().first==2) ) {
    higgsFermionCoupling(q2,part1,part2,part3);
  }
  // gauginos and the Higgs (general case for sleptons)
  else {
    higgsGauginoCoupling(q2,part1,part2,part3);
  }
}
#line 1 "./RPVWWHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVWWHVertex class.
//

#include "RPVWWHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "RPV.h"

using namespace Herwig;

RPVWWHVertex::RPVWWHVertex() : coupLast_(ZERO), q2Last_(ZERO) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void RPVWWHVertex::doinit() {
  // extract some models parameters to decide if we need sneutrinos
  tRPVPtr model = dynamic_ptr_cast<tRPVPtr>(generator()->standardModel());
  if( !model ) throw InitException() << "RPVWWHVertex::doinit() - "
				     << "The pointer to the RPV object is null!"
				     << Exception::abortnow;
  // get the Higgs mixing matrix
  MixingMatrixPtr mix = model->CPevenHiggsMix();
  // possible interactions
  vector<long> higgs(2);
  higgs[0] = 25; higgs[1] = 35;
  if(mix->size().first>2) {
    higgs.push_back(1000012);
    higgs.push_back(1000014);
    higgs.push_back(1000016);
  }
  for(unsigned int ix=0;ix<higgs.size();++ix) {
    addToList( 23, 23,higgs[ix]);
    addToList(-24, 24,higgs[ix]);
  }
  VVSVertex::doinit();
  // SM parameters
  Energy mw = getParticleData(ParticleID::Wplus)->mass();
  Energy mz = getParticleData(ParticleID::Z0)->mass();
  double sw = sqrt(sin2ThetaW());
  double cw = sqrt(1.-sin2ThetaW());
  vector<Energy> vnu = model->sneutrinoVEVs();
  Energy v = 2.*mw/electroMagneticCoupling(sqr(mw))*sw;
  double tanb = model->tanBeta();
  Energy vd = sqrt((sqr(v)-sqr(vnu[0])-sqr(vnu[1])-sqr(vnu[2]))/
		   (1.+sqr(tanb)));
  Energy vu = vd*tanb;
  for(unsigned int ix=0;ix<higgs.size();++ix) {
    complex<Energy> c = vd*(*mix)(ix,0)+vu*(*mix)(ix,1);
    for(size_t iy=2; iy<mix->size().second; ++iy) c += vnu[iy-2]*(*mix)(ix,iy);
    vector<complex<Energy> > coup(2);
    coup[0] = c/v*mw;
    coup[1] = c/v*mz/cw;
    couplings_.push_back(coup);
  }
}

IBPtr RPVWWHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVWWHVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(couplings_,GeV);
}

void RPVWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(couplings_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RPVWWHVertex,Helicity::VVSVertex>
describeHerwigRPVWWHVertex("Herwig::RPVWWHVertex", "HwSusy.so HwRPV.so");

void RPVWWHVertex::Init() {

  static ClassDocumentation<RPVWWHVertex> documentation
    ("The RPVWWHVertex class implements the couplings of a pair of electroweak"
     " gauge bosons to the higgs boson in he R-parity violating MSSM.");

}

void RPVWWHVertex::setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr,
			      tcPDPtr particle3) {
  long bosonID = abs(particle1->id());
  long higgsID =     particle3->id();
  assert( bosonID == ParticleID::Wplus || bosonID == ParticleID::Z0 );
  int ihiggs = higgsID>1000000 ? (higgsID-1000008)/2 : (higgsID-25)/10;
  assert(ihiggs>=0 && ihiggs<=4);
  complex<Energy> fact = bosonID==ParticleID::Wplus ? 
    couplings_[ihiggs][0] : couplings_[ihiggs][1];
  if( q2 != q2Last_ ) {
    q2Last_ = q2;
    coupLast_ = weakCoupling(q2);
  }
  norm(coupLast_*fact*UnitRemoval::InvE);
}
#line 1 "./RPVSSSVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVSSSVertex class.
//

#include "RPVSSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "RPV.h"
#include <cassert>

using namespace Herwig;
using namespace ThePEG::Helicity;

RPVSSSVertex::RPVSSSVertex() : interactions_(0), q2Last_(ZERO), gLast_(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr RPVSSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVSSSVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVSSSVertex::persistentOutput(PersistentOStream & os) const {
  os << interactions_ << stau_ << sbottom_ << stop_ << ounit(scalarScalarScalar_,GeV)
     << ounit(scalarPseudoPseudo_,GeV) << ounit(scalarSup_,GeV) << ounit(scalarSdown_,GeV)
     << ounit(scalarSneutrino_,GeV) << ounit(scalarSlepton_,GeV)
     << ounit(chargedSquark_,GeV) << ounit(chargedSlepton_,GeV)
     << ounit(scalarChargedCharged_,GeV) << ounit(pseudoChargedCharged_,GeV)
     << ounit(pseudoSup_,GeV) << ounit(pseudoSdown_,GeV) << ounit(pseudoSlepton_,GeV);
}

void RPVSSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> interactions_ >> stau_ >> sbottom_ >> stop_ >> iunit(scalarScalarScalar_,GeV)
     >> iunit(scalarPseudoPseudo_,GeV) >> iunit(scalarSup_,GeV) >> iunit(scalarSdown_,GeV)
     >> iunit(scalarSneutrino_,GeV) >> iunit(scalarSlepton_,GeV)
     >> iunit(chargedSquark_,GeV) >> iunit(chargedSlepton_,GeV)
     >> iunit(scalarChargedCharged_,GeV) >> iunit(pseudoChargedCharged_,GeV)
     >> iunit(pseudoSup_,GeV) >> iunit(pseudoSdown_,GeV) >> iunit(pseudoSlepton_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RPVSSSVertex,SSSVertex>
describeHerwigRPVSSSVertex("Herwig::RPVSSSVertex", "HwRPV.so");

void RPVSSSVertex::Init() {

  static ClassDocumentation<RPVSSSVertex> documentation
    ("The RPVSSSVertex class implements the coupling of three"
     " scalar particles in the RPV model");

  static Switch<RPVSSSVertex,unsigned int> interfaceInteractions
    ("Interactions",
     "Which interactions to include",
     &RPVSSSVertex::interactions_, 0, false, false);
  static SwitchOption interfaceInteractionsAll
    (interfaceInteractions,
     "All",
     "Include both triple Higgs and Higgs sfermion interactions",
     0);
  static SwitchOption interfaceInteractionsHiggsHiggsHiggs
    (interfaceInteractions,
     "HiggsHiggsHiggs",
     "Only include triple Higgs boson interactions",
     1);
  static SwitchOption interfaceInteractionsHiggsSfermions
    (interfaceInteractions,
     "HiggsSfermions",
     "Only include Higgs sfermion interactions",
     2);

}

void  RPVSSSVertex::doinit() {
  // extract the model 
  tRPVPtr model = dynamic_ptr_cast<tRPVPtr>(generator()->standardModel());
  if( !model ) throw InitException() << "RPVSSSVertex::doinit() - The"
				     << " pointer to the RPV object is null!"
				     << Exception::abortnow;
  // get the Higgs mixing matrices
  MixingMatrixPtr mixH = model->CPevenHiggsMix() ;
  MixingMatrixPtr mixP = model->CPoddHiggsMix()  ;
  MixingMatrixPtr mixC = model->ChargedHiggsMix();
  // find the codes for the Higgs bosons
  vector<long> pseudo(1,36);
  vector<long> scalar(2);
  vector<long> charged(1,37);
  scalar[0] = 25;
  scalar[1] = 35;
  if(mixH&&mixH->size().first>2) {
    scalar.push_back(1000012);
    scalar.push_back(1000014);
    scalar.push_back(1000016);
  }
  if(mixP&&mixP->size().first>1) {
    pseudo.push_back(1000017);
    pseudo.push_back(1000018);
    pseudo.push_back(1000019);
  }
  if(mixC&&mixC->size().first>2) {
    charged.push_back(-1000011);
    charged.push_back(-1000013);
    charged.push_back(-1000015);
    charged.push_back(-2000011);
    charged.push_back(-2000013);
    charged.push_back(-2000015);
  }
  // triple Higgs interactions
  if(interactions_==0||interactions_==1) {
    // three neutral scalar bosons
    for(unsigned int i1=0;i1<scalar.size();++i1) {
      for(unsigned int i2=0;i2<=i1;++i2) {
	for(unsigned int i3=0;i3<=i2;++i3) {
	  addToList(scalar[i1],scalar[i2],scalar[i3]);
	}
      }
    }
    // one neutral scalar two pseudoscalars
    for(unsigned int i1=0;i1<scalar.size();++i1) {
      for(unsigned int i2=0;i2<pseudo.size();++i2) {
	for(unsigned int i3=0;i3<=i2;++i3) {
	  addToList(scalar[i1],pseudo[i2],pseudo[i3]);
	}
      }
    }
    // one neutral scalar two charged scalars
    for(unsigned int i1=0;i1<scalar.size();++i1) {
      for(unsigned int i2=0;i2<charged.size();++i2) {
	for(unsigned int i3=0;i3<charged.size();++i3) {
	  if(!(abs(charged[i2])>1000000&&abs(charged[i3])>1000000&&
	     abs(charged[i2])%1000000==abs(charged[i3])%1000000))
	    addToList(scalar[i1],charged[i2],-charged[i3]);
	}
      }
    }
    // one pseudo scalar two charged scalars
    if(charged.size()>1) {
      for(unsigned int i1=0;i1<pseudo.size();++i1) {
	for(unsigned int i2=0;i2<charged.size();++i2) {
	  for(unsigned int i3=0;i3<charged.size();++i3) {
	    if(i2==i3) continue;
	    if(!(abs(charged[i2])>1000000&&abs(charged[i3])>1000000&&
		 abs(charged[i2])%1000000==abs(charged[i3])%1000000))
	      addToList(pseudo[i1],charged[i2],-charged[i3]);
	  }
	}
      }
    }
  }
  // sfermion interactions
  if(interactions_==0||interactions_==2) {
    // scalar neutral higgs sfermion sfermion 
    for(unsigned int i = 0; i < scalar.size(); ++i) {
      // squarks
      for(unsigned int j = 1; j < 7; ++j) {
	long lj = 1000000 + j;
	long rj = 2000000 + j;
	//LLbar
	addToList(scalar[i],lj,-lj);
	//RRbar
	addToList(scalar[i],rj,-rj);
	//LRbar
	addToList(scalar[i],lj,-rj);
	//RLbar
	addToList(scalar[i],rj,-lj);
      }
      // sleptons
      for(unsigned int j = 11; j < 17; ++j) {
	long lj = 1000000 + j;
	long rj = 2000000 + j;
	if( j % 2 != 0) {
	  // LL
	  addToList(scalar[i],lj,-lj);
	  // RR
	  addToList(scalar[i],rj,-rj);
	  //LRbar
	  addToList(scalar[i],lj,-rj);
	  //RLbar
	  addToList(scalar[i],rj,-lj);
	}
	// sneutrino only if no mixing
	else if(scalar.size()==2) {
	  addToList(scalar[i],lj,-lj);
	}
      }
    }
    // pseudoscalar sfermion sfermion
    for(unsigned int i = 0; i < pseudo.size(); ++i) {
      // squarks
      for(unsigned int j = 1; j < 7; ++j) {
	long lj = 1000000 + j;
	long rj = 2000000 + j;
	//LRbar
	addToList(pseudo[i],lj,-rj);
	//RLbar
	addToList(pseudo[i],rj,-lj);
      }
      // sleptons
      if(scalar.size()==2) {
	for(unsigned int j = 11; j < 17; j += 2) {
	  long lj = 1000000 + j;
	  long rj = 2000000 + j;
	  addToList(36,lj,-rj);
	  addToList(36,rj,-lj);
	}
      }
    }
    //outgoing H+
    for(unsigned int i=0;i<charged.size();++i) {
      // squarks
      for(long ii = 2; ii < 7; ii += 2) {
	//LL
	addToList(charged[i], 999999 + ii, -1000000 - ii);
	//RR
	addToList(charged[i], 1999999 + ii, -2000000 - ii);
	//RL
	addToList(charged[i], 1999999 + ii, -1000000 - ii);
	//LR
	addToList(charged[i], 999999 + ii, -2000000 - ii);
      }
      if(scalar.size()==2) {
	for(long ii = 11; ii < 17; ii += 2) {
	  addToList(37, 1000000 + ii, -1000001 - ii);
	  addToList(37, 2000000 + ii, -1000001 - ii);
	}
      }
      //outgoing H-
      for(long ii = 2; ii < 7; ii += 2) {
	//LL
	addToList(-charged[i], 1000000 + ii, -999999 - ii);
	//RR
	addToList(-charged[i], 2000000 + ii, -1999999 - ii);
	//RL
	addToList(-charged[i], 1000000 + ii, -1999999 - ii);
	//LR
	addToList(-charged[i], 2000000 + ii, -999999 - ii);
      }
      if(scalar.size()==2) {
	for(long ii = 11; ii < 17; ii += 2) {
	  addToList(-37, 1000001 + ii, -1000000 - ii);
	  addToList(-37, 1000001 + ii, -2000000 - ii);
	}
      }
    }
  }
  SSSVertex::doinit();
  // extract the sfermion mixing matrices
  // sfermion mixing
  stop_ = model->stopMix();
  sbottom_ = model->sbottomMix();
  if(!stop_ || !sbottom_)
    throw InitException() << "RPVSSSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << stop_ << " sbottom: " << sbottom_
			  << Exception::abortnow;
  stau_ = model->stauMix();
  if(!stau_ && (!mixC || mixC->size().first<2))
    throw InitException() << "RPVSSSVertex::doinit() either the stau"
			  << " mixing matrix must be set or the stau"
			  << " included in mixing with the"
			  << " charged Higgs bosons" << Exception::abortnow;
  // couplings of the neutral Higgs bosons 
  Energy mw = getParticleData(ParticleID::Wplus)->mass();
  double sw = sqrt(sin2ThetaW());
  double cw = sqrt(1.-sin2ThetaW());
  // extract the vevs
  vector<Energy> vnu = model->sneutrinoVEVs();
  double g = electroMagneticCoupling(sqr(mw))/sw;
  Energy v = 2.*mw/g;
  double tanb = model->tanBeta();
  vector<Energy> vevs(5);
  vevs[0] = sqrt((sqr(v)-sqr(vnu[0])-sqr(vnu[1])-sqr(vnu[2]))/
		   (1.+sqr(tanb)));
  vevs[1] = vevs[0]*tanb;
  for(unsigned int ix=0;ix<vnu.size();++ix) vevs[2+ix] = vnu[ix];
  // bilinear RPV terms
  // coupling of three scalar Higgs bosons
  scalarScalarScalar_.resize(scalar.size(),vector<vector<complex<Energy > > >
			     (scalar.size(),vector<complex<Energy> >(scalar.size(),complex<Energy>(ZERO))));
  double delta[5]={1.,-1.,1.,1.,1.};
  if(vevs.size()>scalar.size()) vevs.resize(scalar.size());
  for(unsigned int i1=0;i1<scalar.size();++i1) {
    for(unsigned int i2=0;i2<scalar.size();++i2) {
      for(unsigned int i3=0;i3<scalar.size();++i3) {
	scalarScalarScalar_[i1][i2][i3] = ZERO;
	for(unsigned int ix=0;ix<vevs.size();++ix) {
	  for(unsigned int k=0;k<scalar.size();++k) {
	    scalarScalarScalar_[i1][i2][i3] += 
	      delta[ix]*vevs[ix]*(*mixH)(i1,ix)*delta[k]*(*mixH)(i2,k)*(*mixH)(i3,k)+
	      delta[ix]*vevs[ix]*(*mixH)(i2,ix)*delta[k]*(*mixH)(i1,k)*(*mixH)(i3,k)+
	      delta[ix]*vevs[ix]*(*mixH)(i3,ix)*delta[k]*(*mixH)(i1,k)*(*mixH)(i2,k);
	  }
	}
	scalarScalarScalar_[i1][i2][i3] *= -0.25*g/sqr(cw);
      }
    }
  }
  // coupling of scalar Higgs to 2 pseudoscalars 
  scalarPseudoPseudo_.resize(scalar.size(),vector<vector<complex<Energy > > >
			     (pseudo.size(),vector<complex<Energy> >(pseudo.size(),complex<Energy>(ZERO))));
  for(unsigned int i1=0;i1<scalar.size();++i1) {
    for(unsigned int i2=0;i2<pseudo.size();++i2) {
      for(unsigned int i3=0;i3<pseudo.size();++i3) {
	scalarPseudoPseudo_[i1][i2][i3] = ZERO;
	for(unsigned int ix=0;ix<vevs.size();++ix) {
	  for(unsigned int k=0;k<pseudo.size();++k) {
	    scalarPseudoPseudo_[i1][i2][i3] += 
	      delta[k]*delta[ix]*vevs[ix]*(*mixH)(i1,ix)*(*mixP)(i2,k)*(*mixP)(i3,k);
	  }
	}
	scalarPseudoPseudo_[i1][i2][i3] *= -0.25*g/sqr(cw);
      }
    }
  }
  // coupling of scalar Higgs to 2 charged Higgs
  double gp = g*sw/cw;
  Energy mu = model->muParameter();
  vector<Energy> eps = model->epsilon();
  scalarChargedCharged_.resize(scalar.size(),vector<vector<complex<Energy > > >
			       (charged.size(),vector<complex<Energy> >(charged.size(),complex<Energy>(ZERO))));
  double yl[3] = {sqrt(2.)*getParticleData(11)->mass()/vevs[1],
		  sqrt(2.)*getParticleData(13)->mass()/vevs[1],
		  sqrt(2.)*getParticleData(15)->mass()/vevs[1]};
  complex<Energy> Al[3] =  {complex<Energy>(ZERO),complex<Energy>(ZERO),yl[2]*model->tauTrilinear()};
  for(unsigned int i1=0;i1<scalar.size();++i1) {
    for(unsigned int i2=0;i2<charged.size();++i2) {
      for(unsigned int i3=0;i3<charged.size();++i3) {
	complex<Energy> umdi = vevs[0]*(*mixH)(i1,0)-vevs[1]*(*mixH)(i1,1)+vevs[2]*(*mixH)(i1,2)+
	  vevs[3]*(*mixH)(i1,3)+vevs[4]*(*mixH)(i1,4);
	scalarChargedCharged_[i1][i2][i3] =
	  // HH term
	  0.25*sqr(g)*(-vevs[1]*((*mixH)(i1,1)*((*mixC)(i2,1)*(*mixC)(i3,1)+(*mixC)(i2,2)*(*mixC)(i3,2))+
				 (*mixH)(i1,0)*((*mixC)(i2,2)*(*mixC)(i3,1)+(*mixC)(i2,1)*(*mixC)(i3,2)))
		       -vevs[0]*((*mixH)(i1,0)*((*mixC)(i2,1)*(*mixC)(i3,1)+(*mixC)(i2,2)*(*mixC)(i3,2))+
				 (*mixH)(i1,1)*((*mixC)(i2,2)*(*mixC)(i3,1)+(*mixC)(i2,1)*(*mixC)(i3,2)))
		       +(vevs[2]*(*mixH)(i1,2)+vevs[3]*(*mixH)(i1,3)+vevs[4]*(*mixH)(i1,4))*
		       ((*mixC)(i2,1)*(*mixC)(i3,1)-(*mixC)(i2,2)*(*mixC)(i3,2)))-
	  0.25*sqr(gp)*((*mixC)(i2,1)*(*mixC)(i3,1)-(*mixC)(i2,2)*(*mixC)(i3,2))*umdi
	  -0.5*(*mixC)(i2,1)*(*mixC)(i3,1)*(yl[0]*vevs[2]*(*mixH)(i1,2)+yl[1]*vevs[3]*(*mixH)(i1,3)+yl[2]*vevs[4]*(*mixH)(i1,4))
	  // LL term
	  +0.25*(sqr(g)-sqr(gp))*umdi*((*mixC)(i2,2)*(*mixC)(i3,2)+(*mixC)(i2,3)*(*mixC)(i3,3)+(*mixC)(i2,4)*(*mixC)(i3,4))
	  -vevs[0]*(*mixH)(i1,1)*(yl[0]*(*mixC)(i2,2)*(*mixC)(i3,2)+yl[1]*(*mixC)(i2,3)*(*mixC)(i3,3)+yl[2]*(*mixC)(i2,4)*(*mixC)(i3,4))
	  -0.25*sqr(g)*((vevs[2]*(*mixC)(i3,2)+vevs[3]*(*mixC)(i3,3)+vevs[4]*(*mixC)(i3,4))*
			((*mixH)(i1,2)*(*mixC)(i2,2)+(*mixH)(i1,3)*(*mixC)(i2,3)+(*mixH)(i1,4)*(*mixC)(i2,4))+
			(vevs[2]*(*mixC)(i2,2)+vevs[3]*(*mixC)(i2,3)+vevs[4]*(*mixC)(i2,4))*
			((*mixH)(i1,2)*(*mixC)(i3,2)+(*mixH)(i1,3)*(*mixC)(i3,3)+(*mixH)(i1,4)*(*mixC)(i3,4)))
	  // RR term
	  +0.5*sqr(gp)*umdi*((*mixC)(i2,5)*(*mixC)(i3,5)+(*mixC)(i2,6)*(*mixC)(i3,6)+(*mixC)(i2,7)*(*mixC)(i3,7))
	  -vevs[0]*(*mixH)(i1,1)*(yl[0]*(*mixC)(i2,5)*(*mixC)(i3,5)+yl[1]*(*mixC)(i2,6)*(*mixC)(i3,6)+yl[2]*(*mixC)(i2,7)*(*mixC)(i3,7))
	  -0.5*((vevs[2]*yl[0]*(*mixC)(i2,5)+vevs[3]*yl[1]*(*mixC)(i2,6)+vevs[4]*yl[2]*(*mixC)(i2,7))*
		(yl[0]*(*mixH)(i1,2)*(*mixC)(i3,5)+yl[1]*(*mixH)(i1,3)*(*mixC)(i3,6)+yl[2]*(*mixH)(i1,4)*(*mixC)(i3,7))+
		(vevs[2]*yl[0]*(*mixC)(i3,5)+vevs[3]*yl[1]*(*mixC)(i3,6)+vevs[4]*yl[2]*(*mixC)(i3,7))*
		(yl[0]*(*mixH)(i1,2)*(*mixC)(i2,5)+yl[1]*(*mixH)(i1,3)*(*mixC)(i2,6)+yl[2]*(*mixH)(i1,4)*(*mixC)(i2,7)))
	  // LR
	  -(*mixH)(i1,1)/sqrt(2.)*(Al[0]*(*mixC)(i2,2)*(*mixC)(i3,2)+Al[1]*(*mixC)(i2,3)*(*mixC)(i3,3)+Al[2]*(*mixC)(i2,4)*(*mixC)(i3,4))
	  +(*mixH)(i1,2)/sqrt(2.)*mu*(yl[0]*(*mixC)(i2,2)*(*mixC)(i3,2)+yl[1]*(*mixC)(i2,3)*(*mixC)(i3,3)+yl[2]*(*mixC)(i2,4)*(*mixC)(i3,4))
	  // RL
	  -(*mixH)(i1,1)/sqrt(2.)*(Al[0]*(*mixC)(i3,2)*(*mixC)(i2,2)+Al[1]*(*mixC)(i3,3)*(*mixC)(i2,3)+Al[2]*(*mixC)(i3,4)*(*mixC)(i2,4))
	  +(*mixH)(i1,2)/sqrt(2.)*mu*(yl[0]*(*mixC)(i3,2)*(*mixC)(i2,2)+yl[1]*(*mixC)(i3,3)*(*mixC)(i2,3)+yl[2]*(*mixC)(i3,4)*(*mixC)(i2,4))
	  // HL
	  -0.25*sqr(g)*(((*mixH)(i1,2)*(*mixC)(i3,2)+(*mixH)(i1,3)*(*mixC)(i3,3)+(*mixH)(i1,4)*(*mixC)(i3,4))*
			(vevs[0]*(*mixC)(i2,0)+vevs[1]*(*mixC)(i2,1))
			+(vevs[2]*(*mixC)(i3,2)+vevs[3]*(*mixC)(i3,3)+vevs[4]*(*mixC)(i3,4))*
			((*mixH)(i1,0)*(*mixC)(i2,0)+(*mixH)(i1,1)*(*mixC)(i2,1)))
	  +0.5*(*mixH)(i1,0)*(*mixC)(i2,0)*
	  (sqr(yl[0])*vevs[2]*(*mixC)(i3,2)+sqr(yl[1])*vevs[3]*(*mixC)(i3,3)+sqr(yl[2])*vevs[4]*(*mixC)(i3,4))
	  +0.5*vevs[0]*(*mixC)(i2,0)*(sqr(yl[0])*(*mixH)(i1,2)*(*mixC)(i3,2)+
				      sqr(yl[1])*(*mixH)(i1,3)*(*mixC)(i3,3)+
				      sqr(yl[2])*(*mixH)(i1,4)*(*mixC)(i3,4))
	  // LH
	  -0.25*sqr(g)*(((*mixH)(i1,2)*(*mixC)(i2,2)+(*mixH)(i1,3)*(*mixC)(i2,3)+(*mixH)(i1,4)*(*mixC)(i2,4))*
			(vevs[0]*(*mixC)(i3,0)+vevs[1]*(*mixC)(i3,1))
			+(vevs[2]*(*mixC)(i2,2)+vevs[3]*(*mixC)(i2,3)+vevs[4]*(*mixC)(i2,4))*
			((*mixH)(i1,0)*(*mixC)(i3,0)+(*mixH)(i1,1)*(*mixC)(i3,1)))
	  +0.5*(*mixH)(i1,0)*(*mixC)(i3,0)*
	  (sqr(yl[0])*vevs[2]*(*mixC)(i2,2)+sqr(yl[1])*vevs[3]*(*mixC)(i2,3)+sqr(yl[2])*vevs[4]*(*mixC)(i2,4))
	  +0.5*vevs[0]*(*mixC)(i3,0)*(sqr(yl[0])*(*mixH)(i1,2)*(*mixC)(i2,2)+
				      sqr(yl[1])*(*mixH)(i1,3)*(*mixC)(i2,3)+
				      sqr(yl[2])*(*mixH)(i1,4)*(*mixC)(i2,4))
	  // HR
	  +sqrt(0.5)*(eps[0]*yl[0]*(*mixC)(i3,5)+eps[1]*yl[1]*(*mixC)(i3,6)+eps[2]*yl[2]*(*mixC)(i3,7))*
	  ((*mixH)(i1,0)*(*mixC)(i2,1)+(*mixH)(i1,1)*(*mixC)(i2,0))
	  +sqrt(0.5)*(*mixC)(i2,0)*(Al[0]*(*mixH)(i1,2)*(*mixC)(i3,5)+
				    Al[1]*(*mixH)(i1,2)*(*mixC)(i3,6)+Al[2]*(*mixH)(i1,2)*(*mixC)(i3,7))
	  +sqrt(0.5)*mu*(*mixC)(i2,1)*(yl[0]*(*mixH)(i1,2)*(*mixC)(i3,5)+
				       yl[1]*(*mixH)(i1,3)*(*mixC)(i3,6)+
				       yl[2]*(*mixH)(i1,4)*(*mixC)(i3,7))
	  // RH
	  +sqrt(0.5)*(eps[0]*yl[0]*(*mixC)(i2,5)+eps[1]*yl[1]*(*mixC)(i2,6)+eps[2]*yl[2]*(*mixC)(i2,7))*
	  ((*mixH)(i1,0)*(*mixC)(i3,1)+(*mixH)(i1,1)*(*mixC)(i3,0))
	  +sqrt(0.5)*(*mixC)(i3,0)*(Al[0]*(*mixH)(i1,2)*(*mixC)(i2,5)+
				    Al[1]*(*mixH)(i1,2)*(*mixC)(i2,6)+Al[2]*(*mixH)(i1,2)*(*mixC)(i2,7))
	  +sqrt(0.5)*mu*(*mixC)(i3,1)*(yl[0]*(*mixH)(i1,2)*(*mixC)(i2,5)+
				       yl[1]*(*mixH)(i1,3)*(*mixC)(i2,6)+
				       yl[2]*(*mixH)(i1,4)*(*mixC)(i2,7));
	// normalization
	scalarChargedCharged_[i1][i2][i3] /= g;
      }
    }
  }
  // coupling of the pseudoscalar Higgs to two charged Higgs bosons
  pseudoChargedCharged_.resize(pseudo.size(),vector<vector<complex<Energy > > >
  			       (charged.size(),vector<complex<Energy> >(charged.size(),complex<Energy>(ZERO))));
  if(charged.size()<5) {
    for(unsigned int i1=0;i1<pseudo.size();++i1) {
      for(unsigned int i2=0;i2<charged.size();++i2) {
	for(unsigned int i3=0;i3<charged.size();++i3) {
	  pseudoChargedCharged_[i1][i2][i3] = ZERO;
	  for(unsigned int il=0;il<3;++il) {
	    pseudoChargedCharged_[i1][i2][i3] += (Al[il]*(*mixP)(i1,0)-mu*yl[il]*(*mixP)(i1,1))*
	      ((*mixC)(i2,il+2)*(*mixC)(i3,il+5)-(*mixC)(i2,il+5)*(*mixC)(i3,il+2));
	  }
	  // normalization // do not revert to *= , breaks with XCode 5.1
	  pseudoChargedCharged_[i1][i2][i3] = pseudoChargedCharged_[i1][i2][i3] * Complex(0.,-1.)*sqrt(0.5)/g;
	}
      }
    }
  }
  // couplings of the Higgs bosons to the up type squarks
  scalarSup_.resize(scalar.size(),vector<vector<vector<complex<Energy > > > >
		    (3,vector<vector<complex<Energy> > >(2,vector<complex<Energy> >(2,complex<Energy>(ZERO)))));
  for(unsigned int iq=0;iq<3;++iq) {
    double y = sqrt(2.)*getParticleData(2*long(iq+1))->mass()/vevs[1];
    complex<Energy> A = iq!=2 ? complex<Energy>(ZERO) : y*model->topTrilinear();
    Complex mixing[2][2];
    if(iq!=2) {
      mixing[0][0] = 1.;
      mixing[0][1] = 0.;
      mixing[1][0] = 0.;
      mixing[1][1] = 1.;
    }
    else {
      mixing[0][0] = (*stop_)(0,0);
      mixing[0][1] = (*stop_)(0,1);
      mixing[1][0] = (*stop_)(1,0);
      mixing[1][1] = (*stop_)(1,1);
    }
    // loop over Higgs
    for(unsigned int ix=0;ix<2;++ix) {
      for(unsigned int iy=0;iy<2;++iy) {
	for(unsigned int ih=0;ih<scalar.size();++ih) {
	  // first LL term
	  Complex LL1 =-(sqr(g)/4.- sqr(gp)/12.)*mixing[ix][0]*conj(mixing[iy][0]);
	  // first RR term 
	  Complex RR1 = -1./3.*sqr(gp)*mixing[ix][1]*conj(mixing[iy][1]);
	  complex<Energy> factors[5];
	  for(unsigned int iz=0;iz<5;++iz) factors[iz] = delta[iz]*vevs[iz]*(LL1+RR1);
	  factors[0] += y/sqrt(2.)*mu*( mixing[ix][0]*conj(mixing[iy][1]) + 
					mixing[ix][1]*conj(mixing[iy][0]) );
	  factors[1] += -vevs[1]*sqr(y)*( mixing[ix][0]*conj(mixing[iy][0]) + 
					  mixing[ix][1]*conj(mixing[iy][1]) )
	    -A/sqrt(2.)*(mixing[ix][0]*conj(mixing[iy][1]) + 
			 mixing[ix][1]*conj(mixing[iy][0]));
	  for(unsigned int iz=0;iz<3;++iz) {
	    factors[iz+2] += -y*eps[iz]/sqrt(2.)*( mixing[ix][0]*conj(mixing[iy][1]) + 
						    mixing[ix][1]*conj(mixing[iy][0]) );
	  }
	  for(unsigned int iz=0;iz<scalar.size();++iz) {
	    scalarSup_[ih][iq][ix][iy] += (*mixH)(ih,iz)*factors[iz];
	  }
	  scalarSup_[ih][iq][ix][iy] /= g;
	}
      }
    }
  }
  // couplings of the pseudoscalar Higgs bosons to the up type squarks
  pseudoSup_.resize(pseudo.size(),vector<complex<Energy > >(3,complex<Energy>(ZERO)));
  for(unsigned int iq=0;iq<3;++iq) {
    double y = sqrt(2.)*getParticleData(2*long(iq+1))->mass()/vevs[1];
    complex<Energy> A = iq!=2 ? complex<Energy>(ZERO) : y*model->topTrilinear();
    // loop over Higgs
    for(unsigned int ih=0;ih<pseudo.size();++ih) {
      pseudoSup_[ih][iq] = Complex(0.,-1.)*sqrt(0.5)*(A*(*mixP)(ih,1)-y*mu*(*mixP)(ih,0))/g;
    }
  }
  // couplings of the Higgs bosons to the down type squarks
  scalarSdown_.resize(scalar.size(),vector<vector<vector<complex<Energy > > > >
		      (3,vector<vector<complex<Energy> > >(2,vector<complex<Energy> >(2,complex<Energy>(ZERO)))));
  for(unsigned int iq=0;iq<3;++iq) {
    double y = sqrt(2.)*getParticleData(2*long(iq+1)-1)->mass()/vevs[0];
    Complex mixing[2][2];
    complex<Energy> A = iq!=2 ? complex<Energy>(ZERO) : y*model->bottomTrilinear();
    if(iq!=2) {
      mixing[0][0] = 1.;
      mixing[0][1] = 0.;
      mixing[1][0] = 0.;
      mixing[1][1] = 1.;
    }
    else {
      mixing[0][0] = (*sbottom_)(0,0);
      mixing[0][1] = (*sbottom_)(0,1);
      mixing[1][0] = (*sbottom_)(1,0);
      mixing[1][1] = (*sbottom_)(1,1);
    }
    for(unsigned int ix=0;ix<2;++ix) {
      for(unsigned int iy=0;iy<2;++iy) {
	for(unsigned int ih=0;ih<scalar.size();++ih) {
	  // first LL term
	  Complex LL1 = (sqr(g)/4.+sqr(gp)/12.)*mixing[ix][0]*conj(mixing[iy][0]);
	  // first RR term 
	  Complex RR1 = 1./6.*sqr(gp)*mixing[ix][1]*conj(mixing[iy][1]);
	  complex<Energy> factors[5];
	  for(unsigned int iz=0;iz<5;++iz) factors[iz] = delta[iz]*vevs[iz]*(LL1+RR1);
	  factors[0] += -vevs[0]*sqr(y)*(mixing[ix][0]*conj(mixing[iy][0])+
					 mixing[ix][1]*conj(mixing[iy][1]))
	    - A/sqrt(2.)*(mixing[ix][0]*conj(mixing[iy][1])+
			  mixing[ix][1]*conj(mixing[iy][0]));
	  factors[1] +=  y/sqrt(2.)*mu*(mixing[ix][0]*conj(mixing[iy][1])+
					mixing[ix][1]*conj(mixing[iy][0]));
	  for(unsigned int iz=0;iz<scalar.size();++iz) {
	    scalarSdown_[ih][iq][ix][iy] += (*mixH)(ih,iz)*factors[iz];
	  }
	  scalarSdown_[ih][iq][ix][iy] /= g;
	}
      }
    }
  }
  // couplings of the pseudoscalar Higgs bosons to the down type squarks
  pseudoSdown_.resize(pseudo.size(),vector<complex<Energy > >(3,complex<Energy>(ZERO)));
  for(unsigned int iq=0;iq<3;++iq) {
    double y = sqrt(2.)*getParticleData(2*long(iq+1)-1)->mass()/vevs[0];
    complex<Energy> A = iq!=2 ? complex<Energy>(ZERO) : y*model->bottomTrilinear();
    for(unsigned int ih=0;ih<pseudo.size();++ih) {
      pseudoSdown_[ih][iq] = Complex(0.,-1.)*sqrt(0.5)*(A*(*mixP)(ih,0)-mu*y*(*mixP)(ih,1))/g;
    }
  }
  // couplings of the scalar Higgs bosons to the sneutrinos
  if(scalar.size()==2) {
    scalarSneutrino_.resize(scalar.size(),vector<complex<Energy > >(3,complex<Energy>(ZERO)));
    for(unsigned int il=0;il<3;++il) {
      // loop over Higgs
      for(unsigned int ih=0;ih<scalar.size();++ih) {
	// first LL term
	Complex LL1 =-(sqr(g)/4.+sqr(gp)/4.);
	complex<Energy> factors[5];
	for(unsigned int iz=0;iz<5;++iz) factors[iz] = delta[iz]*vevs[iz]*LL1;
	for(unsigned int iz=0;iz<scalar.size();++iz) {
	  scalarSneutrino_[ih][il] += (*mixH)(ih,iz)*factors[iz];
	}
	scalarSneutrino_[ih][il] /= g;
      }
    }
  }
  // couplings of the Higgs bosons to the charged sleptons
  if(charged.size()==1) {
    scalarSlepton_.resize(scalar.size(),vector<vector<vector<complex<Energy > > > >
			  (3,vector<vector<complex<Energy> > >(2,vector<complex<Energy> >(2,complex<Energy>(ZERO)))));
    for(unsigned int il=0;il<3;++il) {
      double y = sqrt(2.)*getParticleData(2*long(il+6)-1)->mass()/vevs[0];
      Complex mixing[2][2];
      complex<Energy> A = il!=2 ? complex<Energy>(ZERO) : y*model->tauTrilinear();
      if(il!=2) {
	mixing[0][0] = 1.;
	mixing[0][1] = 0.;
	mixing[1][0] = 0.;
	mixing[1][1] = 1.;
      }
      else {
	mixing[0][0] = (*stau_)(0,0);
	mixing[0][1] = (*stau_)(0,1);
	mixing[1][0] = (*stau_)(1,0);
	mixing[1][1] = (*stau_)(1,1);
      }
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  for(unsigned int ih=0;ih<scalar.size();++ih) {
	    // first LL term
	    Complex LL1 = (sqr(g)/4.-sqr(gp)/4.)*mixing[ix][0]*conj(mixing[iy][0]);
	    // first RR term 
	    Complex RR1 = 0.5*sqr(gp)*mixing[ix][1]*conj(mixing[iy][1]);
	    complex<Energy> factors[5];
	    for(unsigned int iz=0;iz<5;++iz) factors[iz] = delta[iz]*vevs[iz]*(LL1+RR1);
	    factors[0] += -vevs[0]*sqr(y)*(mixing[ix][0]*conj(mixing[iy][0])+
					   mixing[ix][1]*conj(mixing[iy][1]))
	      - A/sqrt(2.)*(mixing[ix][0]*conj(mixing[iy][1])+
			    mixing[ix][1]*conj(mixing[iy][0]));
	    factors[1] +=  y/sqrt(2.)*mu*(mixing[ix][0]*conj(mixing[iy][1])+
					  mixing[ix][1]*conj(mixing[iy][0]));
	    for(unsigned int iz=0;iz<scalar.size();++iz) {
	      scalarSlepton_[ih][il][ix][iy] += (*mixH)(ih,iz)*factors[iz];
	    }
	    scalarSlepton_[ih][il][ix][iy] /= g;
	  }
	}
      }
    }
  }
  // couplings of the pseudoscalar Higgs bosons to the charged sleptons
  if(charged.size()==1) {
    pseudoSlepton_.resize(pseudo.size(),vector<complex<Energy> >(3,complex<Energy>(ZERO)));
    for(unsigned int il=0;il<3;++il) {
      double y = sqrt(2.)*getParticleData(2*long(il+6)-1)->mass()/vevs[0];
      complex<Energy> A = il!=2 ? complex<Energy>(ZERO) : y*model->tauTrilinear();
      for(unsigned int ih=0;ih<pseudo.size();++ih) {
	pseudoSlepton_[ih][il] = Complex(0.,-1.)*sqrt(0.5)*(A*(*mixP)(ih,0)-mu*y*(*mixP)(ih,1))/g;
      }
    }
  }
  // charged Higgs squarks
  chargedSquark_.resize(charged.size(),vector<vector<vector<complex<Energy > > > >
			(3,vector<vector<complex<Energy> > >(2,vector<complex<Energy> >(2,complex<Energy>(ZERO)))));
  for(unsigned int iq=0;iq<3;++iq) {
    double yd = sqrt(2.)*getParticleData(2*long(iq+1)-1)->mass()/vevs[0];
    double yu = sqrt(2.)*getParticleData(2*long(iq+1))->mass()/vevs[1];
    Complex mixd[2][2],mixu[2][2];
    complex<Energy> Ad = iq!=2 ? complex<Energy>(ZERO) : yd*model->bottomTrilinear();
    complex<Energy> Au = iq!=2 ? complex<Energy>(ZERO) : yu*model->topTrilinear();
    if(iq!=2) {
      mixd[0][0] = mixu[0][0] = 1.;
      mixd[0][1] = mixu[0][1] = 0.;
      mixd[1][0] = mixu[1][0] = 0.;
      mixd[1][1] = mixu[1][1] = 1.;
    }
    else {
      mixd[0][0] = (*sbottom_)(0,0);
      mixd[0][1] = (*sbottom_)(0,1);
      mixd[1][0] = (*sbottom_)(1,0);
      mixd[1][1] = (*sbottom_)(1,1);
      mixu[0][0] = (*stop_)(0,0);
      mixu[0][1] = (*stop_)(0,1);
      mixu[1][0] = (*stop_)(1,0);
      mixu[1][1] = (*stop_)(1,1);
    }
    for(unsigned int ih=0;ih<charged.size();++ih) {
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  // various terms (up first) (down second)
	  complex<Energy> LL(ZERO);
	  for(unsigned int iz=0;iz<vevs.size();++iz)
	    LL += 0.5*sqr(g)*vevs[iz]*(*mixC)(ih,iz);
	  LL += -vevs[0]*(*mixC)(ih,0)*sqr(yd)-vevs[1]*(*mixC)(ih,1)*sqr(yu);
	  LL *= -sqrt(0.5);
	  complex<Energy> RR = 
	    yu*yd*(vevs[1]*(*mixC)(ih,0) + vevs[0]*(*mixC)(ih,1))/sqrt(2.);
	  complex<Energy> LR = (*mixC)(ih,0)*Ad+(*mixC)(ih,1)*yd*mu;
	  complex<Energy> RL = yu*mu*(*mixC)(ih,0) + Au*(*mixC)(ih,1);
	  chargedSquark_[ih][iq][ix][iy] = 
	    mixu[ix][0]*mixd[iy][0]*LL + mixu[ix][1]*mixd[iy][1]*RR +
	    mixu[ix][0]*mixd[iy][1]*LR + mixu[ix][1]*mixd[iy][0]*RL;
	  chargedSquark_[ih][iq][ix][iy] /=g;
	}
      }
    }
  }
  // charged Higgs slepton
  if(charged.size()==1) {
    chargedSlepton_.resize(charged.size(),vector<vector<complex<Energy > > >
			   (3,vector<complex<Energy> >(2,complex<Energy>(ZERO))));
    for(unsigned int il=0;il<3;++il) {
      double y = sqrt(2.)*getParticleData(2*long(il+6)-1)->mass()/vevs[0];
      Complex mixd[2][2];
      complex<Energy> Al = il!=2 ? complex<Energy>(ZERO) : y*model->tauTrilinear();
      if(il!=2) {
	mixd[0][0] = 1.;
	mixd[0][1] = 0.;
	mixd[1][0] = 0.;
	mixd[1][1] = 1.;
      }
      else {
	mixd[0][0] = (*stau_)(0,0);
	mixd[0][1] = (*stau_)(0,1);
	mixd[1][0] = (*stau_)(1,0);
	mixd[1][1] = (*stau_)(1,1);
      }
      for(unsigned int ih=0;ih<charged.size();++ih) {
	for(unsigned int ix=0;ix<2;++ix) {
	  // various terms (charged lepton second)
	  complex<Energy> LL(ZERO);
	  for(unsigned int iz=0;iz<vevs.size();++iz)
	    LL += 0.5*sqr(g)*vevs[iz]*(*mixC)(ih,iz);
	  LL += -vevs[0]*(*mixC)(ih,0)*sqr(y);
	  LL *= -sqrt(0.5);
	  complex<Energy> LR = (*mixC)(ih,0)*Al+(*mixC)(ih,1)*y*mu;
	  chargedSlepton_[ih][il][ix] = mixd[ix][0]*LL + mixd[ix][1]*LR;
	  chargedSlepton_[ih][il][ix] /=g;
	}
      }
    }
  }
}

namespace {

int scalarEigenState(long id) {
  if(id<1000000)
    return (id-25)/10;
  else
    return (id-1000008)/2;
}

int pseudoEigenState(long id) {
  if(id<1000000)
    return 0;
  else
    return (id-1000016);
}

int chargedEigenState(long id) {
  if(abs(id)<1000000)
    return 0;
  else if(abs(id)<2000000)
    return (abs(id)-1000009)/2;
  else
    return (abs(id)-2000003)/2;
}

}

void RPVSSSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  // prefactor
  if( q2 != q2Last_ || gLast_==0. ) {
    q2Last_ = q2;
    gLast_ = weakCoupling(q2);
  }
  // some ordering of the particles
  // sort in order of PDG codes, smallest first
  if(abs(part1->id())>abs(part2->id())) swap(part1,part2);
  if(abs(part1->id())>abs(part3->id())) swap(part1,part3);
  if(abs(part2->id())>abs(part3->id())) swap(part2,part3);
  // make sure squarks 2nd and 3rd
  if(abs(part1->id())%1000000<=6) swap(part1,part2);
  if(abs(part1->id())%1000000<=6) swap(part1,part3);
  // extract particle ids
  long sca1(part1->id()), sca2(part2->id()), sca3(part3->id());
  // Higgs squark couplings
  if(abs(sca2)%1000000<=6) {
    // charged Higgs
    if( part1->charged() ) {
      int ih = chargedEigenState(sca1);
      if(abs(sca2)%2!=0) swap(sca2,sca3);
      unsigned int alpha = abs(sca2)/1000000-1;
      unsigned int beta  = abs(sca3)/1000000-1;
      unsigned int iq = (abs(sca2)%1000000-2)/2;
      norm(gLast_*chargedSquark_[ih][iq][alpha][beta]*UnitRemoval::InvE);
      return;
    }
    // neutral Higgs
    else {
      long sm(0);
      // get the left/right light/heavy state
      unsigned int alpha(abs(sca2)/1000000 - 1), beta(abs(sca3)/1000000 - 1);
      sm = abs(sca2)%1000000;
      bool pseudo = sca1==36 || (sca1>=1000017&&sca1<=1000019);
      int higgs = pseudo ? pseudoEigenState(sca1) : scalarEigenState(sca1);
      complex<Energy> coup;
      if(!pseudo) {
	if( sm % 2 == 0 ) {
	  coup = scalarSup_  [higgs][(sm-2)/2][alpha][beta];
	}
	else {
	  coup = scalarSdown_[higgs][(sm-1)/2][alpha][beta];
	}
      }
      else {
	if( sm % 2 == 0 ) {
	  coup = pseudoSup_  [higgs][(sm-2)/2];
	}
	else {
	  coup = pseudoSdown_[higgs][(sm-1)/2];
	}
	if((alpha==1&&sca2<0)||(beta==1&&sca3<0)) coup *=-1.;
      }
      // set coupling and return
      norm(gLast_*coup*UnitRemoval::InvE);
      return;
    }
  }
  // all neutral
  else if(!part1->charged()&&!part2->charged()&&!part3->charged()) {
    // neutral Higgs sneutrino
    if(scalarScalarScalar_.size()==2&&abs(sca2)>1000000) {
      assert(!(sca1==36 || (sca1>=1000017&&sca1<=1000019)));
      int il = (abs(sca2)-1000012)/2;
      norm(gLast_*scalarSneutrino_[scalarEigenState(sca1)][il]*UnitRemoval::InvE);
      return;
    }
    else {
      if(sca1==36 || (sca1>=1000017&&sca2<=1000019)) swap(sca1,sca2);
      if(sca1==36 || (sca1>=1000017&&sca2<=1000019)) swap(sca1,sca3);
      // 2 pseudoscalar 1 scale
      if(sca2==36 || (sca2>=1000017&&sca2<=1000019)) {
	norm(gLast_*scalarPseudoPseudo_[scalarEigenState(sca1)]
	     [pseudoEigenState(sca2)][pseudoEigenState(sca3)]*UnitRemoval::InvE);
	return;
      }
      // 3 scalars
      else {
	norm(gLast_*scalarScalarScalar_[scalarEigenState(sca1)]
	     [scalarEigenState(sca2)][scalarEigenState(sca3)]*UnitRemoval::InvE);
	return;
      }
    }
  }
  // two charged
  else {
    // put the charged last
    if(!part2->charged()) {
      swap(part1,part2);
      swap(sca1 ,sca2 );
    }
    if(!part3->charged()) {
      swap(part1,part3);
      swap(sca1 ,sca3 );
    }
    // sleptons
    if(scalarChargedCharged_[0].size()<5&&abs(sca2)>1000000) {
      // neutral Higgs charged sleptons
      if((abs(sca2)>=1000011&&abs(sca2)<=1000015&&abs(sca2)%2!=0)&&
	 (abs(sca3)>=1000011&&abs(sca3)<=1000015&&abs(sca3)%2!=0)) {
	long sm(0);
	// get the left/right light/heavy state
	unsigned int alpha(abs(sca2)/1000000 - 1), beta(abs(sca3)/1000000 - 1);
	sm = (abs(sca2)%1000000-11)/2;
	bool pseudo = sca1==36 || (sca1>=1000017&&sca1<=1000019);
	int higgs = pseudo ? pseudoEigenState(sca1) : scalarEigenState(sca1);
	complex<Energy> coup;
	if(!pseudo) {
	  coup = scalarSlepton_[higgs][(sm-1)/2][alpha][beta];
	}
	else {
    	  coup = pseudoSlepton_[higgs][(sm-1)/2];
	  if((alpha==1&&sca2<0)||(beta==1&&sca3<0)) coup *=-1.;
    	}
	// set coupling and return
	norm(gLast_*coup*UnitRemoval::InvE);
	return;
      }
      // charged Higgs
      else {
	if(abs(sca2)<1000000) {
	  swap(part1,part2);
	  swap(sca1 ,sca2 );
	}
	if(abs(sca3)<1000000) {
	  swap(part1,part3);
	  swap(sca1 ,sca3 );
	}
	if(part3->charged()) {
	  swap(part2,part3);
	  swap(sca2 ,sca3 );
	}
	unsigned int il    = (abs(sca2)%1000000-11)/2;
	unsigned int alpha = abs(sca2)/1000000-1;
	norm(gLast_*chargedSlepton_[chargedEigenState(sca1)][il][alpha]*UnitRemoval::InvE);
      }
    }
    else {
      // pseudoscalar charged charged
      if(sca1==36 || (sca1>=1000017&&sca1<=1000019)) {
	if(sca2>0) swap(sca2,sca3);
	norm(gLast_*pseudoChargedCharged_[pseudoEigenState(sca1)]
	     [chargedEigenState(sca2)][chargedEigenState(sca3)]*UnitRemoval::InvE);
      }
      // scalar charged charged
      else {
	norm(gLast_*scalarChargedCharged_[scalarEigenState(sca1)]
	     [chargedEigenState(sca2)][chargedEigenState(sca3)]*UnitRemoval::InvE);
      }
    }
    return;
  }	
}
