#line 1 "./SusyBase.cc"
// -*- C++ -*-
//
// SusyBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SusyBase class.
//

#include "SusyBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;

SusyBase::SusyBase() : readFile_(false), MPlanck_(2.4e18*GeV),
		       gravitino_(false), majoranaNeutrinos_(false),
		       tanBeta_(0), mu_(ZERO), 
		       M1_(ZERO), M2_(ZERO), M3_(ZERO),
		       mH12_(ZERO),mH22_(ZERO),
		       meL_(ZERO),mmuL_(ZERO),mtauL_(ZERO),
		       meR_(ZERO),mmuR_(ZERO),mtauR_(ZERO),
		       mq1L_(ZERO),mq2L_(ZERO),mq3L_(ZERO),
		       mdR_(ZERO),muR_(ZERO),msR_(ZERO),
		       mcR_(ZERO),mbR_(ZERO),mtR_(ZERO),
		       gluinoPhase_(1.), allowedToResetSMMasses_(false)
{}

IBPtr SusyBase::clone() const {
  return new_ptr(*this);
}

IBPtr SusyBase::fullclone() const {
  return new_ptr(*this);
}

void SusyBase::doinit() {
  addVertex(WSFSFVertex_);
  addVertex(WWSFSFVertex_);
  addVertex(NFSFVertex_);
  addVertex(GFSFVertex_);
  addVertex(HSFSFVertex_);
  addVertex(CFSFVertex_);
  addVertex(GSFSFVertex_);
  addVertex(GGSQSQVertex_);
  addVertex(WGSQSQVertex_);
  addVertex(GSGSGVertex_);
  addVertex(NNZVertex_);
  if(NNPVertex_) addVertex(NNPVertex_);
  if(GNGVertex_) addVertex(GNGVertex_);
  addVertex(CCZVertex_);
  addVertex(CNWVertex_);
  addVertex(GOGOHVertex_);
  addVertex(WHHVertex_);
  if(NCTVertex_) addVertex(NCTVertex_);
  if(gravitino_) {
    if(GVNHVertex_) addVertex(GVNHVertex_);
    if(GVNVVertex_) addVertex(GVNVVertex_);
    if(GVFSVertex_) addVertex(GVFSVertex_);
  }
  if(majoranaNeutrinos()) {
    idMap().insert(make_pair(12,17));
    idMap().insert(make_pair(14,18));
    idMap().insert(make_pair(16,19));
  }
  BSMModel::doinit();
}

void SusyBase::persistentOutput(PersistentOStream & os) const {
  os << readFile_ << gravitino_
     << NMix_ << UMix_ << VMix_ << WSFSFVertex_ << WWSFSFVertex_ 
     << NFSFVertex_ << GFSFVertex_ << HSFSFVertex_ << CFSFVertex_ 
     << GSFSFVertex_ << GGSQSQVertex_ << WGSQSQVertex_ << GSGSGVertex_ 
     << NNZVertex_ << NNPVertex_ << CCZVertex_ << CNWVertex_ 
     << GOGOHVertex_ << WHHVertex_ << GNGVertex_ << NCTVertex_
     << GVNHVertex_ << GVNVVertex_ << GVFSVertex_
     << tanBeta_ << ounit(mu_,GeV) 
     << ounit(M1_,GeV) << ounit(M2_,GeV) << ounit(M3_,GeV)
     << ounit(mH12_,GeV2) << ounit(mH22_,GeV2) 
     << ounit(meL_,GeV)  << ounit(mmuL_,GeV) << ounit(mtauL_,GeV) 
     << ounit(meR_,GeV)  << ounit(mmuR_,GeV) << ounit(mtauR_,GeV) 
     << ounit(mq1L_,GeV) << ounit(mq2L_,GeV) << ounit(mq3L_,GeV) 
     << ounit(mdR_,GeV)  << ounit(muR_,GeV)  << ounit(msR_,GeV) 
     << ounit(mcR_,GeV)  << ounit(mbR_,GeV)  << ounit(mtR_,GeV)
     << gluinoPhase_ << ounit(MPlanck_,GeV) << allowedToResetSMMasses_;
}

void SusyBase::persistentInput(PersistentIStream & is, int) {
  is >> readFile_  >> gravitino_
     >> NMix_ >> UMix_ >> VMix_ >> WSFSFVertex_ >> WWSFSFVertex_ 
     >> NFSFVertex_ >> GFSFVertex_ >> HSFSFVertex_ >> CFSFVertex_ 
     >> GSFSFVertex_ >> GGSQSQVertex_ >> WGSQSQVertex_ >> GSGSGVertex_ 
     >> NNZVertex_ >> NNPVertex_ >> CCZVertex_ >> CNWVertex_
     >> GOGOHVertex_ >> WHHVertex_ >> GNGVertex_ >> NCTVertex_
     >> GVNHVertex_ >> GVNVVertex_ >> GVFSVertex_
     >> tanBeta_ >> iunit(mu_,GeV) 
     >> iunit(M1_,GeV) >> iunit(M2_,GeV) >> iunit(M3_,GeV)
     >> iunit(mH12_,GeV2) >> iunit(mH22_,GeV2) 
     >> iunit(meL_,GeV)  >> iunit(mmuL_,GeV) >> iunit(mtauL_,GeV) 
     >> iunit(meR_,GeV)  >> iunit(mmuR_,GeV) >> iunit(mtauR_,GeV) 
     >> iunit(mq1L_,GeV) >> iunit(mq2L_,GeV) >> iunit(mq3L_,GeV) 
     >> iunit(mdR_,GeV)  >> iunit(muR_,GeV)  >> iunit(msR_,GeV) 
     >> iunit(mcR_,GeV)  >> iunit(mbR_,GeV)  >> iunit(mtR_,GeV)
     >> gluinoPhase_ >> iunit(MPlanck_,GeV) >> allowedToResetSMMasses_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SusyBase,BSMModel>
describeHerwigSusyBase("Herwig::SusyBase", "HwSusy.so");

void SusyBase::Init() {

  static ClassDocumentation<SusyBase> documentation
    ("This is the base class for any SUSY model.",
     "SUSY spectrum files follow the Les Houches accord"
     " \\cite{Skands:2003cj,Allanach:2008qq}.",
     " %\\cite{Skands:2003cj}\n"
     "\\bibitem{Skands:2003cj}\n"
     "  P.~Skands {\\it et al.},\n"
     "   ``SUSY Les Houches accord: Interfacing SUSY spectrum calculators, decay\n"
     "  %packages, and event generators,''\n"
     "  JHEP {\\bf 0407}, 036 (2004)\n"
     "  [arXiv:hep-ph/0311123].\n"
     "  %%CITATION = JHEPA,0407,036;%%\n"
     "%\\cite{Allanach:2008qq}\n"
     "\\bibitem{Allanach:2008qq}\n"
     "  B.~Allanach {\\it et al.},\n"
     "  %``SUSY Les Houches Accord 2,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 180}, 8 (2009)\n"
     "  [arXiv:0801.0045 [hep-ph]].\n"
     "  %%CITATION = CPHCB,180,8;%%\n"
     );


  static Reference<SusyBase,Helicity::AbstractVSSVertex> interfaceVertexWSS
    ("Vertex/WSFSF",
     "Reference to Susy W SF SF vertex",
     &SusyBase::WSFSFVertex_, false, false, true, false);
  
  static Reference<SusyBase,Helicity::AbstractVVSSVertex> interfaceVertexWWSS
    ("Vertex/WWSFSF",
     "Reference to Susy W W SF SF vertex",
     &SusyBase::WWSFSFVertex_, false, false, true, false);
  
  static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexNFSF
    ("Vertex/NFSF",
     "Reference to the neutralino-fermion-sfermion vertex",
     &SusyBase::NFSFVertex_, false, false, true, false);

  static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexGFSF
    ("Vertex/GFSF",
     "Reference to the gluino-fermion-sfermion vertex",
     &SusyBase::GFSFVertex_, false, false, true, false);
  
  static Reference<SusyBase,Helicity::AbstractSSSVertex> interfaceVertexHSFSF
    ("Vertex/HSFSF",
     "Reference to the Higgs-fermion-sfermion vertex",
     &SusyBase::HSFSFVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexCFSF
   ("Vertex/CFSF",
      "Reference to the chargino-fermion-sfermion vertex",
      &SusyBase::CFSFVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractVSSVertex> interfaceVertexGSFSF
   ("Vertex/GSFSF",
      "Reference to the gluon-sfermion-sfermion vertex",
      &SusyBase::GSFSFVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractVVSSVertex> interfaceVertexGGSS
     ("Vertex/GGSQSQ",
      "Reference to the gluon-gluon-squark-squark vertex",
      &SusyBase::GGSQSQVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractVVSSVertex> interfaceVertexWGSS
     ("Vertex/WGSQSQ",
      "Reference to the gauge boson-gluon-squark-squark vertex",
      &SusyBase::WGSQSQVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexGSGSG
     ("Vertex/GSGSG",
      "Reference to the gluon-gluino-gluino vertex",
      &SusyBase::GSGSGVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexNNZ
    ("Vertex/NNZ",
     "Reference to Z-~chi_i0-~chi_i0 vertex",
     &SusyBase::NNZVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexNNP
    ("Vertex/NNP",
     "Reference to photon-~chi_i0-~chi_i0 vertex",
     &SusyBase::NNPVertex_, false, false, true, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexGNG
    ("Vertex/GNG",
     "Reference to gluon-~chi_i0-gluino vertex",
     &SusyBase::GNGVertex_, false, false, true, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexCCZ
    ("Vertex/CCZ",
     "Reference to ~chi_i+-~chi_i-Z vertex",
     &SusyBase::CCZVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexCNW
    ("Vertex/CNW",
     "Reference to ~chi_i+-chi_i0-W vertex",
     &SusyBase::CNWVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexGOGOH
   ("Vertex/GOGOH",
    "Reference to the gaugino-gaugino-higgs vertex",
    &SusyBase::GOGOHVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractVSSVertex> interfaceVertexWHH
    ("Vertex/SSWHH",
     "Reference to Susy WHHVertex",
     &SusyBase::WHHVertex_, false, false, true, false);

  static Reference<SusyBase,AbstractFFSVertex> interfaceVertexNCT
    ("Vertex/NCT",
     "Vertex for the flavour violating coupling of the top squark "
     "to the neutralino and charm quark.",
     &SusyBase::NCTVertex_, false, false, true, true, false);

  static Reference<SusyBase,AbstractRFSVertex> interfaceVertexGVNH
    ("Vertex/GVNH",
     "Vertex for the interfaction of the gravitino-neutralino"
     " and Higgs bosons",
     &SusyBase::GVNHVertex_, false, false, true, true, false);

  static Reference<SusyBase,AbstractRFVVertex> interfaceVertexGVNV
    ("Vertex/GVNV",
     "Vertex for the interfaction of the gravitino-neutralino"
     " and vector bosons",
     &SusyBase::GVNVVertex_, false, false, true, true, false);

  static Reference<SusyBase,AbstractRFSVertex> interfaceVertexGVFS
    ("Vertex/GVFS",
     "Vertex for the interfaction of the gravitino-fermion"
     " and sfermion",
     &SusyBase::GVFSVertex_, false, false, true, true, false);

  static Parameter<SusyBase,Energy> interfaceMPlanck
    ("MPlanck",
     "The Planck mass for GMSB models",
     &SusyBase::MPlanck_, GeV, 2.4e18*GeV, 1.e16*GeV, 1.e20*GeV,
     false, false, Interface::limited);

  static Switch<SusyBase,bool> interfaceMajoranaNeutrinos
    ("MajoranaNeutrinos",
     "Whether or not the neutrinos should be treated as Majorana particles",
     &SusyBase::majoranaNeutrinos_, false, false, false);
  static SwitchOption interfaceMajoranaNeutrinosYes
    (interfaceMajoranaNeutrinos,
     "Yes",
     "Neutrinos are Majorana",
     true);
  static SwitchOption interfaceMajoranaNeutrinosNo
    (interfaceMajoranaNeutrinos,
     "No",
     "Neutrinos are Dirac fermions",
     false);

  static Switch<SusyBase,bool> interfaceAllowedToResetSMMasses
    ("AllowedToResetSMMasses",
     "Whether or not to allow SM masses to be reset via the SLHA file",
     &SusyBase::allowedToResetSMMasses_, false, false, false);
  static SwitchOption interfaceAllowedToResetSMMassesNo
    (interfaceAllowedToResetSMMasses,
     "No",
     "Not allowed",
     false);
  static SwitchOption interfaceAllowedToResetSMMassesYes
    (interfaceAllowedToResetSMMasses,
     "Yes",
     "Allowed",
     true);
}

void SusyBase::readSetup(istream & is) {
  string filename = dynamic_ptr_cast<istringstream*>(&is)->str();
  if(majoranaNeutrinos()) {
    idMap().insert(make_pair(12,17));
    idMap().insert(make_pair(14,18));
    idMap().insert(make_pair(16,19));
  }
  filename = StringUtils::stripws(filename);
  if(readFile_)
    throw SetupException() 
      << "A second SLHA file " << filename << " has been opened."
      << "This is probably unintended and as it can cause crashes"
      << " and other unpredictable behaviour it is not allowed."
      << Exception::runerror;
  CFileLineReader cfile;
  cfile.open(filename);
  if( !cfile ) throw SetupException() 
		 << "SusyBase::readSetup - An error occurred in opening the "
		 << "spectrum file \"" << filename << "\". A SUSY model cannot be "
		 << "run without this."
		 << Exception::runerror;
  useMe();
  // read first line and check if this is a Les Houches event file
  cfile.readline();
  bool lesHouches = cfile.find("<LesHouchesEvents");
  bool reading = !lesHouches;
  if(lesHouches) cfile.readline();
  //function pointer for putting all characters to lower case.
  int (*pf)(int) = tolower;
  do {
    string line = cfile.getline();
    // check for start of slha block in SLHA files
    if(lesHouches && !reading) {
      if(line.find("<slha")==0) reading = true;
      if(!cfile.readline()) break;
      continue;
    }
    // ignore comment lines
    if(line[0] == '#') {
      if(!cfile.readline()) break;
      continue;
    }
    // make everything lower case
    transform(line.begin(), line.end(), line.begin(), pf);
    // start of a block
    if(line.find("block") == 0) {
      string name = StringUtils::car(StringUtils::cdr(line), " #");
      name = StringUtils::stripws(name);
      // mixing matrix
      if((name.find("mix")  != string::npos && 
	  name.find("hmix") != 0)) {
	unsigned int row(0),col(0);
	MixingVector vals = readMatrix(cfile,row,col);
	mixings_[name] = make_pair(make_pair(row,col),vals);
      }
      else if(name.find("au") == 0 || name.find("ad") == 0 ||
	      name.find("ae") == 0 ) {
	string test = StringUtils::car(line, "#");
	while (test.find("=")!= string::npos) {
	  test = StringUtils::cdr(test, "=");
	}
	istringstream is(test);
	double scale;
	is >> scale;
	unsigned int row(0),col(0);
	MixingVector vals = readMatrix(cfile,row,col);
	if(scale>1e10) continue;
	mixings_[name] = make_pair(make_pair(row,col),vals);
      }
      else if( name == "spinfo" ) {
	readBlock(cfile,name,line,true);
      }
      else if( name.find("info") == string::npos) {
	readBlock(cfile,name,line,false);
      }
      else {
	if(!cfile.readline()) break;
      }
      continue;
    }
    else if( lesHouches && line.find("</slha") == 0 ) {
      break;
    }
    if(!cfile.readline()) break;
  }
  while(true);
  // extract the relevant parameters
  extractParameters();
  // create the mixing matrices we need
  createMixingMatrices();
  // set the masses, this has to be done after the 
  // mixing matrices have been created
  resetRepositoryMasses();
  // have now read the file
  if(decayFile()=="") decayFile(filename);
  readFile_=true;
}

void SusyBase::readBlock(CFileLineReader & cfile,string name,string linein,
			 bool stringBlock) {
  if(!cfile)
    throw SetupException() 
				    << "SusyBase::readBlock() - The input stream is in a bad state"
				    << Exception::runerror;
  // storage or the parameters
  string test = StringUtils::car(linein, "#");
  ParamMap  storeParam;
  StringMap storeString;
  bool set = true;
  // special for the alpha block
  if(name.find("alpha") == 0 ) {
    double alpha;
    cfile.readline();
    string line = cfile.getline();
    istringstream iss(line);
    iss >> alpha;
    storeParam.insert(make_pair(1,alpha));
    parameters_[name]=storeParam;
    return;
  }
  // extract the scale from the block if present
  if(test.find("=")!= string::npos) { 
    while(test.find("=")!=string::npos)
      test= StringUtils::cdr(test,"=");
    istringstream is(test);
    double scale;
    is >> scale;
    // only store the lowest scale block
    if(parameters_.find(name)!=parameters_.end()) {
      set = scale < parameters_[name][-1];
    }
    else {
      storeParam.insert(make_pair(-1,scale));
    }
  }
  while(cfile.readline()) {
    string line = cfile.getline();
    // skip comments
    if(line[0] == '#') continue;
    // reached the end
    if( line[0] == 'B' || line[0] == 'b' ||
	line[0] == 'D' || line[0] == 'd' ||
	line[0] == '<' ) {
      cfile.resetline();
      break;
    }
    istringstream is(line);
    long index;
    if(!stringBlock) {
      double value;
      if(name.find("rvlam")!= string::npos||
	 name=="rvt" || name=="rvtp" || name=="rvtpp") {
	int i,j,k;
	is >> i >> j >> k >> value;
	index = i*100+j*10+k;
      }
      else {
	is >> index >> value;
      }
      storeParam.insert(make_pair(index, value));
    }
    else {
      string value;
      is >> index >> value;
      storeString.insert(make_pair(index,value));
    }
  }
  if(set) {
    if(stringBlock) info_      [name] = storeString;
    else            parameters_[name] = storeParam ;
  }
}

const MixingVector
SusyBase::readMatrix(CFileLineReader & cfile, 
		     unsigned int & row, unsigned int & col) {
  if(!cfile)
    throw SetupException() 
      << "SusyBase::readMatrix() - The input stream is in a bad state."
      << Exception::runerror;
  unsigned int rowmax(0), colmax(0);
  MixingVector values;
  while(cfile.readline()) {
    string line = cfile.getline();
    // skip comments
    if(line[0] == '#') continue;
    // reached the end
    if( line[0] == 'B' || line[0] == 'b' ||
	line[0] == 'D' || line[0] == 'd' ||
	line[0] == '<' ) {
      cfile.resetline();
      break;
    }
    istringstream is(line);
    unsigned int index1, index2;
    double real(0.), imag(0.);   
    line = StringUtils::stripws(line);
    if(!line.empty()) {
      is >> index1 >> index2 >> real >> imag;
      values.push_back(MixingElement(index1,index2,Complex(real, imag)));
      if(index1 > rowmax) rowmax = index1;
      if(index2 > colmax) colmax = index2;
    }
  }
  col=colmax;
  row=rowmax;
  return values;
}

void SusyBase::createMixingMatrix(MixingMatrixPtr & matrix,
				  string name, const MixingVector & values,
				  MatrixSize size) {
  matrix = new_ptr(MixingMatrix(size.first,size.second));
  for(unsigned int ix=0; ix < values.size(); ++ix)
    (*matrix)(values[ix].row-1,values[ix].col-1) = values[ix].value;
  // test against stupid mixing matrices  
  for(unsigned int ix=0;ix<matrix->size().first;++ix) {
    Complex sum(0.);
    for(unsigned int iy=0;iy<matrix->size().second;++iy) {
      sum += norm((*matrix)(ix,iy));
    }
    if(abs(sum-1.)>1e-4) {
      cerr << "The sum of the mod squares of row " << ix+1
	   << " of the " << name << " block does not sum to 1. \n"
	   << "sum = " << sum.real() << ". We strongly suggest you check your SLHA file.\n";
    }
  }

  vector<long> ids;
  if(name == "nmix") {
    ids.resize(4);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035; 
  }
  else if(name == "nmnmix") {
    ids.resize(5);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035;
    ids[4] = 1000045;
  }
  else if(name == "rvnmix") {
    ids.resize(7);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035;
    if(!majoranaNeutrinos()) {
      ids[4] = 12; ids[5] = 14; ids[6] = 16;
    }
    else { 
      ids[4] = 17; ids[5] = 18; ids[6] = 19;
    }
  }
  else if(name == "rpvmix") {
    ids.resize(5);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035;
    ids[4] = 1000045;
  }
  else if(name == "umix" || name == "vmix") {
    ids.resize(2);
    ids[0] = 1000024; ids[1] = 1000037;
  }
  else if(name == "rvumix" || name == "rvvmix" ) {
    ids.resize(5);
    ids[3] = 1000024; ids[4] = 1000037; 
    ids[0] = -11; ids[1] = -13; ids[2] = -15;
  }
  else if(name == "stopmix") {
    ids.resize(2);
    ids[0] = 1000006; ids[1] = 2000006;
  }
  else if(name == "sbotmix") {
    ids.resize(2);
    ids[0] = 1000005; ids[1] = 2000005;
  }
  else if(name == "staumix") {
    ids.resize(2);
    ids[0] = 1000015; ids[1] = 2000015;
  }
  else if(name == "nmhmix") {
    ids.resize(3);
    ids[0] = 25; ids[1] = 35; ids[2] = 45;
  }
  else if(name == "nmamix") {
    ids.resize(2);
    ids[0] = 36; ids[1] = 46;
  }
  else if(name == "rvamix") {
    // some programs, i.e. spheno include the goldstone even though the standard says they shouldn't
    // obeys standard
    if(size.first==4) {
      ids.resize(4);
    }
    // includes goldstone
    else if(size.first==5) {
      ids.resize(5);
      ids[4] = 0;
    }
    else {
      throw Exception() << "SusyBase::createMixingMatrix() "
			<< "rvamix matrix must have either 4 rows (obeying standard) "
			<< "or 5 including the goldstone, not " << size.first << Exception::runerror;
    }
    ids[0] = 36     ; ids[1] = 1000017;
    ids[2] = 1000018; ids[3] = 1000019;
  }
  else if(name == "rvhmix") {
    ids.resize(5);
    ids[0] = 25; ids[1] = 35;
    ids[2] = 1000012; ids[3] = 1000014; ids[4] = 1000016;
  }
  else if(name == "rvlmix") {
    // some programs, i.e. spheno include the goldstone even though the standard says they shouldn't
    // obeys standard
    if(size.first==7) {
      ids.resize(7);
    }
    // includes goldstone
    else if(size.first==8) {
      ids.resize(8);
      ids[7] = 0;
    }
    else {
      throw Exception() << "SusyBase::createMixingMatrix() "
			<< "rvlmix matrix must have either 7 rows (obeying standard) "
			<< "or 8 including the goldstone, not " << size.first << Exception::runerror;
    }
    ids[0] = 37;
    ids[1] = -1000011;
    ids[2] = -1000013;
    ids[3] = -1000015;
    ids[4] = -2000011;
    ids[5] = -2000013;
    ids[6] = -2000015;
  }
  else if(name == "usqmix" ) {
    ids.resize(6);
    ids[0] = 1000002; ids[1] = 1000004; ids[2] = 1000006;
    ids[3] = 2000002; ids[4] = 2000004; ids[5] = 2000006;
  }
  else if(name == "dsqmix" ) {
    ids.resize(6);
    ids[0] = 1000001; ids[1] = 1000003; ids[2] = 1000005;
    ids[3] = 2000001; ids[4] = 2000003; ids[5] = 2000005;
  }
  else {
    throw SetupException() << "SusyBase::createMixingMatrix() "
			   << "cannot find correct title for mixing matrix "
			   << name << Exception::runerror;
  }
  matrix->setIds(ids);
}

void SusyBase::resetRepositoryMasses() {
  map<string,ParamMap>::const_iterator fit=parameters_.find("mass");
  if(fit==parameters_.end()) 
    throw Exception() << "BLOCK MASS not found in input file"
		      << " can't set masses of SUSY particles"
		      << Exception::runerror;
  ParamMap theMasses = fit->second;
  for(ParamMap::iterator it = theMasses.begin(); it != theMasses.end(); 
      ++it) {
    long id = it->first;
    map<long,long>::const_iterator dit = idMap().find(id);
    if(dit!=idMap().end()) id = dit->second;
    double mass = it->second;
    //a negative mass requires an adjustment to the 
    //associated mixing matrix by a factor of i
    if(mass < 0.0) adjustMixingMatrix(id);
    PDPtr part = getParticleData(id);
    if(!part) throw SetupException() 
      << "SusyBase::resetRepositoryMasses() - Particle with PDG code " << id  
      << " not found." << Exception::warning;
    if(abs(id)<=5||abs(id)==23||abs(id)==24||
       (abs(id)>=11&&abs(id)<=16)) {
      if(allowedToResetSMMasses_) {
	cerr << "SusyBase::resetRepositoryMasses() Resetting mass of " 
	     << part->PDGName() << " using SLHA "
	     << "file,\nthis can affect parts of the Standard Model simulation and"
	     << " is strongly discouraged.\n";
	// reset the masses
	resetMass(it->first,it->second*GeV,part);
      }
      else {
	cerr << "SusyBase::resetRepositoryMasses() You have tried to Reset the mass of " 
	     << part->PDGName() << " using an SLHA "
	     << "file,\nthis can affect parts of the Standard Model simulation and"
	     << " is not allowed by default.\n If you really want to be this stupid"
	     << " set AllowedToResetSMMasses to Yes for this model.\n";
      }
    }
    // reset the masses for not SM particles
    else {
      resetMass(it->first,it->second*GeV,part);
    }
    // switch on gravitino interactions?
    gravitino_ |= id== ParticleID::SUSY_Gravitino;
  }
  theMasses.clear();
}

void SusyBase::adjustMixingMatrix(long id) {
  //get correct mixing matrix
  switch(id) {
  case 1000021 :
    gluinoPhase_ = Complex(0.,1.);
    break;
  case 1000022 :
  case 1000023 :
  case 1000025 :
  case 1000035 : 
  case 1000045 :
  case 12 : case 17 : 
  case 14 : case 18 : 
  case 16 : case 19 :
    if(NMix_) {
      if(id>20||(id<=16&&NMix_->size().first>4))
	 NMix_->adjustPhase(id);
    }
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The neutralino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    break;
  case 1000024 :
  case 1000037 : 
    if(UMix_)
      UMix_->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The U-Type chargino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    if(VMix_) VMix_->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The V-Type chargino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    break;
  default : 
    throw SetupException() 
      << "SusyBase::adjustMixingMatrix - Trying to adjust mixing matrix "
      << "phase for a particle that does not have a mixing matrix "
      << "associated with it. " << id << " must have a negative mass in "
      << "the spectrum file, this should only occur for particles that mix."
      << Exception::runerror;
  }
}

void SusyBase::createMixingMatrices() {
  map<string,pair<MatrixSize, MixingVector > >::const_iterator it;
  for(it=mixings_.begin();it!=mixings_.end();++it) {
    string name=it->first;
    // create the gaugino mixing matrices
    if(name == "nmix" || name == "nmnmix" || name == "rvnmix" ) {
      createMixingMatrix(NMix_,name,it->second.second,it->second.first);
    }
    else if (name == "umix" || name == "rvumix" ) {
      createMixingMatrix(UMix_,name,it->second.second,it->second.first);
    }
    else if (name == "vmix" || name == "rvvmix" ) {
      createMixingMatrix(VMix_,name,it->second.second,it->second.first);
    }
  }
}

void SusyBase::extractParameters(bool checkmodel) {
  map<string,ParamMap>::const_iterator pit;
  ParamMap::const_iterator it;
  // try and get tan beta from hmin and extpar first
  // extract tan beta
  tanBeta_ = -1.;
  if(tanBeta_<0.) {
    pit=parameters_.find("hmix");
    if(pit!=parameters_.end()) {
      it = pit->second.find(2);
      if(it!=pit->second.end()) tanBeta_ = it->second;
    }
  }
  if(tanBeta_<0.) {
    pit=parameters_.find("extpar");
    if(pit!=parameters_.end()) {
      it = pit->second.find(25);
      if(it!=pit->second.end()) tanBeta_ = it->second;
    }
  }
  // otherwise from minpar
  if(tanBeta_<0.) {
    pit=parameters_.find("minpar");
    if(pit!=parameters_.end()) { 
      it = pit->second.find(3);
      if(it!=pit->second.end()) tanBeta_ = it->second;
    }
  }
  if(tanBeta_<0.) 
    throw Exception() << "SusyBase::extractParameters() "
		      << "Can't find tan beta in BLOCK MINPAR"
		      << " or BLOCK EXTPAR " << Exception::runerror;
  if(tanBeta_==0.)
    throw Exception() << "Tan(beta) = 0 in SusyBase::extractParameters()"
		      << Exception::runerror;
  // extract parameters from hmix
  pit=parameters_.find("hmix");
  if(pit==parameters_.end()) {
    if(generator())
      generator()->logWarning(Exception("SusyBase::extractParameters() BLOCK HMIX not found setting mu to zero\n",
					Exception::warning));
    else
      cerr << "SusyBase::extractParameters() BLOCK HMIX not found setting mu to zero\n";
     mu_=ZERO;
  }
  else {
    mu_=findValue(pit,1,"HMIX","mu")*GeV;
  }
  pit = parameters_.find("msoft");
  if( pit == parameters_.end() )
    throw Exception() << "SusyBase::extractParameters() "
		      << "BLOCK MSOFT not found in " 
		      << "SusyBase::extractParameters()"
		      << Exception::runerror;
  M1_    = findValue(pit,1 ,"MSOFT","M_1"   )*GeV;
  M2_    = findValue(pit,2 ,"MSOFT","M_2"   )*GeV;
  M3_    = findValue(pit,3 ,"MSOFT","M_3"   )*GeV;
  mH12_  = findValue(pit,21,"MSOFT","m_H1^2")*GeV2;
  mH22_  = findValue(pit,22,"MSOFT","m_H2^2")*GeV2;
  meL_   = findValue(pit,31,"MSOFT","M_eL"  )*GeV;
  mmuL_  = findValue(pit,32,"MSOFT","M_muL" )*GeV;
  mtauL_ = findValue(pit,33,"MSOFT","M_tauL")*GeV; 
  meR_   = findValue(pit,34,"MSOFT","M_eR"  )*GeV;
  mmuR_  = findValue(pit,35,"MSOFT","M_muR" )*GeV;
  mtauR_ = findValue(pit,36,"MSOFT","M_tauR")*GeV; 
  mq1L_  = findValue(pit,41,"MSOFT","M_q1L" )*GeV;
  mq2L_  = findValue(pit,42,"MSOFT","M_q2L" )*GeV;
  mq3L_  = findValue(pit,43,"MSOFT","M_q3L" )*GeV; 
  muR_   = findValue(pit,44,"MSOFT","M_uR"  )*GeV;
  mcR_   = findValue(pit,45,"MSOFT","M_cR"  )*GeV;
  mtR_   = findValue(pit,46,"MSOFT","M_tR"  )*GeV;
  mdR_   = findValue(pit,47,"MSOFT","M_dR"  )*GeV;
  msR_   = findValue(pit,48,"MSOFT","M_sR"  )*GeV;
  mbR_   = findValue(pit,49,"MSOFT","M_bR"  )*GeV;
  if(checkmodel) {
    throw Exception() << "The SusyBase class should not be used as a "
		      << "Model class, use one of the models which inherit"
		      << " from it" << Exception::runerror;
  }
}
#line 1 "./MSSM.cc"
// -*- C++ -*-
//
// MSSM.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MSSM class.
//

#include "MSSM.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void MSSM::persistentOutput(PersistentOStream & os) const {
  os << theStopMix << theSbotMix << theStauMix << theAlpha 
     << ounit(theAtop,GeV) << ounit(theAbottom,GeV) << ounit(theAtau,GeV) 
     << theHiggsMix << HiggsAMix_ << HiggsPMix_ << createDiagonalMixing_;
}

void MSSM::persistentInput(PersistentIStream & is, int) {
  is >> theStopMix >> theSbotMix >> theStauMix >> theAlpha 
     >> iunit(theAtop,GeV) >> iunit(theAbottom,GeV) >> iunit(theAtau,GeV) 
     >> theHiggsMix >> HiggsAMix_ >> HiggsPMix_ >>  createDiagonalMixing_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MSSM,SusyBase>
describeMSSM("Herwig::MSSM", "HwSusy.so");

void MSSM::Init() {

  static ClassDocumentation<MSSM> documentation
    ("The MSSM class is the base class for the MSSM model.",
     "MSSM Feynman rules were taken from \\cite{Haber:1984rc,Gunion:1984yn}.",
     " %\\cite{Haber:1984rc}\n"
     "\\bibitem{Haber:1984rc}\n"
     "  H.~E.~Haber and G.~L.~Kane,\n"
     "  %``The Search For Supersymmetry: Probing Physics Beyond The Standard Model,''\n"
     "  Phys.\\ Rept.\\  {\\bf 117}, 75 (1985).\n"
     "  %%CITATION = PRPLC,117,75;%%\n"
     "%\\cite{Gunion:1984yn}\n"
     "\\bibitem{Gunion:1984yn}\n"
     "  J.~F.~Gunion and H.~E.~Haber,\n"
     "  %``Higgs Bosons In Supersymmetric Models. 1,''\n"
     "  Nucl.\\ Phys.\\  B {\\bf 272}, 1 (1986)\n"
     "  [Erratum-ibid.\\  B {\\bf 402}, 567 (1993)].\n"
     "  %%CITATION = NUPHA,B272,1;%%\n"
    );

  static Switch<MSSM,bool> interfaceCreateDiagonalMixingMatrices
    ("CreateDiagonalMixingMatrices",
     "Create diagonal stop, sbottom and stau mixings if not present.",
     &MSSM::createDiagonalMixing_, false, false, false);
  static SwitchOption interfaceCreateDiagonalMixingMatricesNo
    (interfaceCreateDiagonalMixingMatrices,
     "No",
     "Don't create them",
     false);
  static SwitchOption interfaceCreateDiagonalMixingMatricesYes
    (interfaceCreateDiagonalMixingMatrices,
     "Yes",
     "Create them if needed",
     true);

}

void MSSM::createMixingMatrices() {
  useMe();
  map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    // create the stop, sbottom and stau mixing matrices
    if(name == "stopmix"  ){
      createMixingMatrix(theStopMix,name,it->second.second,it->second.first);
    }
    else if (name == "sbotmix" ) {
      createMixingMatrix(theSbotMix,name,it->second.second,it->second.first);
    }
    else if (name == "staumix") {
      createMixingMatrix(theStauMix,name,it->second.second,it->second.first);
    }
    // Higgs mixing matrix in extended models
    else if (name == "nmhmix" || name == "rvhmix") {
      createMixingMatrix(theHiggsMix,name,it->second.second,it->second.first);
    }
  }
  // create stop, sbottom and stau mixing if needed and absent
  if(createDiagonalMixing_) {
    // stop
    if(!theStopMix) {
      theStopMix = new_ptr(MixingMatrix(2,2));
      (*theStopMix)(0,0) = 1.;
      (*theStopMix)(1,1) = 1.;
    }
    // sbottom
    if(!theSbotMix) {
      theSbotMix = new_ptr(MixingMatrix(2,2));
      (*theSbotMix)(0,0) = 1.;
      (*theSbotMix)(1,1) = 1.;
    }
    // stau
    if(!theStauMix) {
      theStauMix = new_ptr(MixingMatrix(2,2));
      (*theStauMix)(0,0) = 1.;
      (*theStauMix)(1,1) = 1.;
    }
  }
  // neutral higgs mixing if not already set
  if(!theHiggsMix) {
    MixingVector hmix;
    hmix.push_back(MixingElement(2,1, cos(theAlpha)));
    hmix.push_back(MixingElement(2,2, sin(theAlpha)));
    hmix.push_back(MixingElement(1,1,-sin(theAlpha)));
    hmix.push_back(MixingElement(1,2, cos(theAlpha)));
    vector<long> ids(2);
    ids[0] = 25; ids[1] = 35;
    theHiggsMix = new_ptr(MixingMatrix(2,2));
    (*theHiggsMix).setIds(ids);
    for(unsigned int ix=0; ix < hmix.size(); ++ix)
      (*theHiggsMix)(hmix[ix].row-1,hmix[ix].col-1) = hmix[ix].value;
    hmix.clear();
    double beta = atan(tanBeta());
    hmix.push_back(MixingElement(1,1,sin(beta)));
    hmix.push_back(MixingElement(1,2,cos(beta)));
    ids.clear();
    ids.resize(1,37);
    HiggsPMix_ = new_ptr(MixingMatrix(1,2));
    (*HiggsPMix_).setIds(ids);
    for(unsigned int ix=0; ix < hmix.size(); ++ix)
      (*HiggsPMix_)(hmix[ix].row-1,hmix[ix].col-1) = hmix[ix].value;
    ids.clear();
    ids.resize(1,36);
    HiggsAMix_ = new_ptr(MixingMatrix(1,2));
    (*HiggsAMix_).setIds(ids);
    for(unsigned int ix=0; ix < hmix.size(); ++ix)
      (*HiggsAMix_)(hmix[ix].row-1,hmix[ix].col-1) = hmix[ix].value;
  }
  // base class for neutralinos and charginos
  SusyBase::createMixingMatrices();
}

void MSSM::adjustMixingMatrix(long id) {
  switch (id) {
  case 1000006 :
  case 2000006 :
    if(theStopMix)
      theStopMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The stop mixing matrix pointer is null!" 
			     << Exception::runerror;
    break;
  case 1000005 :
  case 2000005 :
    if(theSbotMix)
      theSbotMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The sbottom mixing matrix pointer is null!" 
			     << Exception::runerror;
    break;
  case 1000015 :
  case 2000015 :
    if(theStauMix)
      theStauMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The stau mixing matrix pointer is null!" 
			     << Exception::runerror;
    break;
  default :
    SusyBase::adjustMixingMatrix(id);
    break;
  }
}

void MSSM::extractParameters(bool checkmodel) {
  map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
  // trilinear couplings
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    MixingVector::const_iterator vit;
    if(name=="au") {
      theAtop=ZERO;
      for(vit=it->second.second.begin();vit!=it->second.second.end();++vit) {
	if(vit->row==3&&vit->col==3) theAtop=vit->value*GeV;
      }
    }
    else if(name=="ad") {
      theAbottom=ZERO;
      for(vit=it->second.second.begin();vit!=it->second.second.end();++vit) {
	if(vit->row==3&&vit->col==3) theAbottom=vit->value*GeV;
      }
    }
    else if(name=="ae") {
      theAtau=ZERO;
      for(vit=it->second.second.begin();vit!=it->second.second.end();++vit) {
	if(vit->row==3&&vit->col==3) theAtau=vit->value*GeV;
      }
    }
  }
  // neutralino and chargino paramters in the base class
  SusyBase::extractParameters(false);
  // check the model
  map<string,ParamMap>::const_iterator pit;
  pit=parameters().find("modsel");
  if(pit==parameters().end()) return;
  // nmssm or mssm
  ParamMap::const_iterator jt;
  jt = pit->second.find(3);
  int inmssm = jt!=pit->second.end() ? int(jt->second) : 0; 
  // RPV
  jt = pit->second.find(4);
  int irpv = jt!=pit->second.end() ? int(jt->second) : 0;
  // CPV
  jt = pit->second.find(5);
  int icpv = jt!=pit->second.end() ? int(jt->second) : 0;
  // flavour violation
  jt = pit->second.find(6);
  int ifv = jt!=pit->second.end() ? int(jt->second) : 0;
  // the higgs mixing angle 
  theAlpha=0.;
  bool readAlpha = false;
  pit=parameters().find("alpha");
  if(pit!=parameters().end()) {
    ParamMap::const_iterator it = pit->second.find(1);
    if(it!=pit->second.end()) {
      readAlpha = true;
      theAlpha=it->second;
    }
  }
  if(inmssm==0&&irpv==0&&!readAlpha) 
    throw Exception() << "In the MSSM model BLOCK ALPHA which must be"
		      << " present in the SLHA file is missing"
		      << Exception::runerror;
  if(checkmodel) {
    if(inmssm!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				    << " used but NMSSM read in " 
				    << Exception::runerror;
    if(irpv!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				  << " used but RPV read in " 
				  << Exception::runerror; 
    if(icpv!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				  << " used but CPV read in " 
				  << Exception::runerror; 
    if(ifv!=0) throw Exception() << "R-parity, CP and flavour conserving MSSM model"
				 << " used but flavour violation read in " 
				 << Exception::runerror;
  }
}
#line 1 "./MixingMatrix.cc"
// -*- C++ -*-
//
// MixingMatrix.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MixingMatrix class.
//

#include "MixingMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void MixingMatrix::persistentOutput(PersistentOStream & os) const {
  os << mixingMatrix_ << ids_ << size_;
}

void MixingMatrix::persistentInput(PersistentIStream & is, int) {
  is >> mixingMatrix_ >> ids_ >> size_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MixingMatrix,Interfaced>
describeMixingMatrix("Herwig::MixingMatrix", "HwSusy.so");

void MixingMatrix::Init() {

  static ClassDocumentation<MixingMatrix> documentation
    ("The MixingMatrix class implements the storage of the SUSY mixing "
     "matrices.");

}

void MixingMatrix::adjustPhase(long id) { 
  unsigned int irow(0);
  while(irow < size().first && ids_[irow] != id) 
    ++irow;
  for(unsigned int c = 0; c < size_.second; ++c)
    mixingMatrix_[irow][c] *= Complex(0., 1.);
}

ostream & Herwig::operator<<(ostream & os,const MixingMatrix & mix) {
  for(unsigned int ix=0;ix<mix.size().first;++ix) {
    for(unsigned int iy=0;iy<mix.size().second;++iy) {
      os << mix(ix,iy) << "\t";
    }
    os << "\n";
  }
  return os;
}
#line 1 "./SSCFSVertex.cc"
// -*- C++ -*-
//
// SSCFSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSCFSVertex class.
//

#include "SSCFSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSCFSVertex::SSCFSVertex(): _sb(0.),_cb(0.),_mw(ZERO),
			    _q2last(0.*GeV2), _couplast(0.),
			    _leftlast(0.),_rightlast(0.),
			    _id1last(0), _id2last(0), _id3last(0),
			    yukawa_(1) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SSCFSVertex::doinit() {
  long chargino[2] = {1000024, 1000037};
  for(unsigned int ic = 0; ic < 2; ++ic) {
    //quarks 
    for(long ix = 1; ix < 7; ++ix) {
      if( ix % 2 == 0 ) {
	addToList(-chargino[ic],ix,-(999999+ix));
	
	addToList(-chargino[ic],ix,-(1999999+ix));
	
	addToList(-ix,chargino[ic],(999999+ix));
	
	addToList(-ix,chargino[ic],(1999999+ix));
      }
      else {
	addToList(-chargino[ic],-ix,(1000001+ix));

	addToList(-chargino[ic],-ix,2000001+ix);

	addToList(chargino[ic],ix,-(1000001+ix));

	addToList(chargino[ic],ix,-(2000001+ix));
      }
    }
    //leptons
    for(long ix = 11; ix < 17; ++ix) {
      if( ix % 2 == 0 ) {
	addToList(-chargino[ic],ix,-(999999+ix));
      
	addToList(-chargino[ic],ix,-(1999999+ix));

	addToList(-ix,chargino[ic],(999999+ix));

	addToList(-ix,chargino[ic],(1999999+ix));	
      }
      else {
	addToList(-chargino[ic],-ix,1000001+ix);

	addToList(chargino[ic],ix,-(1000001+ix));
      }
    }
  } 
  FFSVertex::doinit();
  _theSS = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  //mixing matrices
  _stop = _theSS->stopMix();
  _sbot = _theSS->sbottomMix();
  _stau = _theSS->stauMix();
  _umix = _theSS->charginoUMix();
  _vmix = _theSS->charginoVMix();

  if(!_stop || !_stau || !_sbot || !_umix || !_vmix)
    throw InitException() << "SSCFSVertex:: doinit  - " 
			  << "A mixing matrix pointer is null."
			  << " stop: " << _stop << " sbottom: " << _sbot
			  << " stau: " << _stau << " U: " << _umix
			  << " V:" << _vmix
			  << Exception::abortnow;

  _mw = getParticleData(24)->mass();
  double tb = _theSS->tanBeta();
  _sb = tb/sqrt(1 + sqr(tb));
  _cb = sqrt(1.- sqr(_sb));
}


void SSCFSVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS << _sb << _cb << ounit(_mw,GeV) << _stop 
     << _sbot << _stau << _umix << _vmix << yukawa_;
}

void SSCFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS  >> _sb >> _cb >> iunit(_mw,GeV) >> _stop
     >> _sbot >> _stau >> _umix >> _vmix >> yukawa_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSCFSVertex,FFSVertex>
describeHerwigSSCFSVertex("Herwig::SSCFSVertex", "HwSusy.so");

void SSCFSVertex::Init() {

  static ClassDocumentation<SSCFSVertex> documentation
    ("The implementation of the coupling of the charginos to fermion-"
     "sfermions.");

  static Switch<SSCFSVertex,unsigned int> interfaceYukawa
    ("Yukawa",
     "Whether or not to include the Yukawa type couplings",
     &SSCFSVertex::yukawa_, true, false, false);
  static SwitchOption interfaceYukawaYes
    (interfaceYukawa,
     "Yes",
     "Include the terms",
     1);
  static SwitchOption interfaceYukawaNo
    (interfaceYukawa,
     "No",
     "Don't include them",
     0);
  static SwitchOption interfaceYukawa3rdGen
    (interfaceYukawa,
     "ThirdGeneration",
     "Only include them for the third generation",
     2);

}

void SSCFSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  long isc(abs(part3->id())), ism(abs(part1->id())), 
    ichg(abs(part2->id()));
  tcPDPtr smfermion = part1;
  if( ism / 1000000 == 1 )  {
    swap( ism, ichg);
    smfermion = part2;
  }
  //overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _q2last=q2;
    _couplast = -weakCoupling(q2);
  }
  norm(_couplast);


  if( ichg != _id1last || ism != _id2last || isc != _id3last ) {
    _id1last = ichg;
    _id2last = ism;
    _id3last = isc;
    // determine chargino and sfermion eigenstates
    unsigned int alpha(isc/1000000 - 1);
    unsigned int ch = (ichg == 1000024 ) ? 0 : 1;

    Complex ul1 = (*_umix)(ch,0);
    Complex ul2 = (*_umix)(ch,1);
    Complex vl1 = (*_vmix)(ch,0);
    Complex vl2 = (*_vmix)(ch,1);

    if( ism >= 11 && ism <= 16 ) {
      long lept = ( ism % 2 == 0 ) ? ism - 1 : ism;
      double y = 0.;
      if(yukawa_==1 || (lept==15 && yukawa_==2))
	y = double(_theSS->mass(q2, getParticleData(lept))/_mw/sqrt(2)/_cb);


      if( ism == 12 || ism == 14 ) {
	_leftlast = Complex(0., 0.);
	if( alpha == 0 )
	  _rightlast = ul1;
	else
	  _rightlast = -y*ul2;
      }
      else if( ism == 16 ) {
	_leftlast = Complex(0., 0.);
	_rightlast = ul1*(*_stau)(alpha, 0) - y*(*_stau)(alpha,1)*ul2;
      }
      else if( ism == 11 || ism == 13 || ism == 15 ) {
	_leftlast = -y*conj(ul2);
	_rightlast = vl1;
      }
    }
    else {
      double yd(0.), yu(0.);
      if(yukawa_==1 || ((ism==5 || ism==6 ) && yukawa_==2)) {
	if( ism % 2 == 0) {
	  yu = _theSS->mass(q2, getParticleData(ism))/_mw/sqrt(2)/_sb;
	  yd = _theSS->mass(q2, getParticleData(ism - 1))/_mw/sqrt(2)/_cb;
	}
	else {
	  yu = _theSS->mass(q2, getParticleData(ism + 1))/_mw/sqrt(2)/_sb;
	  yd = _theSS->mass(q2, getParticleData(ism))/_mw/sqrt(2)/_cb;
	}
      }
      //heavy quarks
      if( ism == 5 ) {
	_leftlast = -yd*conj(ul2)*(*_stop)(alpha,0);
	_rightlast = vl1*(*_stop)(alpha, 0) - yu*vl2*(*_stop)(alpha,1);
      }
      else if( ism == 6 ) {
	_leftlast = -yu*conj(vl2)*(*_sbot)(alpha,0);
	_rightlast = ul1*(*_sbot)(alpha, 0) - yd*ul2*(*_sbot)(alpha,1);
      }
      else {
	if( alpha == 0 ) {
	  _leftlast = (ism % 2 == 0) ? -yu*conj(vl2) : -yd*conj(ul2);
	  _rightlast = (ism % 2 == 0) ? ul1 : vl1;
	}
	else {
	  _leftlast = Complex(0.);
	  _rightlast = (ism % 2 == 0) ? -yd*ul2 : -yu*vl2;
	}
      }
    }
  }//end of coupling calculation

  //determine the helicity order of the vertex
  if( smfermion->id() < 0 ) {
    left(conj(_rightlast));
    right(conj(_leftlast));
  }
  else {
    left(_leftlast);
    right(_rightlast);
  }
}
#line 1 "./SSGFSVertex.cc"
// -*- C++ -*-
//
// SSGFSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGFSVertex class.
//

#include "SSGFSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGFSVertex::SSGFSVertex() :_q2last(0.*GeV2),_couplast(0.), 
			    _id1last(0), _id2last(0) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TFUND);
}

void SSGFSVertex::persistentOutput(PersistentOStream & os) const {
  os << _stop << _sbottom << gluinoPhase_;
}

void SSGFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _stop >> _sbottom >> gluinoPhase_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSGFSVertex,FFSVertex>
describeHerwigSSGFSVertex("Herwig::SSGFSVertex", "HwSusy.so");

void SSGFSVertex::Init() {

  static ClassDocumentation<SSGFSVertex> documentation
    ("The SSGFSVertex implements coupling of the gluinos to the "
     "squarks and quarks");

}

void SSGFSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2, tcPDPtr part3) {
  tcPDPtr ferm;
  long isc(0);
  if(abs(part1->id()) == 1000021) {
    if(part2->iSpin() == PDT::Spin1Half) {
      ferm = part2;
      isc = abs(part3->id());
    }
    else {
      ferm = part3;
      isc = abs(part2->id());
    }
  }
  else if(abs(part2->id()) == 1000021) {
    if(part1->iSpin() == PDT::Spin1Half) {
      ferm = part1;
      isc = abs(part3->id());
    }
    else {
      ferm = part3;
      isc = abs(part1->id());
    }
  }
  else if(abs(part3->id()) == 1000021) {
    if(part1->iSpin() == PDT::Spin1Half) {
      ferm = part1;
      isc = abs(part2->id());
    }
    else {
      ferm = part2;
      isc = abs(part1->id());
    }
  }
  else throw HelicityConsistencyError()
    << "SSGFSVertex::setCoupling() - There is no gluino in this vertex!"
    << part1->id() << " " << part2->id() << " " << part3->id()
    << Exception::runerror;
  long iferm = abs(ferm->id());
  assert(iferm >=1 && iferm <=6);

  if(q2 != _q2last  || _couplast==0.) {
    _couplast = -strongCoupling(q2)*sqrt(2.);
    _q2last = q2;
  }
  if(iferm != _id1last || isc != _id2last) { 
    _id1last = iferm;
    _id2last = isc;
    unsigned int eig = (isc/1000000) - 1;
    if(iferm == 6) {
      _leftlast  = -(*_stop)(eig,1)*conj(gluinoPhase_);
      _rightlast =  (*_stop)(eig,0)*     gluinoPhase_ ;
    }
    else if(iferm == 5){
      _leftlast  = -(*_sbottom)(eig,1)*conj(gluinoPhase_);
      _rightlast =  (*_sbottom)(eig,0)*     gluinoPhase_ ;
    }
    else {
      if(eig == 0) { 
	_leftlast  =  0.;
	_rightlast =  gluinoPhase_;
      }
      else {
	_leftlast  = -conj(gluinoPhase_);
	_rightlast =  0.;
      }
    }
  }
  norm(_couplast);
  //arrange l/r couplings
  if(ferm->id() < 0) {
    left(conj(_rightlast));
    right(conj(_leftlast));
  }
  else {
    left(_leftlast);
    right(_rightlast);
  }
}

void SSGFSVertex::doinit() {
  for(long ix=1;ix<7;++ix) {
    addToList(1000021, ix, -(ix+1000000));
    addToList(1000021, ix, -(ix+2000000));
    addToList(1000021, -ix, (ix+1000000));
    addToList(1000021, -ix, (ix+2000000));
  }
  FFSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());

  _stop = model->stopMix();
  _sbottom = model->sbottomMix();
  gluinoPhase_ = model->gluinoPhase();
  if(!_stop || !_sbottom)
    throw InitException() << "SSGFSVertex::doinit() - "
			  << "There is a null mixing matrix pointer. "
			  << "stop: " << _stop << " sbottom: " << _sbottom 
			  << Exception::abortnow;
}
#line 1 "./SSHSFSFVertex.cc"
// -*- C++ -*-
//
// SSHSFSFVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSHSFSFVertex class.
//

#include "SSHSFSFVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSHSFSFVertex::SSHSFSFVertex() : theMix(3), theTriC(9, complex<Energy>(ZERO)), 
				 theSinA(0.0),
				 theCosA(0.0), theSinB(0.0), theCosB(0.0),
				 theTanB(0.0), theSinAB(0.0), theCosAB(0.0),
				 theMw(ZERO), theMz(ZERO), theMu(ZERO), 
				 theSw(0.0), theCw(0.0), theCoupLast(ZERO),
				 theq2Last(ZERO), theHLast(0), theSF1Last(0),
				 theSF2Last(0) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SSHSFSFVertex::doinit() {
  int higgs = 25;
  //h0,H0
  for(unsigned int i = 0; i < 2; ++i) {
    if( i == 1 ) higgs = 35;
    //quarks
    for(unsigned int j = 1; j < 7; ++j) {
      long lj = 1000000 + j;
      long rj = 2000000 + j;
      //LLbar
      addToList(higgs,lj,-lj);
      //RRbar
      addToList(higgs,rj,-rj);
      //LRbar
      addToList(higgs,lj,-rj);
      //RLbar
      addToList(higgs,rj,-lj);
    }
    for(unsigned int j = 11; j < 17; ++j) {
      long lj = 1000000 + j;
      long rj = 2000000 + j;
      addToList(higgs,lj,-lj);
      if( j % 2 != 0) {
	addToList(higgs,rj,-rj);
	//LRbar
	addToList(higgs,lj,-rj);
	//RLbar
	addToList(higgs,rj,-lj);
      }
    }
  }
  //A0
  for(unsigned int j = 1; j < 7; ++j) {
    long lj = 1000000 + j;
    long rj = 2000000 + j;
    //LRbar
    addToList(36,lj,-rj);
    //RLbar
    addToList(36,rj,-lj);
  }
  for(unsigned int j = 11; j < 17; j += 2) {
    long lj = 1000000 + j;
    long rj = 2000000 + j;
    addToList(36,lj,-rj);
    addToList(36,rj,-lj);
  }
  //outgoing H+
  for(long ii = 2; ii < 7; ii += 2) {
    //LL
    addToList(37, 999999 + ii, -1000000 - ii);
    //RR
    addToList(37, 1999999 + ii, -2000000 - ii);
    //RL
    addToList(37, 1999999 + ii, -1000000 - ii);
    //LR
    addToList(37, 999999 + ii, -2000000 - ii);
  }
  for(long ii = 11; ii < 17; ii += 2) {
    addToList(37, 1000000 + ii, -1000001 - ii);
    addToList(37, 2000000 + ii, -1000001 - ii);
  }
  //outgoing H-
  for(long ii = 2; ii < 7; ii += 2) {
    //LL
    addToList(-37, 1000000 + ii, -999999 - ii);
    //RR
    addToList(-37, 2000000 + ii, -1999999 - ii);
    //RL
    addToList(-37, 1000000 + ii, -1999999 - ii);
    //LR
    addToList(-37, 2000000 + ii, -999999 - ii);
  }
  for(long ii = 11; ii < 17; ii += 2) {
    addToList(-37, 1000001 + ii, -1000000 - ii);
    addToList(-37, 1000001 + ii, -2000000 - ii);}
  SSSVertex::doinit();
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() << "SSHSFSFVertex::doinit - A problem occurred"
			  << "while trying to cast the SM pointer to "
			  << "a Susy one." << Exception::abortnow;
  //mixing vector should have been sized correctly already
  assert( theMix.size() == 3 );
  theMix[0] = theMSSM->stopMix();
  theMix[1] = theMSSM->sbottomMix();
  theMix[2] = theMSSM->stauMix();

  if(!theMix[0] || !theMix[1] || !theMix[2])
    throw InitException() << "SSHSFSFVertex::doinit -  "
			  << "A null mixing matrix pointer. stop: " << theMix[0] 
			  << " sbottom: " << theMix[1] << " stau: " << theMix[2]
			  << Exception::abortnow;
  //trilinear vector should have been sized correctly already
  assert( theTriC.size() == 9 );
  //vector has been zeroed in constructor
  theTriC[4]=theMSSM->bottomTrilinear().real();
  theTriC[5]=theMSSM->topTrilinear().real();
  theTriC[8]=theMSSM->tauTrilinear().real();
  //Masses
  theMw = getParticleData(ParticleID::Wplus)->mass();
  theMz = getParticleData(ParticleID::Z0)->mass();
  //parameters
  theSinA = sin(theMSSM->higgsMixingAngle());
  theCosA = cos(theMSSM->higgsMixingAngle());
  theTanB = theMSSM->tanBeta();
  theMu = theMSSM->muParameter();
  theSinB = theTanB/sqrt(1. + sqr(theTanB));
  theCosB = sqrt( 1. - sqr(theSinB) );
  theSinAB = theSinA*theCosB + theCosA*theSinB;
  theCosAB = theCosA*theCosB - theSinA*theSinB;

  theSw = sqrt(sin2ThetaW());
  theCw = sqrt(1. - sin2ThetaW());
}

void SSHSFSFVertex::persistentOutput(PersistentOStream & os) const {
  os << theMix << theSinA << theCosA << theSinB << theMSSM
     << theCosB << theTanB << ounit(theMu, GeV) << theSinAB << theCosAB 
     << ounit(theMw,GeV) << ounit(theMz,GeV) << theSw << theCw 
     << ounit(theTriC,GeV);

}

void SSHSFSFVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMix >>  theSinA >> theCosA >> theSinB >> theMSSM
     >> theCosB >> theTanB >> iunit(theMu, GeV) >> theSinAB >> theCosAB 
     >> iunit(theMw,GeV) >> iunit(theMz,GeV) >> theSw >> theCw
     >> iunit(theTriC,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSHSFSFVertex,SSSVertex>
describeHerwigSSHSFSFVertex("Herwig::SSHSFSFVertex", "HwSusy.so");

void SSHSFSFVertex::Init() {

  static ClassDocumentation<SSHSFSFVertex> documentation
    ("The coupling of an MSSM Higgs to a pair of sfermions.");

}

void SSHSFSFVertex::setCoupling(Energy2 q2, tcPDPtr part1, 
				tcPDPtr part2, tcPDPtr part3) {
  // extract particle ids
  long higgs(part1->id()), sq1(part2->id()), sq2(part3->id());
  // higgs first
  if(abs(sq1)<100) swap(higgs,sq1);
  if(abs(sq2)<100) swap(higgs,sq2);
  // squark second, antisquark third
  if(sq1<0) swap(sq1,sq2);
  assert( higgs == 25 || higgs == 35 || 
	  higgs == 36 || abs(higgs) == 37 );
  sq2 *=-1; 
  assert(sq1>0&&sq2>0);
  
  // running coupling
  if( q2 != theq2Last || thegLast==0.) {
    thegLast = weakCoupling(q2);
    theq2Last = q2;
  }
  
  if( higgs == theHLast && sq1 == theSF1Last && sq2 == theSF2Last) {
    norm(thegLast*theCoupLast*UnitRemoval::InvE);
    return;
  }
  theHLast = higgs;
  theSF1Last = sq1;
  theSF2Last = sq2;
  if( abs(higgs) == ParticleID::Hplus ) 
    chargedHiggs(q2,sq1, sq2);
  else {
    long sm(0);
    unsigned int alpha(sq1/1000000 - 1), beta(sq2/1000000 - 1);
    if( sq1/1000000 == 2 )
      sm = sq1 - 2000000;
    else
      sm = sq1 - 1000000;
    if( sm < 7  ) {
      if( sm % 2 == 0 ) 
	upSF(q2,higgs, sm, alpha, beta);
      else 
	downSF(q2,higgs, sm, alpha, beta);
    }
    else 
      leptonSF(q2,higgs, sm, alpha, beta);
  
  }
 
  norm(thegLast*theCoupLast*UnitRemoval::InvE);
}

void SSHSFSFVertex::downSF(Energy2 q2, long higgs, long sm, 
			   unsigned int alpha, unsigned int beta) {
   assert( sm < 7 && sm % 2 != 0);
   Energy fmass = theMSSM->mass(q2,getParticleData(sm));
   double mfacta = 0.5*fmass/theMw;
   
   if( higgs == ParticleID::A0 ) {
     theCoupLast = -Complex(0.,1.)*mfacta*(theTriC[sm - 1]*theTanB + theMu);
     if(sm == 5) {
      // do not revert to *=, breaks with XCode 5.1
       theCoupLast = theCoupLast * (
	 (*theMix[1])(alpha, 1) * (*theMix[1])(beta , 0) -
	 (*theMix[1])(alpha, 0) * (*theMix[1])(beta , 1) );
     }
     else if(alpha<beta) theCoupLast *= -1.;
     return;
   }
   Energy mfactb = sqr(fmass)/theMw/theCosB;
   Energy facta = theMz/theCw;
   double factb = theMSSM->ed()*theSw*theSw;
   //mixing parameters
   Complex q1a(0.), q1b(0.), q2a(0.), q2b(0.);  
   if( sm == 1 || sm == 3) {
     q1a = (alpha == 0) ? 1.0 : 0.0;
     q1b = (beta == 0) ? 1.0 : 0.0;
     q2a = (alpha == 0) ? 0.0 : 1.0;
     q2b = (beta == 0) ? 0.0 : 1.0;
   }
   else {
     q1a = (*theMix[1])(alpha, 0);
     q1b = (*theMix[1])(beta, 0);
     q2a = (*theMix[1])(alpha, 1);
     q2b = (*theMix[1])(beta, 1);
   }
   Complex fbrac = (q1a*q1b*(0.5 + factb) - factb*q2a*q2b);
   Complex sbrac = (q1a*q1b + q2a*q2b);
   Complex tbrac = (q2a*q1b + q1a*q2b);
   if( higgs == ParticleID::h0) {
     theCoupLast = -facta*theSinAB*fbrac + mfactb*theSinA*sbrac
       + mfacta*(theTriC[sm - 1]*theSinA + theMu*theCosA)*tbrac/theCosB;
   }
   else if( higgs == ParticleID::H0 ) {
     theCoupLast = facta*theCosAB*fbrac - mfactb*theCosA*sbrac 
       + mfacta*(theMu*theSinA - theTriC[sm - 1]*theCosA)*tbrac/theCosB;
   }
   else
     throw HelicityConsistencyError() 
      << "SSHSFSFVertex::downSF - Unrecognised higgs particle " 
      << higgs << Exception::warning;
   
 }

void SSHSFSFVertex::upSF(Energy2 q2, long higgs, long sm, 
			 unsigned int alpha, unsigned int beta) {
  assert( sm < 7 && sm % 2 == 0);
  Energy fmass = theMSSM->mass(q2,getParticleData(sm));
  double mfacta = 0.5*fmass/theMw;
  if( higgs == ParticleID::A0 ){
    theCoupLast = -Complex(0.,1.)*mfacta*(theTriC[sm - 1]/theTanB + theMu);
    if(sm == 6) {
    // do not revert to *=, breaks with XCode 5.1
      theCoupLast = theCoupLast * (
	(*theMix[0])(alpha, 1) * (*theMix[0])(beta , 0) -
	(*theMix[0])(alpha, 0) * (*theMix[0])(beta , 1) );
    }
    else if(alpha<beta) theCoupLast *= -1.;
    return;
  }
  Energy mfactb = sqr(fmass)/theMw/theSinB;
  Energy facta = theMz/theCw;
  double factb = theMSSM->eu()*sqr(theSw);
  //mixing parameters
  Complex q1a(0.), q1b(0.), q2a(0.), q2b(0.);  
  if( sm == 2 || sm == 4) {
    q1a = (alpha == 0) ? 1.0 : 0.0;
    q1b = (beta  == 0) ? 1.0 : 0.0;
    q2a = (alpha == 0) ? 0.0 : 1.0;
    q2b = (beta  == 0) ? 0.0 : 1.0;
  }
  else {
    q1a = (*theMix[0])(alpha, 0);
    q1b = (*theMix[0])(beta , 0);
    q2a = (*theMix[0])(alpha, 1);
    q2b = (*theMix[0])(beta , 1);
  }
  Complex fbrac = (q1a*q1b*(0.5 - factb) + factb*q2a*q2b);
  Complex sbrac = (q1a*q1b + q2a*q2b);
  Complex tbrac = (q2a*q1b + q1a*q2b);
  if( higgs == ParticleID::h0) {
    theCoupLast = facta*theSinAB*fbrac - mfactb*theCosA*sbrac
      - mfacta*(theTriC[sm - 1]*theCosA + theMu*theSinA)*tbrac/theSinB;

  }
  else if( higgs == ParticleID::H0 ) {
    theCoupLast = -facta*theCosAB*fbrac - mfactb*theSinA*sbrac 
      - mfacta*(theTriC[sm - 1]*theSinA - theMu*theCosA )*tbrac/theSinB;

  }
  else
    assert(false);
}

void SSHSFSFVertex::leptonSF(Energy2 q2, long higgs, long sm, 
			     unsigned int alpha, unsigned int beta) {
  assert( sm >= 11 && sm <= 16 ); 
  Energy facta = theMz/theCw;
  if( sm % 2 == 0 ) {
    theCoupLast = complex<Energy>(facta/2.);
    if( higgs == ParticleID::h0) 
      theCoupLast *=  theSinAB;
    else
      theCoupLast *= -theCosAB;
    return;
  }
  Energy fmass = theMSSM->mass(q2,getParticleData(sm));
  double mfacta = fmass/2./theMw;
  if( higgs == ParticleID::A0 ) {
    theCoupLast = -Complex(0.,1.)*mfacta*(theTriC[(sm + 1)/2]*theTanB + theMu);
    if(sm == 15) {
      // do not revert to *=, breaks with XCode 5.1
      theCoupLast = theCoupLast * (
	(*theMix[2])(alpha, 1) * (*theMix[2])(beta , 0) -
	(*theMix[2])(alpha, 0) * (*theMix[2])(beta , 1) );
    }
    else if(alpha<beta) theCoupLast *= -1.;
    return;
  }
  Energy mfactb = fmass*fmass/theMw/theCosB;
  double factb = theSw*theSw;
   //mixing parameters
   Complex l1a(0.), l1b(0.), l2a(0.), l2b(0.);  
   if( sm == 11 || sm == 13) {
     l1a = (alpha == 0) ? 1.0 : 0.0;
     l1b = (beta == 0) ? 1.0 : 0.0;
     l2a = (alpha == 0) ? 0.0 : 1.0;
     l2b = (beta == 0) ? 0.0 : 1.0;
   }
   else {
     l1a = (*theMix[2])(alpha, 0);
     l1b = (*theMix[2])(beta, 0);
     l2a = (*theMix[2])(alpha, 1);
     l2b = (*theMix[2])(beta, 1);
   }
   Complex fbrac = (l1a*l1b*(0.5 - factb) + factb*l2a*l2b);
   Complex sbrac = (l1a*l1b + l2a*l2b);
   Complex tbrac = (l2a*l1b + l1a*l2b);
   if( higgs == ParticleID::h0) {
     theCoupLast = -facta*theSinAB*fbrac + mfactb*theSinA*sbrac
       + mfacta*(theTriC[(sm + 1)/2]*theSinA + theMu*theCosA)*tbrac/theCosB;
   }
   else if( higgs == ParticleID::H0 ) {
     theCoupLast = facta*theCosAB*fbrac - mfactb*theCosA*sbrac 
       + mfacta*(theMu*theSinA - theTriC[(sm + 1)/2]*theCosA)*tbrac/theCosB;
   }
   else
     throw HelicityConsistencyError() 
      << "SSHSFSFVertex::leptonSF - Unrecognised higgs particle " 
      << higgs << Exception::warning;
  
  
}

void SSHSFSFVertex::chargedHiggs(Energy2 q2, long id1, long id2) {
  //have id1 as up-type
  if( id1 % 2 != 0)
    swap(id1, id2);
  unsigned int beta(0);
  if( id2/1000000 == 2 )
    beta = 1;
  else 
    beta = 0;

  long smd = (beta == 0) ? id2 - 1000000 : id2 - 2000000;
  Energy mfd = theMSSM->mass(q2,getParticleData(smd));
  Energy2 facta = sqr(theMw)*2.*theSinB*theCosB;
  if( smd == 11 || smd == 13 || smd == 15) {
    Complex l1b = (beta == 0) ? 1.0 : 0.0;
    Complex l2b = (beta == 0) ? 0.0 : 1.0;
    if( smd == 15 ) {
      l1b = (*theMix[2])(beta, 0);
      l2b = (*theMix[2])(beta, 1);
    }
    theCoupLast = ( l1b*(mfd*mfd*theTanB - facta) 
		    + l2b*mfd*(theTriC[(smd + 1)/2]*theTanB + theMu)
		   )/theMw/sqrt(2.);
  }
  else {
    unsigned int alpha(0);
    if( id1/1000000 == 2 )
    alpha = 1;
  else 
    alpha = 0;

  long smu = (alpha == 0) ? id1 - 1000000 : id1 - 2000000;
  Energy mfu = theMSSM->mass(q2,getParticleData(smu));
  Complex q1a(0.0), q1b(0.0), q2a(0.0), q2b(0.0);
  if( smu == 2 || smu == 4 ) {
    q1a = (alpha == 0) ? 1.0 : 0.0;
    q2a = (alpha == 0) ? 0.0 : 1.0;
  }
  else {
    q1a = (*theMix[0])(alpha,0);
    q2a = (*theMix[0])(alpha,1);
  }
  if( smd == 1 || smd == 3 ) {
    q1b = (beta == 0) ? 1.0 : 0.0;
    q2b = (beta == 0) ? 0.0 : 1.0;
  }
  else {
    q1b = (*theMix[1])(beta,0);
    q2b = (*theMix[1])(beta,1);
  }
  
  theCoupLast = ( q1a*q1b*(mfd*mfd*theTanB + mfu*mfu/theTanB - facta)
		  + q2a*q2b*mfu*mfd*(theTanB + (1./theTanB))
		  + q1a*q2b*mfd*(theTriC[smd - 1]*theTanB + theMu)
		  + q2a*q1b*mfu*(theMu + theTriC[smu-1]/theTanB)
		 )/theMw/sqrt(2.);
  }
}
#line 1 "./SSNFSVertex.cc"
// -*- C++ -*-
//
// SSNFSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSNFSVertex class.
//

#include "SSNFSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSNFSVertex::SSNFSVertex() :  _sw(0.), _cw(0.), _mw(), 
			     _sb(0.), _cb(0.), _q2last(), _couplast(0.),
			     _leftlast(0.), _rightlast(0.), _id1last(0), 
			     _id2last(0),
			      yukawa_(1) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SSNFSVertex::persistentOutput(PersistentOStream & os) const {
  os << _stop << _sbot << _stau << _nmix << _theSS  << _sw << _cw 
     << ounit(_mw,GeV) << _sb << _cb << yukawa_;
}

void SSNFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _stop >> _sbot >> _stau >> _nmix >> _theSS >> _sw >> _cw 
     >> iunit(_mw,GeV) >> _sb >> _cb >> yukawa_;
}

void SSNFSVertex::doinit() {
  long neut[5] = {1000022, 1000023, 1000025, 1000035, 1000045};
  for(unsigned int nl = 0; nl < 5; ++nl) {
    //quarks
    for(long ix=1;ix<7;++ix){
      addToList( neut[nl],  ix, -(1000000+ix) );
      addToList( neut[nl],  ix, -(2000000+ix) );
      addToList( neut[nl], -ix,  (1000000+ix) );
      addToList( neut[nl], -ix,  (2000000+ix) );
    }
    //leptons
    for(long ix=11;ix<17;++ix) {
      addToList( neut[nl],  ix, -(1000000+ix) );
      addToList( neut[nl], -ix,  (1000000+ix) );
     
      if( ix % 2 != 0 ) {
	addToList( neut[nl],  ix, -(2000000+ix) );
	addToList( neut[nl], -ix,  (2000000+ix) );
      }
    }
    
  }
  FFSVertex::doinit();
  _theSS = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  if(!_theSS)
    throw InitException() << "SSGSSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;

  _stop = _theSS->stopMix();
  _sbot = _theSS->sbottomMix();
  _stau = _theSS->stauMix();
  _nmix = _theSS->neutralinoMix();
  if(!_stop || !_stau || !_sbot || !_nmix)
    throw InitException() << "SSNFSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << _stop << " sbottom: "
			  << _sbot << " stau: " << _stau 
			  << " N: " << _nmix << Exception::abortnow;

  _sw = sqrt(sin2ThetaW());
  _mw = getParticleData(24)->mass();
  double tb = _theSS->tanBeta();
  _cw = sqrt(1. - sqr(_sw));
  _sb = tb/sqrt(1 + sqr(tb));
  _cb = sqrt(1 - sqr(_sb));
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSNFSVertex,FFSVertex>
describeHerwigSSNFSVertex("Herwig::SSNFSVertex", "HwSusy.so");

void SSNFSVertex::Init() {

  static ClassDocumentation<SSNFSVertex> documentation
    ("The SSNFSVertex implements the coupling of a neutralino to "
     "a fermion-sfermion");

  static Switch<SSNFSVertex,unsigned int> interfaceYukawa
    ("Yukawa",
     "Whether or not to include the Yukawa type couplings",
     &SSNFSVertex::yukawa_, 1, false, false);
  static SwitchOption interfaceYukawaYes
    (interfaceYukawa,
     "Yes",
     "Include the terms",
     1);
  static SwitchOption interfaceYukawaNo
    (interfaceYukawa,
     "No",
     "Don't include them",
     0);
  static SwitchOption interfaceYukawa3rdGen
    (interfaceYukawa,
     "ThirdGeneration",
     "Only include for the third generation",
     2);
}

void SSNFSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  long isc(abs(part3->id())), ism(abs(part1->id())),
    ineut(abs(part2->id()));
  tcPDPtr smfermion = part1;
  if( ism / 1000000 == 1 )  {
    swap( ism, ineut);
    smfermion = part2;
  }
  
  if(q2!=_q2last || _couplast==0.) {
    _couplast = -sqrt(2)*weakCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);

  if( ineut != _id1last || ism != _id2last || isc != _id3last ) {
    _id1last = ineut;
    _id2last = ism;
    _id3last = isc;
    // determine neutralino and squark eigenstates
    unsigned int alpha(isc/1000000 - 1), nl(0);
    switch( ineut ) {
    case 1000022 : nl = 0;
      break;
    case 1000023 : nl = 1;
      break;
    case 1000025 : nl = 2;
      break;
    case 1000035 : nl = 3;
      break;
    case 1000045 : nl = 4;
      break;
    default : assert(false);
    }
    // common primed neutralino matrices
    Complex n2prime = (*_nmix)(nl,1)*_cw - (*_nmix)(nl,0)*_sw;
    //handle neutrinos first
    if( ism == 12 || ism == 14 || ism == 16 ) {
      _leftlast = Complex(0., 0.);
      _rightlast = n2prime/2./_cw;
    }
    else {
      Complex n1prime = (*_nmix)(nl,0)*_cw + (*_nmix)(nl,1)*_sw;
      tcPDPtr smf = getParticleData(ism);
      double qf = smf->charge()/eplus;
      Complex bracketl = qf*_sw*( conj(n1prime) - _sw*conj(n2prime)/_cw );
      double y = 0.;
      if(yukawa_==1 || ((ism==5 || ism==6 || ism==15) && yukawa_==2))
	y = double(_theSS->mass(q2, smf)/2./_mw);
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
    }
  }
  //determine the helicity order of the vertex
  if( smfermion->id() < 0 ) {
    left(conj(_rightlast));
    right(conj(_leftlast));
  }
  else {
    left(_leftlast);
    right(_rightlast);
  }
}
#line 1 "./SSWSSVertex.cc"
// -*- C++ -*-
//
// SSWSSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWSSVertex class.
//

#include "SSWSSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSWSSVertex::SSWSSVertex():_sw(0.), _cw(0.), _q2last(),_couplast(0.), 
			   _ulast(0), _dlast(0), _gblast(0),
			   _factlast(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SSWSSVertex::doinit() {
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
 
  //LL-sleptons
  for(long ix=1000011;ix<1000016;ix+=2) {
    addToList(-24,-ix,ix+1);
  }
  //2-L stau
  addToList(-24,-2000015,1000016);
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

  //LL-sleptons
  for(long ix=1000011;ix<1000016;ix+=2) {
    addToList(24,ix,-ix-1);
  }
  //2-L stau
  addToList(24,2000015,-1000016);
  
  //---Z0----
//LL-sleptons
  for(long ix=1000011;ix<1000017;++ix) {
    addToList(23,ix,-ix);
  }
  //RR-sleptons
  for(long ix=2000011;ix<2000016;ix+=2) {
    addToList(23,ix,-ix);
  }
  //L-Rbar stau
  addToList(23,1000015,-2000015);
  //Lbar-R stau
  addToList(23,-1000015,2000015);
   
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
  
  //----gamma----
  //sleptons
  for(long ix=1000011;ix<1000016;ix+=2) {
    addToList(22,ix,-ix);
  }
  for(long ix=2000011;ix<2000016;ix+=2) {
    addToList(22,ix,-ix);
  }
  //squarks
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(22,ix,-ix);
  }
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(22,ix,-ix);
  }
  VSSVertex::doinit();
  tMSSMPtr theSS = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  if(!theSS)
    throw InitException() << "SSWSSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt( 1. - sqr(_sw) );
  _stop = theSS->stopMix();
  _sbottom = theSS->sbottomMix();
  _stau = theSS->stauMix();
  if(!_stop || !_stau || !_sbottom)
    throw InitException() << "SSWSSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << _stop << " sbottom: " << _sbottom
			  << " stau: " << _stau << Exception::abortnow;
}

void SSWSSVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw  << _cw << _stau << _stop << _sbottom;
}

void SSWSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _stau >> _stop >> _sbottom;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSWSSVertex,VSSVertex>
describeHerwigSSWSSVertex("Herwig::SSWSSVertex", "HwSusy.so");

void SSWSSVertex::Init() {

  static ClassDocumentation<SSWSSVertex> documentation
    ("This is the implementation of the coupling of an SM boson "
     "a pair of sfermions");
  
}

void SSWSSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  long boson(abs(part1->id()));
  assert( boson == ParticleID::Wplus || boson == ParticleID::Z0 ||
	  boson == ParticleID::gamma );
  long sf1(abs(part2->id())),sf2(abs(part3->id()));

  assert( (sf1 >= 1000001 && sf1 <= 1000006) 
	  || (sf1 >= 1000011 && sf1 <= 1000016)
	  || (sf1 >= 2000001 && sf1 <= 2000006)
	  || (sf1 >= 2000011 && sf1 <= 2000016) );
  
  assert( (sf2 >= 1000001 && sf2 <= 1000006) 
	  || (sf2 >= 1000011 && sf2 <= 1000016)
	  || (sf2 >= 2000001 && sf2 <= 2000006)
	  || (sf2 >= 2000011 && sf2 <= 2000016) );

  if( sf1 % 2 != 0 ) swap(sf1, sf2);
  if( sf1 != _ulast || sf2 != _dlast || boson != _gblast) {
    _gblast = boson;
    _ulast = sf1;
    _dlast = sf2;
    //photon is simplest
    if( boson == ParticleID::gamma )
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
      if( boson == ParticleID::Wplus ) {
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
  if( q2 != _q2last || _couplast==0. ) {
    _q2last = q2;
    _couplast = electroMagneticCoupling(q2);
  }
  if(part2->id()>0)
    norm(-_couplast*_factlast);
  else
    norm(+_couplast*_factlast);
}

#line 1 "./SSWWSSVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWWSSVertex class.
//

#include "SSWWSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SSWWSSVertex::SSWWSSVertex() : sw_(0.), cw_(0.), q2last_(), couplast_(0.), 
			       ulast_(0), dlast_(0), gblast1_(0), gblast2_(0),
			       factlast_(0.) {
  colourStructure(ColourStructure::DELTA);
}

IBPtr SSWWSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SSWWSSVertex::fullclone() const {
  return new_ptr(*this);
}

void SSWWSSVertex::persistentOutput(PersistentOStream & os) const {
  os << sw_  << cw_ << stau_ << stop_ << sbottom_;
}

void SSWWSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> sw_ >> cw_ >> stau_ >> stop_ >> sbottom_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSWWSSVertex,VVSSVertex>
describeHerwigSSWWSSVertex("Herwig::SSWWSSVertex", "HwSusy.so");

void SSWWSSVertex::Init() {

  static ClassDocumentation<SSWWSSVertex> documentation
    ("The SSWWSSVertex class implements the coupling of two"
     " weak bosons and two scalar fermions");

}

void SSWWSSVertex::doinit() {
  // gamma gamma, gamma Z0 and Z0 Z0 PAIRS
  for(long ib1=22;ib1<24;++ib1) {
    for(long ib2=22;ib2<24;++ib2) {
      // sleptons
      long istep = ib1==23&&ib2==23 ? 1 : 2;
      for(long ix=1000011;ix<1000017;ix+=istep) {
	addToList(ib1,ib2,ix,-ix);
      }
      for(long ix=2000011;ix<2000017;ix+=2) {
	addToList(ib1,ib2,ix,-ix);
      }
      // squarks
      for(long ix=1000001;ix<1000007;++ix) {
	addToList(ib1,ib2,ix,-ix);
      }
      for(long ix=2000001;ix<2000007;++ix) {
	addToList(ib1,ib2,ix,-ix);
      }
      // L/R mixing if Z
      if(ib2==23) {
	//L-Rbar stau
	addToList(ib1,ib2,1000015,-2000015);
	//Lbar-R stau
	addToList(ib1,ib2,-1000015,2000015);
	//L-Rbar stop
	addToList(ib1,ib2,1000006,-2000006);
	//Lbar-R stop
	addToList(ib1,ib2,-1000006,2000006);
	//L-Rbar sbottom
	addToList(ib1,ib2,1000005,-2000005);
	//Lbar-R sbottom
	addToList(ib1,ib2,-1000005,2000005);
      }
    }
    // W- gamma/Z0
    // LL-squarks
    for(long ix=1000001;ix<1000006;ix+=2) {
      addToList(-24,ib1,ix+1,-ix);
    }
    // 1-2 stop sbottom
    addToList(-24,ib1,1000006,-2000005);
    // 2-1 stop sbottom
    addToList(-24,ib1,2000006,-1000005);
    // 2-2 stop sbottom
    addToList(-24,ib1,2000006,-2000005);
    // LL-sleptons
    for(long ix=1000011;ix<1000016;ix+=2) {
      addToList(-24,ib1,-ix,ix+1);
    }
    //2-L stau
    addToList(-24,ib1,-2000015,1000016);
    // W+ gamma/Z0
    // LL-squarks
    for(long ix=1000001;ix<1000006;ix+=2) {
      addToList(24,ib1,-(ix+1),ix);
    }
    // 1-2 stop sbottom
    addToList(24,ib1,-1000006,2000005);
    // 2-1 stop sbottom
    addToList(24,ib1,-2000006,1000005);
    // 2-2 stop sbottom
    addToList(24,ib1,-2000006,2000005);
    // LL-sleptons
    for(long ix=1000011;ix<1000016;ix+=2) {
      addToList(24,ib1,ix,-ix-1);
    }
    //2-L stau
    addToList(24,ib1,2000015,-1000016);
  }
  // finally WW pairs
  // sleptons
  for(long ix=1000011;ix<=1000017;++ix) {
    addToList(24,-24,ix,-ix);
  }
  addToList(24,-24,1000015,-2000015);
  addToList(24,-24,2000015,-1000015);
  // squarks
  for(long ix=1000001;ix<=1000007;++ix) {
    addToList(24,-24,ix,-ix);
  }
  addToList(24,-24,1000005,-2000005);
  addToList(24,-24,2000005,-1000005);
  addToList(24,-24,1000006,-2000006);
  addToList(24,-24,2000006,-1000006);
  // couplings etc
  orderInGem(2);
  orderInGs(0);
  VVSSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "SSWWSSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  sw_ = sqrt(sin2ThetaW());
  cw_ = sqrt( 1. - sqr(sw_) );
  stop_ = model->stopMix();
  sbottom_ = model->sbottomMix();
  stau_ = model->stauMix();
  if(!stop_ || !stau_ || !sbottom_)
    throw InitException() << "SSWWSSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << stop_ << " sbottom: " << sbottom_
			  << " stau: " << stau_ << Exception::abortnow;
}

void SSWWSSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			       tcPDPtr part2,tcPDPtr part3,
			       tcPDPtr part4) {
  long boson[2] = {abs(part1->id()),abs(part2->id())};
  assert( boson[0] == ParticleID::Wplus || boson[0] == ParticleID::Z0 ||
	  boson[0] == ParticleID::gamma );
  assert( boson[1] == ParticleID::Wplus || boson[1] == ParticleID::Z0 ||
	  boson[1] == ParticleID::gamma );
  long sf1(abs(part3->id())),sf2(abs(part4->id()));

  assert( (sf1 >= 1000001 && sf1 <= 1000006) 
  	  || (sf1 >= 1000011 && sf1 <= 1000016)
  	  || (sf1 >= 2000001 && sf1 <= 2000006)
  	  || (sf1 >= 2000011 && sf1 <= 2000016) );
  
  assert( (sf2 >= 1000001 && sf2 <= 1000006) 
  	  || (sf2 >= 1000011 && sf2 <= 1000016)
  	  || (sf2 >= 2000001 && sf2 <= 2000006)
  	  || (sf2 >= 2000011 && sf2 <= 2000016) );

  if( sf1 % 2 != 0 ) swap(sf1, sf2);
  if(boson[0]>boson[1]) swap(boson[0],boson[1]);
  if( sf1 != ulast_ || sf2 != dlast_ ||
      boson[0] != gblast1_ || boson[1] != gblast2_) {
    gblast1_ = boson[0];
    gblast2_ = boson[1];
    ulast_ = sf1;
    dlast_ = sf2;
    double ef1 = getParticleData(sf1)->charge()/eplus;
    double ef2 = getParticleData(sf2)->charge()/eplus;
    // photon pairs are simplest
    if(gblast1_ == ParticleID::gamma &&
       gblast2_ == ParticleID::gamma ) {
      factlast_ = 2.*sqr( ef1 );
    }
    // otherwise we need the helicity states
    else {
      //determine which helicity state
      unsigned int alpha(sf1/1000000 - 1), beta(sf2/1000000 - 1);
      //mixing factors
      Complex m1a(0.), m1b(0.);
      if( sf1 == ParticleID::SUSY_t_1 || sf1 == ParticleID::SUSY_t_2 )
  	m1a = (*stop_)(alpha, 0);
      else if( sf1 == ParticleID::SUSY_b_1 || sf1 == ParticleID::SUSY_b_2 )
  	m1a = (*sbottom_)(alpha, 0);
      else if( sf1 == ParticleID::SUSY_tau_1minus || 
  	       sf1 == ParticleID::SUSY_tau_2minus )
  	m1a = (*stau_)(alpha, 0);
      else
  	m1a = (alpha == 0) ? Complex(1.) : Complex(0.);
      
      if( sf2 == ParticleID::SUSY_t_1 || sf2 == ParticleID::SUSY_t_2 )
  	m1b = (*stop_)(beta, 0);
      else if( sf2 == ParticleID::SUSY_b_1 || sf2 == ParticleID::SUSY_b_2 )
  	m1b = (*sbottom_)(beta, 0);
      else if( sf2 == ParticleID::SUSY_tau_1minus || 
  	       sf2 == ParticleID::SUSY_tau_2minus )
  	m1b = (*stau_)(beta, 0);
      else
  	m1b = (beta == 0) ? Complex(1.) : Complex(0.);
      // if either boson is a W
      if(gblast2_==ParticleID::Wplus) {
	// WW
	if(gblast1_==ParticleID::Wplus) {
	  factlast_ = 0.5*m1a*m1b/sqr(sw_);
	}
	// gamma W
	else if(gblast1_==ParticleID::gamma) {
	  factlast_ = sqrt(0.5)*m1a*m1b/sw_*(ef1+ef2);
	}
	// Z0 W
	else if(gblast1_==ParticleID::Z0) {
	  factlast_ = -sqrt(0.5)*m1a*m1b/cw_*(ef1+ef2);
	}
      }
      else {
	// compute the Z coupling
	factlast_=1.;
	if(gblast1_==ParticleID::Z0||gblast2_==ParticleID::Z0) {
	  if( sf1 == ParticleID::SUSY_nu_eL || sf1 == ParticleID::SUSY_nu_muL ||
	      sf1 == ParticleID::SUSY_nu_tauL ) {
	    factlast_ = 0.5/cw_/sw_;
	  }
	  else {
	    double lmda =  sf2 % 2 == 0 ? 1. : -1.;
	    factlast_ = lmda*m1a*m1b;
	    if( alpha == beta) factlast_ -= 2.*ef1*sqr(sw_);
	    factlast_ *= 0.5/cw_/sw_; 
	  }
	}
	// photon Z
	if(gblast1_ == ParticleID::gamma &&
	   gblast2_ == ParticleID::Z0 ) {
	  factlast_ *= 2.*ef1; 
	}
	// Z pairs
	else if(gblast1_ == ParticleID::Z0 &&
		gblast2_ == ParticleID::Z0 ) {
	  factlast_ *= 2.*factlast_;
	}
      }
    }
  }
  if( q2 != q2last_ || couplast_==0. ) {
    q2last_ = q2;
    couplast_ = sqr(electroMagneticCoupling(q2));
  }
  norm(couplast_*factlast_);
}

#line 1 "./SSWGSSVertex.cc"
// -*- C++ -*-
//
// SSWGSSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWGSSVertex class.
//

#include "SSWGSSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSWGSSVertex::SSWGSSVertex() : _sw(0.), _cw(0.), _q2last(),_emcouplast(0.),
			       _scouplast(0.), _ulast(0), _dlast(0),
			       _gblast(0), _factlast(0.)  {
  colourStructure(ColourStructure::SU3TFUND);
}

void SSWGSSVertex::doinit() {
  //W-
  //LL-squarks
  for(long ix=1000001;ix<1000006;ix+=2) {
    addToList(-24,21,ix+1,-ix);
  }
  //1-2 stop sbottom
  addToList(-24,21,1000006,-2000005);
  //2-1 stop sbottom
  addToList(-24,21,2000006,-1000005);
  //2-2 stop sbottom
  addToList(-24,21,2000006,-2000005);

  //W+
  for(long ix=1000001;ix<1000006;ix+=2) {
    addToList(24,21,-(ix+1),ix);
  }
  //1-2 stop sbottom
  addToList(24,21,-1000006,2000005);
  //2-1 stop sbottom
  addToList(24,21,-2000006,1000005);
  //2-2 stop sbottom
  addToList(24,21,-2000006,2000005);
  
  //---Z0----   
  //LL squarks
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(23,21,ix,-ix);
  }
  //RR squarks
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(23,21,ix,-ix);
  }
  //L-Rbar stop
  addToList(23,21,1000006,-2000006);
  //Lbar-R stop
  addToList(23,21,-1000006,2000006);

  //L-Rbar sbottom
  addToList(23,21,1000005,-2000005);
  //Lbar-R sbottom
  addToList(23,21,-1000005,2000005);
  
  //----gamma----
  //squarks
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(22,21,ix,-ix);
  }
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(22,21,ix,-ix);
  }
  orderInGem(1);
  orderInGs(1);

  VVSSVertex::doinit();
  tMSSMPtr theSS = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  if(!theSS)
    throw InitException() << "SSWGSSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt( 1. - sqr(_sw) );
  _stop = theSS->stopMix();
  _sbottom = theSS->sbottomMix();
  if(!_stop || !_sbottom)
    throw InitException() << "SSWGSSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << _stop << " sbottom: " << _sbottom
			  << Exception::abortnow;
}


void SSWGSSVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw  << _cw << _stop << _sbottom;
}

void SSWGSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _stop >> _sbottom;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSWGSSVertex,VVSSVertex>
describeHerwigSSWGSSVertex("Herwig::SSWGSSVertex", "HwSusy.so");


void SSWGSSVertex::Init() {

  static ClassDocumentation<SSWGSSVertex> documentation
    ("This implements the gluon-gluon-squark-squark vertex.");
}

void SSWGSSVertex::setCoupling(Energy2 q2,   tcPDPtr part1,
			       tcPDPtr part2,tcPDPtr part3,tcPDPtr part4) { 

  long boson(abs(part1->id()));
  long gluon(abs(part2->id()));
  if (gluon > boson) swap(gluon, boson);

  if( boson != ParticleID::Wplus && boson != ParticleID::Z0 && 
      boson != ParticleID::gamma ) {
    throw HelicityConsistencyError()
      << "SSWGSSVertex::setCoupling() - Vector particle in this "
      << "vertex is not a W/Z/gamma. " << boson << Exception::warning;
    norm(0.);
  }
  
  if( gluon != ParticleID::g ) {
    throw HelicityConsistencyError()
      << "SSWGSSVertex::setCoupling() - Vector particle in this "
      << "vertex is not a gluon. " << gluon << Exception::warning;
    norm(0.);
  }
  long sq1(abs(part3->id())),sq2(abs(part4->id()));

  if( sq1 % 2 != 0 ) swap(sq1, sq2);
  if( sq1 != _ulast || sq2 != _dlast || boson != _gblast) {
    _gblast = boson;
    _ulast = sq1;
    _dlast = sq2;
    //photon is simplest
    if( boson == ParticleID::gamma )
      _factlast = -2.*getParticleData(sq1)->charge()/eplus;

    else {
     //determine which helicity state
      unsigned int alpha(sq1/1000000 - 1), beta(sq2/1000000 - 1);
      //mixing factors
      Complex m1a(0.), m1b(0.);
      if( sq1 == ParticleID::SUSY_t_1 || sq1 == ParticleID::SUSY_t_2 )
	m1a = (*_stop)(alpha, 0);
      else if( sq1 == ParticleID::SUSY_b_1 || sq1 == ParticleID::SUSY_b_2 ) 
	m1a = (*_sbottom)(alpha, 0);
      else
	m1a = (alpha == 0) ? Complex(1.) : Complex(0.);

      if( sq2 == ParticleID::SUSY_t_1 || sq2 == ParticleID::SUSY_t_2 )
	m1b = (*_stop)(beta, 0);
      else if( sq2 == ParticleID::SUSY_b_1 || sq2 == ParticleID::SUSY_b_2 )
	m1b = (*_sbottom)(beta, 0);
      else
	m1b = (beta == 0) ? Complex(1.) : Complex(0.);
    
      //W boson
      if( boson == ParticleID::Wplus ) {
	_factlast = -1.*m1a*m1b*sqrt(2)/_sw;
      }
      //Z boson
      else {
	double lmda(1.);
	if( sq2 % 2 == 0 ) lmda = -1.;
	_factlast = lmda*m1a*m1b;
	if( alpha == beta) {
	    double ef = getParticleData(sq1)->charge()/eplus;
	    _factlast += 2.*ef*sqr(_sw);
	}
	_factlast *= 1./_cw/_sw; 	
      }
    }
  }
  if( q2 != _q2last || _emcouplast==0. || _scouplast==0. ) {
    _q2last = q2;
    _emcouplast = electroMagneticCoupling(q2);
    _scouplast = strongCoupling(q2);
  }
  norm(-_emcouplast*_scouplast*_factlast);
}
#line 1 "./SSGSSVertex.cc"
// -*- C++ -*-
//
// SSGSSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGSSVertex class.
//

#include "SSGSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGSSVertex::SSGSSVertex() : _couplast(0.),_q2last(ZERO) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TFUND);
}

void SSGSSVertex::doinit() {
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(21,ix,-ix);
  }
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(21,ix,-ix);
  }
  VSSVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SSGSSVertex,Helicity::VSSVertex>
describeSSGSSVertex("Herwig::SSGSSVertex", "HwSusy.so");

void SSGSSVertex::Init() {

  static ClassDocumentation<SSGSSVertex> documentation
    ("The SSGSSVertex class implements the coupling"
     " of the gluon to the squarks");

}

void SSGSSVertex::setCoupling(Energy2 q2,
#ifndef NDEBUG
			      tcPDPtr part1,
#else
			      tcPDPtr,
#endif
			      tcPDPtr part2,
#ifndef NDEBUG
				tcPDPtr part3) {
#else
				tcPDPtr ) {
#endif
  assert(part1->id()==ParticleID::g);
#ifndef NDEBUG
  long isf = abs(part2->id());
#endif
  assert( (isf >= 1000001 && isf <= 1000006) || 
	  (isf >= 2000001 && isf <= 2000006) );
  assert(part2->id()==-part3->id());
  if(q2 != _q2last || _couplast == 0.) {
    _couplast = strongCoupling(q2);
    _q2last = q2;
  }
  if(part2->id()>0)
    norm(-_couplast);
  else
    norm( _couplast);
}
#line 1 "./SSGGSQSQVertex.cc"
// -*- C++ -*-
//
// SSGGSQSQVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGGSQSQVertex class.
//

#include "SSGGSQSQVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGGSQSQVertex::SSGGSQSQVertex() : q2last_(),couplast_(0.) {
  colourStructure(ColourStructure::SU3TTFUNDS);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SSGGSQSQVertex,Helicity::VVSSVertex>
describeSSGGSQSQVertex("Herwig::SSGGSQSQVertex", "HwSusy.so");

void SSGGSQSQVertex::Init() {

  static ClassDocumentation<SSGGSQSQVertex> documentation
    ("This implements the gluon-gluon-squark-squark vertex.");

}

void SSGGSQSQVertex::setCoupling(Energy2 q2, tcPDPtr, tcPDPtr, tcPDPtr,
				 tcPDPtr) { 
  if(q2 != q2last_ || couplast_ == 0.) {
    couplast_ = sqr(strongCoupling(q2));
    q2last_ = q2;
  }
  norm(couplast_);
}

void SSGGSQSQVertex::doinit() {
  //L-L squarks
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(21,21,ix,-ix);
  }
  //R-R squarks
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(21,21,ix,-ix);
  }
  orderInGs(2);
  orderInGem(0);
  VVSSVertex::doinit();
}
#line 1 "./SSGSGSGVertex.cc"
// -*- C++ -*-
//
// SSGSGSGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGSGSGVertex class.
//

#include "SSGSGSGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGSGSGVertex::SSGSGSGVertex() : _couplast(0.),_q2last(ZERO) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3F);
}

// Static variable needed for the type description system in ThePEG.
DescribeNoPIOClass<SSGSGSGVertex,FFVVertex>
describeHerwigSSGSGSGVertex("Herwig::SSGSGSGVertex", "HwSusy.so");

void SSGSGSGVertex::Init() {

  static ClassDocumentation<SSGSGSGVertex> documentation
    ("This class implements the gluon-gluino-gluino vertex");

}

void SSGSGSGVertex::setCoupling(Energy2 q2,
#ifndef NDEBUG
				tcPDPtr part1,
#else
				tcPDPtr ,
#endif
#ifndef NDEBUG
				tcPDPtr part2,
#else
				tcPDPtr ,
#endif
#ifndef NDEBUG
				tcPDPtr part3) {
#else
				tcPDPtr ) {
#endif
  assert(part1->id()==ParticleID::SUSY_g &&
	 part2->id()==ParticleID::SUSY_g &&
	 part3->id() == ParticleID::g);
  if(q2 != _q2last || _couplast==0.) {
    _couplast = strongCoupling(q2);
    _q2last = q2;
  }
  norm(_couplast);
  left(1.);
  right(1.);
}

void SSGSGSGVertex::doinit() {
  addToList(1000021, 1000021, 21);
  FFVVertex::doinit();
}
#line 1 "./SSNNZVertex.cc"
// -*- C++ -*-
//
// SSNNZVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSNNZVertex class.
//

#include "SSNNZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "MixingMatrix.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSNNZVertex::SSNNZVertex() : _sw(0.), _cw(0.), _id1last(0), 
			     _id2last(0), _q2last(), _couplast(0.),
			     _leftlast(0.), _rightlast(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void SSNNZVertex::doinit() {
  long neu[] = { 1000022, 1000023, 1000025, 1000035, 1000045 };
  for(unsigned int i = 0; i < 5; ++i)
    for(unsigned int j = 0; j < 5; ++j)
      addToList(neu[i], neu[j], 23);
  FFVVertex::doinit();
  tSusyBasePtr theSS = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!theSS)
    throw InitException() << "SSNNZVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  
  _theN  = theSS->neutralinoMix();
  if(!_theN)
    throw InitException() << "SSNNZVertex::doinit - The neutralino "
			  << "mixing matrix pointer is null." 
			  << Exception::abortnow;
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt(1 - _sw*_sw);
}

void SSNNZVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _cw << _theN;
}

void SSNNZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _theN;
  _id1last = 0;
  _id2last = 0;
  _q2last = ZERO;
  _couplast = 0.;
  _leftlast = 0.;
  _rightlast = 0.;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSNNZVertex,Helicity::FFVVertex>
describeSSNNZVertex("Herwig::SSNNZVertex", "HwSusy.so");

void SSNNZVertex::Init() {

  static ClassDocumentation<SSNNZVertex> documentation
    ("The coupling of a Z-boson to a pair of neutralinos");

}

void SSNNZVertex::setCoupling(Energy2 q2,tcPDPtr part1,
#ifndef NDEBUG
			      tcPDPtr part2,tcPDPtr part3) {
#else
			      tcPDPtr part2,tcPDPtr) {
#endif
  assert(part3->id() == ParticleID::Z0);
  long ic1 = part2->id();
  long ic2 = part1->id();
  assert(ic1 == ParticleID::SUSY_chi_10 || ic1 == ParticleID::SUSY_chi_20 ||
	 ic1 == ParticleID::SUSY_chi_30 || ic1 == ParticleID::SUSY_chi_40 ||
	 ic1 == 1000045);
  assert(ic2 == ParticleID::SUSY_chi_10 || ic2 == ParticleID::SUSY_chi_20 ||
	 ic2 == ParticleID::SUSY_chi_30 || ic2 == ParticleID::SUSY_chi_40 ||
	 ic2 == 1000045);
  if(q2 != _q2last || _couplast==0.) {
    _q2last = q2;
    _couplast = weakCoupling(q2)/_cw;
  }
  if(ic1 != _id1last || ic2 != _id2last) {
    _id1last = ic1;
    _id2last = ic2;
    unsigned int neu1(ic1 - 1000022), neu2(ic2 - 1000022);
    if(neu1 > 1) {
      if(ic1 == 1000025)
	neu1 = 2;
      else if(ic1 == 1000035)
	neu1 = 3;
      else 
	neu1 = 4;
    }
    if(neu2 > 1) {
      if(ic2 == 1000025)
	neu2 = 2;
      else if(ic2 == 1000035)
	neu2 = 3;
      else
	neu2 = 4;
    }
    _leftlast = 0.5*( (*_theN)(neu1, 3)*conj((*_theN)(neu2, 3)) -
		      (*_theN)(neu1, 2)*conj((*_theN)(neu2, 2)) );
    _rightlast = -conj(_leftlast);
  }
  norm(_couplast);
  left(_leftlast);
  right(_rightlast);
}
#line 1 "./SSCCZVertex.cc"
// -*- C++ -*-
//
// SSCCZVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSCCZVertex class.
//

#include "SSCCZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSCCZVertex::SSCCZVertex() : _sw2(0.), _cw(0.), _couplast(0.),
			     _q2last(), _id1last(0), _id2last(0),
			     _leftlast(0.), _rightlast(0.), _gblast(0) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void SSCCZVertex::doinit() {
  addToList(-1000024, 1000024, 23);
  addToList(-1000024, 1000037, 23);

  addToList(-1000037, 1000024, 23);
  addToList(-1000037, 1000037, 23);

  //photon
  addToList(-1000024, 1000024, 22);
  addToList(-1000037, 1000037, 22);
  FFVVertex::doinit();
  tSusyBasePtr theSS = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!theSS) 
    throw InitException() << "SSCCZVertex::doinit - The model pointer "
				     << "is null! "
				     << Exception::abortnow;
  _sw2 = sin2ThetaW();
  _cw = sqrt(1. - _sw2);
  _theU = theSS->charginoUMix();
  _theV = theSS->charginoVMix();
  if(!_theU || !_theV)
    throw InitException() << "SSCCZVertex::doinit - "
			  << "A mixing matrix pointer is null.  U: " 
			  << _theU << "  V: " << _theV
			  << Exception::abortnow;
}

void SSCCZVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw2 << _cw << _theU << _theV;
}

void SSCCZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw2 >> _cw >> _theU >> _theV;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSCCZVertex,Helicity::FFVVertex>
describeSSCCZVertex("Herwig::SSCCZVertex", "HwSusy.so");

void SSCCZVertex::Init() {

  static ClassDocumentation<SSCCZVertex> documentation
    ("This class implements the coupling of a Z/gamma to a pair of"
     " charginos. ");

}

void SSCCZVertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
			      tcPDPtr part3) {
  long ichar1(part1->id()), ichar2(part2->id()), boson(part3->id());
  assert( boson == ParticleID::gamma || boson == ParticleID::Z0);
  assert( abs(ichar1) == 1000024 || abs(ichar1) == 1000037); 
  assert( abs(ichar2) == 1000024 || abs(ichar2) == 1000037); 
  if(_q2last != q2||_couplast==0.) {
    _q2last = q2;
    _couplast = electroMagneticCoupling(q2);
  }
  norm(_couplast);
  if(boson != _gblast || ichar1 != _id1last || ichar2 != _id2last) {
    _gblast = boson;
    _id1last = ichar1;
    _id2last = ichar2;
    if( boson == ParticleID::Z0 ) {
      unsigned int ic1(0), ic2(0);
      if(abs(ichar1) == 1000037) ic1 = 1;
      if(abs(ichar2) == 1000037) ic2 = 1;
      _leftlast = -(*_theV)(ic1, 0)*conj((*_theV)(ic2, 0)) - 
	0.5*(*_theV)(ic1, 1)*conj((*_theV)(ic2, 1));
      _rightlast = -conj((*_theU)(ic1, 0))*(*_theU)(ic2, 0) - 
	0.5*conj((*_theU)(ic1, 1))*(*_theU)(ic2, 1);
      if(abs(ichar1) == abs(ichar2)) {
	_leftlast += _sw2;
	_rightlast += _sw2;
      }
      _leftlast /= sqrt(_sw2)*_cw;
      _rightlast /= sqrt(_sw2)*_cw;
    }
    else {
      if(abs(ichar1) == abs(ichar2)) {
	_leftlast  = -1.;
	_rightlast = -1.;
      }
      else {
	_leftlast  = 0.;
	_rightlast = 0.;
      }
    }
    if(ichar1>0) {
      Complex temp = _leftlast;
      _leftlast  = -_rightlast;
      _rightlast = -temp;
    }
  }
  left(_leftlast);
  right(_rightlast);
}
#line 1 "./SSCNWVertex.cc"
// -*- C++ -*-
//
// SSCNWVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSCNWVertex class.
//

#include "SSCNWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSCNWVertex::SSCNWVertex() : _sw(0.),  _couplast(0.), _q2last(ZERO), 
			     _id1last(0), _id2last(0), _leftlast(0.),
			     _rightlast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void SSCNWVertex::doinit() {
  long neu[] = { 1000022, 1000023, 1000025, 1000035, 1000045 };
  long cha[] = { 1000024, 1000037 };
  // sign == -1 outgoing W-, sign == +1 outgoing W+
  for(int sign = -1; sign < 2; sign += 2)
    for(unsigned int ine = 0; ine < 5; ++ine)
      for(unsigned int ic = 0; ic < 2; ++ic)
	addToList(-sign*cha[ic], neu[ine], sign*24);
  FFVVertex::doinit();
  tSusyBasePtr theSS = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!theSS)
    throw InitException() << "SSCNWVertex::doinit() - The model pointer is null!"
			  << Exception::abortnow;
  _sw = sqrt(sin2ThetaW());
  
  _theN = theSS->neutralinoMix();
  _theU = theSS->charginoUMix();
  _theV = theSS->charginoVMix();

  if(!_theN || !_theU || ! _theV)
    throw InitException() << "SSCNWVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " N: " << _theN << " U: " << _theU << " V: "
			  << _theV << Exception::abortnow;
}

void SSCNWVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _theN << _theU << _theV;
}

void SSCNWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _theN >> _theU >> _theV;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSCNWVertex,Helicity::FFVVertex>
describeSSCNWVertex("Herwig::SSCNWVertex", "HwSusy.so");

void SSCNWVertex::Init() {

  static ClassDocumentation<SSCNWVertex> documentation
    ("This class implements the coupling of a W boson to a "
     "neutralino and a chargino");

}

void SSCNWVertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
#ifndef NDEBUG
			      tcPDPtr part3) {
#else
			      tcPDPtr) {
#endif
  assert(abs(part3->id()) == ParticleID::Wplus);
  long neu, cha;
  if(part1->charged()) {
    cha = part1->id();
    neu = part2->id();
  }
  else {
    cha = part2->id();
    neu = part1->id();
  }
  assert((abs(cha) == 1000024 || abs(cha) == 1000037) && 
    (neu == 1000022 || neu == 1000023 || 
     neu == 1000025 || neu == 1000035 || 
     neu == 1000045) );
  if(q2 != _q2last||_couplast==0.) {
    _q2last = q2;
    _couplast = weakCoupling(q2);
  }
  norm(_couplast);
  if(cha != _id1last || neu != _id2last) {
    _id1last = cha;
    _id2last = neu;
    unsigned int eigc = abs(cha) == 1000037 ? 1 : 0;
    unsigned int eign(0);
    if     (neu == 1000023) eign = 1;
    else if(neu == 1000025) eign = 2;
    else if(neu == 1000035) eign = 3;
    else if(neu == 1000045) eign = 4;
    _leftlast = (*_theN)(eign, 1)*conj((*_theV)(eigc, 0)) - 
      ( (*_theN)(eign, 3)*conj((*_theV)(eigc, 1))/sqrt(2));
    _rightlast = conj((*_theN)(eign, 1))*(*_theU)(eigc, 0) +
      ( conj((*_theN)(eign, 2))*(*_theU)(eigc, 1)/sqrt(2));
  }
  Complex ltemp = _leftlast;
  Complex rtemp = _rightlast;
  // conjugate if +ve chargino
  if(cha>0) {
    ltemp = conj(ltemp);
    rtemp = conj(rtemp);
  }
  if((part1->id()==cha&&cha>0)||(part2->id()==cha&&cha<0)) {
    Complex temp = ltemp;
    ltemp  = -rtemp;
    rtemp = -temp;
  }
  left (ltemp);
  right(rtemp);
}
#line 1 "./SSFFHVertex.cc"
// -*- C++ -*-
//
// SSFFHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSFFHVertex class.
//

#include "SSFFHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSFFHVertex::SSFFHVertex() : thetanb(0.0), theMw(ZERO), 
			     theSa(0.0), theSb(0.0),
			     theCa(0.0), theCb(0.0),
			     theFLast(make_pair(0,0)), theGlast(0.),
			     theq2last(), theMassLast(make_pair(ZERO,ZERO)) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SSFFHVertex::doinit() {
  int higgs[] = { 25, 35, 36 };
  for ( long h = 0; h < 3; ++h ) {
    //neutral higgs
    // quarks
    for(long ix=1;ix<7;++ix)
      addToList(-ix,ix,higgs[h]);
    // charged leptons
    for(long ix=11;ix<16;ix+=2)
      addToList(-ix,ix,higgs[h]);
  }
  for(long ix=1;ix<6;ix+=2) {
    //outgoing H+
    addToList(-ix-1,  ix, 37);
    //outgoing H-
    addToList(-ix  ,ix+1,-37);
  }
  for(long ix=11;ix<16;ix+=2) {
    //outgoing H+
    addToList(-ix-1,  ix, 37);
    //outgoing H-
    addToList(-ix  ,ix+1,-37);
  }
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSFFHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  
  theMw = getParticleData(ParticleID::Wplus)->mass();
  thetanb = theMSSM->tanBeta();
  theSb = thetanb/sqrt(1. + sqr(thetanb));
  theCb = sqrt( 1. - sqr(theSb) );
  theSa = sin(theMSSM->higgsMixingAngle());
  theCa = sqrt(1. - sqr(theSa));
  
  FFSVertex::doinit();
}

void SSFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << theMSSM << thetanb << ounit(theMw,GeV) << theSa
     << theSb << theCa << theCb;
}

void SSFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMSSM >> thetanb >> iunit(theMw,GeV) >> theSa
     >> theSb >> theCa >> theCb;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSFFHVertex,FFSVertex>
describeHerwigSSFFHVertex("Herwig::SSFFHVertex", "HwSusy.so");

void SSFFHVertex::Init() {

  static ClassDocumentation<SSFFHVertex> documentation
    ("The coupling of the higgs bosons to SM fermions in the MSSM");

}

void SSFFHVertex::setCoupling(Energy2 q2, tcPDPtr particle1,
			      tcPDPtr particle2,tcPDPtr particle3) {
  long f1ID(abs(particle1->id())), f2ID(abs(particle2->id())), 
    higgsID(particle3->id());
  // check higgs
  assert( higgsID == ParticleID::h0 ||     higgsID  == ParticleID::H0 || 
	  higgsID == ParticleID::A0 || abs(higgsID) == ParticleID::Hplus );
  // check fermions
  assert(!( ((f1ID > 6 && f1ID < 11) || f1ID > 16 ) ||
	    ((f2ID > 6 && f1ID < 11) || f2ID > 16 ) ));
  if( q2 != theq2last || theGlast==0.) {
    theGlast = weakCoupling(q2);
  }
  if( q2 != theq2last || theFLast.first  != f1ID) {
    theMassLast.first  = theMSSM->mass(q2,particle1);
    theFLast.first  = f1ID;
  }
  if( q2 != theq2last || theFLast.second != f2ID) {
    theMassLast.second = theMSSM->mass(q2,particle2);
    theFLast.second = f2ID;
  }
  theq2last = q2;
  Complex coup(0.);
  Complex lcoup(1.),rcoup(1.);
  if( higgsID == ParticleID::h0 || higgsID == ParticleID::H0 ||
      higgsID == ParticleID::A0 ) {
    coup = 0.5*theMassLast.first/theMw;
    if( higgsID == ParticleID::h0 ) {
      if( f1ID % 2 == 0 )
	coup *= -theCa/theSb;
      else
	coup *=  theSa/theCb;
    }
    else if( higgsID == ParticleID::H0 ) {
      if( f1ID % 2 == 0 )
	coup *= -theSa/theSb;
      else
	coup *= -theCa/theCb;
    }
    else {
      if( f1ID % 2 == 0 )
	coup /= thetanb; 
      else
	coup *= thetanb;
      coup *= Complex(0.,-1.);
      rcoup = -1.;
    }
  }
  //H+
  else {
    if( f1ID % 2 == 0 ) {
      lcoup = theMassLast.first /thetanb/theMw;
      rcoup = theMassLast.second*thetanb/theMw;
    }
    else {
      lcoup = theMassLast.second/thetanb/theMw;
      rcoup = theMassLast.first *thetanb/theMw;
    }
    coup = sqrt(0.5);
    if( higgsID > 0 ) swap(lcoup,rcoup);
  }
  norm(theGlast*coup);
  left (lcoup);
  right(rcoup);
}

#line 1 "./SSGOGOHVertex.cc"
// -*- C++ -*-
//
// SSGOGOHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGOGOHVertex class.
//

#include "SSGOGOHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGOGOHVertex::SSGOGOHVertex() : theMw(), theSij(2, {0., 0.}),
				 theQij(2, {0., 0.}),
				 theQijLp(4, {0., 0.}),
				 theQijRp(4, {0., 0.}),
				 theSijdp(4, {0., 0., 0., 0.}),
				 theQijdp(4, {0., 0., 0., 0.}),
				 theSa(0.0), theSb(0.0),
				 theCa(0.0), theCb(0.0), theCoupLast(0.0),
				 theLLast(0.0), theRLast(0.0), theHLast(0),
				 theID1Last(0), theID2Last(0), theq2last() {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void SSGOGOHVertex::doinit() {
  long neu[4] = {1000022, 1000023, 1000025, 1000035};
  long chg[2] = {1000024, 1000037};
  long higgs[3] =  {25, 35, 36};
  for(unsigned int i = 0; i < 4; ++i) {
    //neutralinos
    for(unsigned int j = 0; j < 4; ++j) {
      for(unsigned int k = 0; k < 4; ++k) {
	if( i < 3 ) {
	  if(k<=j)
	    addToList(neu[j], neu[k], higgs[i]);
	  //charginos
	  if( j < 2 && k < 2 ) {
	    addToList(-chg[j], chg[k], higgs[i]);
	  }
	} 
	else {
	  if( k == 2 ) break;
	  addToList(-chg[k], neu[j], 37);
	  addToList( chg[k], neu[j],-37);
	}
      }
    }
  }
  FFSVertex::doinit();
  
  tMSSMPtr theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSGOGOHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  
  theMw = getParticleData(ParticleID::Wplus)->mass();
  double theSw = sqrt(sin2ThetaW());
  double tw = theSw/sqrt(1. - theSw*theSw);
  double tanb = theMSSM->tanBeta();
  theSb = tanb/sqrt(1. + sqr(tanb));
  theCb = sqrt( 1. - sqr(theSb) );
  theSa = sin(theMSSM->higgsMixingAngle());
  theCa = sqrt(1. - sqr(theSa));
  MixingMatrix nmix = *theMSSM->neutralinoMix();
  MixingMatrix umix = *theMSSM->charginoUMix();
  MixingMatrix vmix = *theMSSM->charginoVMix();

  for(unsigned int i = 0; i < 4; ++i) {
    for(unsigned int j = 0; j < 4; ++j) {
      if( i < 2 && j < 2 ) { 
	theQij[i][j] = vmix(i,0)*umix(j,1)/sqrt(2);
	theSij[i][j] = vmix(i,1)*umix(j,0)/sqrt(2);
      }
      if( j < 2 ) {
	theQijLp[i][j] = conj(nmix(i, 3)*vmix(j,0) 
			      + (nmix(i,1) + nmix(i,0)*tw)*vmix(j,1)/sqrt(2));
	theQijRp[i][j] = nmix(i, 2)*umix(j,0) 
	  - (nmix(i,1) + nmix(i,0)*tw)*umix(j,1)/sqrt(2);
      }
      theQijdp[i][j] = 0.5*( nmix(i,2)*( nmix(j,1) - tw*nmix(j,0) )
			     + nmix(j,2)*( nmix(i,1) - tw*nmix(i,0) ) );
      theSijdp[i][j] = 0.5*( nmix(i,3)*( nmix(j,1) - tw*nmix(j,0) )
			     + nmix(j,3)*( nmix(i,1) - tw*nmix(i,0) ) );
    }
  }
}

void SSGOGOHVertex::persistentOutput(PersistentOStream & os) const {
  os << theSij << theQij << theQijLp << theQijRp << theSijdp
     << theQijdp << ounit(theMw,GeV) << theSa << theSb << theCa 
     << theCb;
}

void SSGOGOHVertex::persistentInput(PersistentIStream & is, int) {
  is >> theSij >> theQij >> theQijLp >> theQijRp >> theSijdp
     >> theQijdp >> iunit(theMw,GeV) >> theSa >> theSb >> theCa 
     >> theCb;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSGOGOHVertex,FFSVertex>
describeHerwigSSGOGOHVertex("Herwig::SSGOGOHVertex", "HwSusy.so");

void SSGOGOHVertex::Init() {

  static ClassDocumentation<SSGOGOHVertex> documentation
    ("The coupling of the higgs bosons to SM fermions in the MSSM");

}

/// \todo fixme
void SSGOGOHVertex::setCoupling(Energy2 q2, tcPDPtr particle1, 
				tcPDPtr particle2,tcPDPtr particle3) {
  long f1ID(particle1->id()), f2ID(particle2->id()), higgsID(particle3->id());
  assert(higgsID == ParticleID::h0 ||     higgsID  == ParticleID::H0 ||
	 higgsID == ParticleID::A0 || abs(higgsID) == ParticleID::Hplus);
  if( f1ID < 0 ) swap(f1ID, f2ID);
  
  if( q2 != theq2last || theCoupLast == 0. ) {
    theCoupLast = weakCoupling(q2);
    theq2last = q2;
  }
  if( higgsID == theHLast && f1ID == theID1Last && f2ID == theID2Last) {
    norm(theCoupLast);
    left(theLLast);
    right(theRLast);
    return;
  }
  theHLast = higgsID;
  theID1Last = f1ID;
  theID2Last = f2ID;
  
  if( higgsID == ParticleID::h0 ) {
    //charginos
    if( abs(f2ID) == ParticleID::SUSY_chi_1plus ||
	abs(f2ID) == ParticleID::SUSY_chi_2plus ) {
      unsigned int ei = (abs(f1ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      unsigned int ej = (abs(f2ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      theLLast =  conj(theQij[ej][ei])*theSa - conj(theSij[ej][ei])*theCa;
      theRLast = theQij[ei][ej]*theSa - theSij[ei][ej]*theCa;
    }
    //neutralinos
    else {
      unsigned int ei(f1ID - ParticleID::SUSY_chi_10), 
	ej(f2ID - ParticleID::SUSY_chi_10);
      if( ei > 1 )
	ei = ( ei == 13 ) ? 3 : 2;
      if( ej > 1 )
	ej = ( ej == 13 ) ? 3 : 2;
      theLLast = conj(theQijdp[ej][ei])*theSa + conj(theSijdp[ej][ei])*theCa;
      theRLast = theQijdp[ei][ej]*theSa + theSijdp[ei][ej]*theCa ;
    }
    
  }
  else if( higgsID == ParticleID::H0 ) {
    //charginos
    if( abs(f2ID) == ParticleID::SUSY_chi_1plus ||
	abs(f2ID) == ParticleID::SUSY_chi_2plus ) {
      unsigned int ei = (abs(f1ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      unsigned int ej = (abs(f2ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      theLLast =  -conj(theQij[ej][ei])*theCa - conj(theSij[ej][ei])*theSa;
      theRLast = -theQij[ei][ej]*theCa - theSij[ei][ej]*theSa;
    }
    //neutralinos
    else {
      unsigned int ei(f1ID - ParticleID::SUSY_chi_10), 
	ej(f2ID - ParticleID::SUSY_chi_10);
      if( ei > 1 )
	ei = ( ei == 13 ) ? 3 : 2;
      if( ej > 1 )
	ej = ( ej == 13 ) ? 3 : 2;
      
      theLLast = -conj(theQijdp[ej][ei])*theCa + conj(theSijdp[ej][ei])*theSa;
      theRLast = -theQijdp[ei][ej]*theCa + theSijdp[ei][ej]*theSa;
    }
  }
  else if( higgsID == ParticleID::A0 ) {
    if( abs(f2ID) == ParticleID::SUSY_chi_1plus ||
	abs(f2ID) == ParticleID::SUSY_chi_2plus ) {
      unsigned int ei = (abs(f1ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      unsigned int ej = (abs(f2ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;

      theLLast = Complex(0.,1.)*( conj(theQij[ej][ei])*theSb 
				  + conj(theSij[ej][ei])*theCb );
      theRLast = -Complex(0.,1.)*(theQij[ei][ej]*theSb + theSij[ei][ej]*theCb);
    }
    //neutralinos
    else {
      unsigned int ei(f1ID - ParticleID::SUSY_chi_10), 
	ej(f2ID - ParticleID::SUSY_chi_10);
      if( ei > 1 )
	ei = ( ei == 13 ) ? 3 : 2;
      if( ej > 1 )
	ej = ( ej == 13 ) ? 3 : 2;

      theLLast = Complex(0.,1.)*( conj(theQijdp[ej][ei])*theSb 
				  - conj(theSijdp[ej][ei])*theCb );
      theRLast = -Complex(0.,1.)*(theQijdp[ei][ej]*theSb - theSijdp[ei][ej]*theCb);
    }
  }
  //charged higgs
  else {
    unsigned int ei(0),ej(0);
    long chg(f2ID), neu(f1ID);
    if( abs(neu) == ParticleID::SUSY_chi_1plus || 
	abs(neu) == ParticleID::SUSY_chi_2plus ) swap(chg, neu);
    ej = ( abs(chg) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
    ei = neu - ParticleID::SUSY_chi_10;
    if( ei > 1 ) ei = ( ei == 13 ) ? 3 : 2;
    theLLast = -theQijLp[ei][ej]*theCb;
    theRLast = -theQijRp[ei][ej]*theSb;
    if( higgsID < 0 ) {
      Complex tmp = theLLast;
      theLLast = conj(theRLast);
      theRLast = conj(tmp);
    }
  }
  norm(theCoupLast);
  left(theLLast);
  right(theRLast);
}

#line 1 "./SSWWHVertex.cc"
// -*- C++ -*-
//
// SSWWHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWWHVertex class.
//

#include "SSWWHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSWWHVertex::SSWWHVertex() : theh0Wfact(ZERO), theH0Wfact(ZERO), 
			     theh0Zfact(ZERO), theH0Zfact(ZERO),
			     theCoupLast(ZERO), theElast(0.0),
			     theq2last(ZERO), theHlast(0), 
			     theGBlast(0) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void SSWWHVertex::doinit() {
  addToList(23,23,25);
  addToList(-24,24,25);
  addToList(23,23,35);
  addToList(-24,24,35);
  VVSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !model )
    throw InitException() 
      << "SSWWHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;

  Energy mw = getParticleData(ParticleID::Wplus)->mass();
  Energy mz = getParticleData(ParticleID::Z0)->mass();
  double sw = sqrt(sin2ThetaW());
  double sinalp = sin(model->higgsMixingAngle());
  double cosalp = sqrt(1. - sqr(sinalp));
  double tanbeta = model->tanBeta();
  double sinbeta = tanbeta/sqrt(1. + sqr(tanbeta));
  double cosbeta = sqrt( 1. - sqr(sinbeta) );
  double sinbma = sinbeta*cosalp - cosbeta*sinalp;
  double cosbma = cosbeta*cosalp + sinbeta*sinalp;
  
  theh0Wfact = mw*sinbma/sw;
  theH0Wfact = mw*cosbma/sw;
  theh0Zfact = mz*sinbma/sw/sqrt(1. - sw*sw);
  theH0Zfact = mz*cosbma/sw/sqrt(1. - sw*sw);
}

void SSWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theh0Wfact,GeV) << ounit(theH0Wfact,GeV) 
     << ounit(theh0Zfact,GeV) << ounit(theH0Zfact,GeV);
}

void SSWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theh0Wfact,GeV) >> iunit(theH0Wfact,GeV) 
     >> iunit(theh0Zfact,GeV) >> iunit(theH0Zfact,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSWWHVertex,VVSVertex>
describeHerwigSSWWHVertex("Herwig::SSWWHVertex", "HwSusy.so");

void SSWWHVertex::Init() {

  static ClassDocumentation<SSWWHVertex> documentation
    ("This is the coupling of a pair of SM gauge bosons"
     "to the higgs particles in the MSSM");

}

void SSWWHVertex::setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr,
			      tcPDPtr particle3) {
  long bosonID = abs(particle1->id());
  long higgsID =     particle3->id();
  assert( higgsID == ParticleID::h0    || higgsID == ParticleID::H0 );
  assert( bosonID == ParticleID::Wplus || bosonID == ParticleID::Z0 );
  if( higgsID != theHlast || bosonID != theGBlast) {
    if( higgsID == ParticleID::h0 )
      theCoupLast = (bosonID == ParticleID::Z0) ? theh0Zfact : theh0Wfact;
    else
      theCoupLast = (bosonID == ParticleID::Z0) ? theH0Zfact : theH0Wfact;
  }
  if( q2 != theq2last ) {
    theq2last = q2;
    theElast = electroMagneticCoupling(q2);
  }
  norm(theElast*theCoupLast*UnitRemoval::InvE);
}
#line 1 "./SSWWHHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWWHHVertex class.
//

#include "SSWWHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "MSSM.h"

using namespace Herwig;

SSWWHHVertex::SSWWHHVertex()  : couplast_(0.), q2last_(ZERO) {
  orderInGs(0);
  orderInGem(2);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr SSWWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SSWWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void SSWWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << coup_;
}

void SSWWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> coup_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSWWHHVertex,VVSSVertex>
describeHerwigSSWWHHVertex("Herwig::SSWWHHVertex", "HwSusy.so");

void SSWWHHVertex::Init() {

  static ClassDocumentation<SSWWHHVertex> documentation
    ("The SSWWHHVertex class implements the coupling of two Higgs bosons and"
     "two electroweak vector bosons in the MSSM.");

}

void SSWWHHVertex::doinit() {
  int id[3]={25,35,36};
  for(unsigned int ix=0;ix<3;++ix) {
    addToList( 24,-24,id[ix],id[ix]);
    addToList( 23, 23,id[ix],id[ix]);
    addToList( 22, 24,id[ix],-37);
    addToList( 22,-24,id[ix], 37);
    addToList( 23, 24,id[ix],-37);
    addToList( 23,-24,id[ix], 37);
  }
  addToList( 24,-24, 37,-37);
  addToList( 23, 23, 37,-37);
  addToList( 22, 23, 37,-37);
  addToList( 22, 22, 37,-37);
  VVSSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !model )
    throw Exception() 
      << "SSWWHHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  coup_.resize(11);
  double sw2 = sin2ThetaW();
  double sw  = sqrt(sw2);
  double cw2 = 1.-sw2;
  double cw  = sqrt(cw2);
  double c2w = cw2-sw2;
  double sinalp = sin(model->higgsMixingAngle());
  double cosalp = sqrt(1. - sqr(sinalp));
  double tanbeta = model->tanBeta();
  double sinbeta = tanbeta/sqrt(1. + sqr(tanbeta));
  double cosbeta = sqrt( 1. - sqr(sinbeta) );
  double sinbma = sinbeta*cosalp - cosbeta*sinalp;
  double cosbma = cosbeta*cosalp + sinbeta*sinalp;
  // WWHH
  coup_[0] = 0.5/sw2;
  // ZZH0H0
  coup_[1] = 0.5/sw2/cw2;
  // ZZH+H-
  coup_[2] = 0.5*sqr(c2w)/cw2/sw2;
  // Z W h0 H+
  coup_[3] =-0.5/cw*cosbma;
  // Z W H0 H+
  coup_[4] = 0.5/cw*sinbma;
  // Z W A0 H+
  coup_[5] =-Complex(0.,0.5)/cw;
  // A A H+H-
  coup_[6] = 2.;
  // A Z H+H-
  coup_[7] = c2w/sw/cw;
  // A W h0 H+
  coup_[8] = 0.5*cosbma/sw;
  // A W H0 H+
  coup_[9] =-0.5*sinbma/sw;
  // A W A0 H+
  coup_[10] = Complex(0.,0.5)/sw;
}

void SSWWHHVertex::setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
			       tcPDPtr part3, tcPDPtr part4) {
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = sqr(electroMagneticCoupling(q2));
    q2last_=q2;
  }
  int ibos1 = part1->id(), ibos2 = part2->id();
  int isca1 = part3->id(), isca2 = part4->id();
  if(abs(ibos1)==abs(ibos2)) {
    if(abs(ibos1)==ParticleID::Wplus) {
      norm(couplast_*coup_[0]);
    }
    else if(ibos1==ParticleID::Z0) {
      if(abs(isca1)==ParticleID::Hplus)
	norm(couplast_*coup_[2]);
      else
	norm(couplast_*coup_[1]);
    }
    else if(ibos1==ParticleID::gamma) {
	norm(couplast_*coup_[6]);
    }
    else
      assert(false);
  }
  else if(abs(ibos1)==ParticleID::Wplus ||
	  abs(ibos2)==ParticleID::Wplus) {
    if(abs(ibos1)==ParticleID::gamma ||
       abs(ibos2)==ParticleID::gamma) {
      if(abs(isca1)==ParticleID::h0 ||
	 abs(isca2)==ParticleID::h0) {
	norm(couplast_*coup_[8]);
      }
      else if(abs(isca1)==ParticleID::H0 ||
	      abs(isca2)==ParticleID::H0) {
	norm(couplast_*coup_[9]);
      }
      else if(abs(isca1)==ParticleID::A0 ||
	      abs(isca2)==ParticleID::A0) {
	if(isca1==ParticleID::Hplus ||
	   isca2==ParticleID::Hplus) {
	  norm(couplast_*     coup_[10] );
	}
	else {
	  norm(couplast_*conj(coup_[10]));
	}
      }
      else
	assert(false);
    }
    else {
      if(abs(isca1)==ParticleID::h0 ||
	 abs(isca2)==ParticleID::h0) {
	norm(couplast_*coup_[3]);
      }
      else if(abs(isca1)==ParticleID::H0 ||
	      abs(isca2)==ParticleID::H0) {
	norm(couplast_*coup_[4]);
      }
      else if(abs(isca1)==ParticleID::A0 ||
	      abs(isca2)==ParticleID::A0) {
	if(isca1==ParticleID::Hplus ||
	   isca2==ParticleID::Hplus) {
	  norm(couplast_*     coup_[5] );
	}
	else {
	  norm(couplast_*conj(coup_[5]));
	}
      }
      else
	assert(false);
    }
  }
  else {
    norm(couplast_*coup_[7]);
  }
}
#line 1 "./SSWHHVertex.cc"
// -*- C++ -*-
//
// SSWHHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWHHVertex class.
//

#include "SSWHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSWHHVertex::SSWHHVertex() : 
  theSw(0.), theS2w(0.), theC2w(0.), thesbma(0.), thecbma(0.), 
  theq2last(ZERO), theElast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void SSWHHVertex::doinit() {
  addToList(22,37,-37);
  addToList(23,36,25);
  addToList(23,36,35);
  addToList(23,37,-37);
  //outgoing W+
  addToList(24,-37,25);
  addToList(24,-37,35);
  addToList(24,-37,36);
  //outgoing W-
  addToList(-24,37,25);
  addToList(-24,37,35);
  addToList(-24,37,36);
  VSSVertex::doinit();
  tMSSMPtr theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSWHHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  theSw = sqrt(sin2ThetaW());
  double cw = sqrt(1. - sqr(theSw));
  theS2w = 2.*theSw*cw;
  theC2w = cw*cw - theSw*theSw;

  double sina = sin(theMSSM->higgsMixingAngle());
  double cosa =  sqrt(1. - sqr(sina));
  double tanb = theMSSM->tanBeta();
  double sinb = tanb/sqrt(1. + sqr(tanb));
  double cosb = sqrt( 1. - sqr(sinb) );
  thesbma = sinb*cosa - sina*cosb;
  thecbma = cosa*cosb + sina*sinb;
}

void SSWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << theSw << theS2w << theC2w << thesbma << thecbma;
}

void SSWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> theSw >> theS2w >> theC2w >> thesbma >> thecbma;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSWHHVertex,VSSVertex>
describeHerwigSSWHHVertex("Herwig::SSWHHVertex", "HwSusy.so");

void SSWHHVertex::Init() {

  static ClassDocumentation<SSWHHVertex> documentation
    ("The coupling of a pair of higgs bosons and a SM gauge boson");

}

void SSWHHVertex::setCoupling(Energy2 q2, tcPDPtr particle1,
			      tcPDPtr particle2, tcPDPtr particle3) {
  long gboson = particle1->id();
  assert(     gboson  == ParticleID::Z0    ||
	      gboson  == ParticleID::gamma || 
	  abs(gboson) == ParticleID::Wplus );
  long h1ID = particle2->id();
  long h2ID = particle3->id();
  Complex coup(0.);
  if( gboson == ParticleID::Z0 ) {
    if( abs(h1ID) == ParticleID::Hplus ) {
      coup = theC2w/theS2w;
      if(h1ID<0) coup *= -1.;
    }
    else if( h1ID == ParticleID::h0 ||  
	     h2ID == ParticleID::h0 ) {
      coup = Complex(0., 1.)*thecbma/theS2w;
    }
    else {
      coup =-Complex(0., 1.)*thesbma/theS2w;
    }
    if(h2ID==ParticleID::A0) coup *=-1.;
  }
  else if( gboson == ParticleID::gamma ) {
    coup = 1.;
    if(h1ID<0) coup *= -1.;
  }
  else {
    long higgs = abs(h1ID) == ParticleID::Hplus ? h2ID : h1ID;
    if( higgs == ParticleID::h0 ) {
      coup =  0.5*thecbma/theSw;
    }
    else if( higgs == ParticleID::H0) 
      coup = -0.5*thesbma/theSw;
    else 
      coup = -Complex(0., 0.5)/theSw;
    if(abs(h2ID) == ParticleID::Hplus ) coup *= -1.;
    if(gboson<0&&higgs!=ParticleID::A0) coup *= -1.;
  }
  if( q2 != theq2last || theElast==0.) {
    theq2last = q2;
    theElast = electroMagneticCoupling(q2);
  }
  norm(theElast*coup);
}
#line 1 "./SSHHHVertex.cc"
// -*- C++ -*-
//
// SSHHHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSHHHVertex class.
//

#include "SSHHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSHHHVertex::SSHHHVertex() : theMw(ZERO), theZfact(ZERO), theSw(0.),
			     theSbpa(0.), theCbpa(0.), theSbma(0.),
			     theCbma(0.), theS2a(0.), theC2a(0.),
			     theS2b(0.), theC2b(0.), theElast(0.),
			     theq2last(ZERO) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void SSHHHVertex::doinit() {
  long sec = 35;
  for(long h = 25; h < 36; h += 10) {
    //self-coupling
    addToList(h, h, h);
    //first-second
    addToList(h,sec,sec);
    //pseudo-scalar
    addToList(h, 36, 36);
    //charged higgs
    addToList(h, 37,-37);
    
    sec = 25;
  }
  SSSVertex::doinit();
  tMSSMPtr theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSHHHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  
  theMw = getParticleData(ParticleID::Wplus)->mass();
  theSw = sqrt(sin2ThetaW());
  theZfact = getParticleData(ParticleID::Z0)->mass()/2./
    theSw/sqrt(1. - sqr(theSw));
  
  double tanbeta = theMSSM->tanBeta();
  double sinbeta = tanbeta/sqrt(1. + sqr(tanbeta));
  double cosbeta = sqrt(1. - sqr(sinbeta));
  double sinalpha = sin(theMSSM->higgsMixingAngle());
  double cosalpha = sqrt( 1. - sqr(sinalpha) );
  
  theS2a = 2.*sinalpha*cosalpha;
  theS2b = 2.*sinbeta*cosbeta;
  theC2a = cosalpha*cosalpha - sinalpha*sinalpha;
  theC2b = cosbeta*cosbeta - sinbeta*sinbeta;
  theSbpa = sinbeta*cosalpha + sinalpha*cosbeta;
  theCbpa = cosbeta*cosalpha - sinbeta*sinalpha;
  theSbma = sinbeta*cosalpha - sinalpha*cosbeta;
  theCbma = cosbeta*cosalpha + sinbeta*sinalpha;

}

void SSHHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMw,GeV) << ounit(theZfact,GeV) << theSw 
     << theSbpa << theCbpa << theSbma << theCbma << theS2a << theC2a 
     << theS2b << theC2b; 
}

void SSHHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMw,GeV) >> iunit(theZfact,GeV) >> theSw 
     >> theSbpa >> theCbpa >> theSbma >> theCbma >> theS2a >> theC2a 
     >> theS2b >> theC2b;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSHHHVertex,SSSVertex>
describeHerwigSSHHHVertex("Herwig::SSHHHVertex", "HwSusy.so");

void SSHHHVertex::Init() {

  static ClassDocumentation<SSHHHVertex> documentation
    ("This is the coupling of a higgs to a pair of higgs bosons "
     "in the MSSM.");

}

void SSHHHVertex::setCoupling(Energy2 q2, tcPDPtr particle1,
			      tcPDPtr particle2,tcPDPtr particle3) {
  long ids[3] = { abs(particle1->id()), abs(particle2->id()),
		  abs(particle3->id()) };
  long h1(0), h2(0), h3(0), hc(0);
  for(unsigned int i = 0; i < 3; ++i) {
    if( ids[i] == ParticleID::h0) ++h1;
    else if( ids[i] == ParticleID::H0) ++h2;
    else if( ids[i] == ParticleID::A0) ++h3;
    else if( ids[i] == ParticleID::Hplus) ++hc;
    else assert(false);
  }
  assert(h1 + h2 + h3 + hc == 3);
  
  complex<Energy> coupling;
  if( h1 == 3 || h2 == 3 ) {
    coupling = -3.*theZfact*theC2a;
    if( h1 == 3 )
      coupling *= theSbpa;
    else
      coupling *= theCbpa;
  }
  else if( h1 == 1 ) {
    if( h2 == 2 )
      coupling = theZfact*( 2.*theS2a*theCbpa + theSbpa*theC2a );
    else if( h3 == 2 )
      coupling = -theZfact*theC2b*theSbpa;
    else if( hc == 2 )
      coupling = -theMw*theSbma/theSw - theZfact*theC2b*theSbpa;
    else assert(false);
  }
  else if( h2 == 1 ) {
    if( h1 == 2 )
      coupling = -theZfact*( 2.*theS2a*theSbpa - theCbpa*theC2a );
    else if( h3 == 2 )
      coupling = theZfact*theC2b*theCbpa;
    else if( hc == 2 )
      coupling = -theMw*theCbma/theSw + theZfact*theC2b*theCbpa;
    else assert(false);
  }
  
  if( q2 != theq2last || theElast==0. ) {
    theq2last = q2;
    theElast = electroMagneticCoupling(q2);
  }
  norm(theElast*coupling*UnitRemoval::InvE);
}
			      
#line 1 "./SSHGGVertex.cc"
// -*- C++ -*-
//
// SSHGGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSHGGVertex class.
//

#include "SSHGGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Looptools/clooptools.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSHGGVertex::SSHGGVertex() : theIncludeTriLinear(true),
			     thePseudoScalarTreatment(false),
			     theSw(0.), theMw(), theZfact(),
			     theQt1L(0.), theQt1R(0.), theQt1LR(0.), 
			     theQt2L(0.), theQt2R(0.), theQt2LR(0.),
			     theQb1L(0.), theQb1R(0.), theQb1LR(0.), 
			     theQb2L(0.), theQb2R(0.), theQb2LR(0.),
			     theSqmass(4,ZERO),
			     theTanB(0.),theSinA(0.), 
			     theCosA(0.), theSinB(0.), theCosB(0.), 
			     theSinApB(0.), theCosApB(0.), theCouplast(0.), 
			     theq2last(), theHaveCoeff(false), theLastID(0) {
  orderInGs(2);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

void SSHGGVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(21,21,25);
  addToList(21,21,35);
  addToList(21,21,36);
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM ) 
    throw InitException()
      << "SSHGGVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  theMw = getParticleData(ParticleID::Wplus)->mass();
  thetop = getParticleData(ParticleID::t);
  thebot = getParticleData(ParticleID::b);
  theSw = sqrt(sin2ThetaW());
  theZfact = getParticleData(ParticleID::Wplus)->mass()/(1. - sqr(theSw));
  
  theSinA = sin(theMSSM->higgsMixingAngle());
  theCosA = sqrt(1. - sqr(theSinA));
  theTanB = theMSSM->tanBeta();
  theSinB = theTanB/sqrt(1. + sqr(theTanB));
  theCosB = sqrt( 1. - sqr(theSinB) );
  theSinApB = theSinA*theCosB + theCosA*theSinB;
  theCosApB = theCosA*theCosB - theSinA*theSinB;
  
  MixingMatrixPtr stop = theMSSM->stopMix();
  MixingMatrixPtr sbot = theMSSM->sbottomMix();
  theQt1L  = (*stop)(0,0)*(*stop)(0,0);
  theQt1R  = (*stop)(0,1)*(*stop)(0,1);
  theQt1LR = (*stop)(0,1)*(*stop)(0,0) + (*stop)(0,1)*(*stop)(0,0);
  theQt2L  = (*stop)(1,0)*(*stop)(1,0);
  theQt2R  = (*stop)(1,1)*(*stop)(1,1);
  theQt2LR = (*stop)(1,1)*(*stop)(1,0) + (*stop)(1,0)*(*stop)(1,1);
  theQb1L  = (*sbot)(0,0)*(*sbot)(0,0);
  theQb1R  = (*sbot)(0,1)*(*sbot)(0,1);
  theQb1LR = (*sbot)(0,1)*(*sbot)(0,0) + (*sbot)(0,1)*(*sbot)(0,0);
  theQb2L  = (*sbot)(1,0)*(*sbot)(1,0);
  theQb2R  = (*sbot)(1,1)*(*sbot)(1,1);
  theQb2LR = (*sbot)(1,1)*(*sbot)(1,0) + (*sbot)(1,0)*(*sbot)(1,1);
  
  assert( theSqmass.size() == 4 );
  theSqmass[0] = getParticleData(ParticleID::SUSY_b_1)->mass();
  theSqmass[1] = getParticleData(ParticleID::SUSY_t_1)->mass();
  theSqmass[2] = getParticleData(ParticleID::SUSY_b_2)->mass();
  theSqmass[3] = getParticleData(ParticleID::SUSY_t_2)->mass();

  VVSLoopVertex::doinit();
  if(loopToolsInitialized()) Looptools::ltexi();
}


void SSHGGVertex::persistentOutput(PersistentOStream & os) const {
  os << theMSSM << theSw << ounit(theMw,GeV) << ounit(theZfact,GeV) 
     << theQt1L << theQt1R << theQt1LR << theQt2L << theQt2R << theQt2LR
     << theQb1L << theQb1R << theQb1LR << theQb2L << theQb2R << theQb2LR
     << thetop << thebot << theTanB << theIncludeTriLinear
     << theSinA << theCosA << theSinB << theCosB << theSinApB << theCosApB
     << ounit(theSqmass, GeV) << thePseudoScalarTreatment;
}

void SSHGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMSSM >> theSw >> iunit(theMw,GeV) >> iunit(theZfact,GeV) 
     >> theQt1L >> theQt1R >> theQt1LR >> theQt2L >> theQt2R >> theQt2LR 
     >> theQb1L >> theQb1R >> theQb1LR >> theQb2L >> theQb2R >> theQb2LR 
     >> thetop >> thebot >> theTanB >> theIncludeTriLinear
     >> theSinA >> theCosA >> theSinB >> theCosB >> theSinApB >> theCosApB
     >> iunit(theSqmass, GeV) >> thePseudoScalarTreatment;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSHGGVertex,VVSLoopVertex>
describeHerwigSSHGGVertex("Herwig::SSHGGVertex", "HwSusy.so");

void SSHGGVertex::Init() {
  
  static ClassDocumentation<SSHGGVertex> documentation
    ("This class implements the higgs-gluon-gluon effective "
     "vertex in the MSSM including stop, sbottom and top quarks " 
     "loops.");

  static Switch<SSHGGVertex,bool> interfaceIncludeTriLinear
    ("IncludeTriLinear",
     "Whether or not to include the A term squark trilinear couplings",
     &SSHGGVertex::theIncludeTriLinear, true, false, false);
  static SwitchOption interfaceIncludeTriLinearYes
    (interfaceIncludeTriLinear,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncludeTriLinearNo
    (interfaceIncludeTriLinear,
     "No",
     "Don't include them",
     false);

  static Switch<SSHGGVertex,bool> interfacePseudoScalarTreatment
    ("PseudoScalarTreatment",
     "Whether to treat the pseudoscalar as pseudoscalar or scalar, for testing only",
     &SSHGGVertex::thePseudoScalarTreatment, false, false, false);
  static SwitchOption interfacePseudoScalarTreatmentPseudoScalar
    (interfacePseudoScalarTreatment,
     "PseudoScalar",
     "Treat as a pseudoscalar, the right physics",
     false);
  static SwitchOption interfacePseudoScalarTreatmentScalar
    (interfacePseudoScalarTreatment,
     "Scalar",
     "Treat as a scalar for testing",
     true);
}
  
void SSHGGVertex::setCoupling(Energy2 q2, tcPDPtr particle2,
			      tcPDPtr particle3, tcPDPtr particle1) {
  long higgs(abs(particle1->id()));
  assert( higgs == ParticleID::h0 || higgs == ParticleID::H0 ||
	  higgs == ParticleID::A0 );
  assert(particle2->id() == ParticleID::g && particle3->id() == ParticleID::g );
  if( q2 != theq2last || theCouplast == 0. || higgs != theLastID ) {
    Looptools::clearcache();
    theCouplast = weakCoupling(q2)*sqr(strongCoupling(q2));
    Energy mt = theMSSM->mass(q2, thetop);    
    Energy mb = theMSSM->mass(q2, thebot);
    masses.clear();
    type.clear();
    if( higgs == ParticleID::h0 || higgs == ParticleID::H0 ) {
      setNParticles(6);
      masses.insert(masses.begin(), theSqmass.begin(), theSqmass.end());
      masses.push_back(mt);
      masses.push_back(mb);
      type.resize(6, PDT::Spin0);
      type[4] = PDT::Spin1Half;
      type[5] = PDT::Spin1Half;
      couplings.resize(6, make_pair(0., 0.));  
      complex<Energy> brac1 = theZfact*(0.5 + theMSSM->ed()*sqr(theSw));
      complex<Energy> brac2 = theZfact*(0.5 - theMSSM->eu()*sqr(theSw));
      complex<Energy> brac3 = theZfact*theMSSM->ed()*sqr(theSw);
      complex<Energy> brac4 = theZfact*theMSSM->eu()*sqr(theSw);
      Energy Trib=theMSSM->bottomTrilinear().real();
      Energy Trit=theMSSM->topTrilinear().real();
      Energy theMu = theMSSM->muParameter();
      if( higgs == ParticleID::h0 ) {
	// lightest sbottom
	complex<Energy> trilinear = theIncludeTriLinear ? 
	  theQb1LR*0.5*mb/theMw*(Trib*theSinA + theMu*theCosA)/theCosB : Energy();
	Complex coup = 0.5*UnitRemoval::InvE*
	  (theQb1L *(   sqr(mb)*theSinA/theMw/theCosB - theSinApB*brac1) +
	   theQb1R *(   sqr(mb)*theSinA/theMw/theCosB + theSinApB*brac3)+trilinear);

	couplings[0] = make_pair(coup, coup);
	// lightest stop
	trilinear = theIncludeTriLinear ? 
	 -theQt1LR*0.5*mt/theMw*(Trit*theCosA + theMu*theSinA)/theSinB : Energy();
	coup = Complex(0.5*UnitRemoval::InvE*
		       (theQt1L *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac2) +
			theQt1R *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac4)+trilinear));

	couplings[1] = make_pair(coup, coup);
	// heavier sbottom
	trilinear = theIncludeTriLinear ? 
	  theQb2LR*0.5*mb/theMw*(Trib*theSinA + theMu*theCosA)/theCosB : Energy();
	coup = Complex(0.5*UnitRemoval::InvE*
		       (theQb2L *(   sqr(mb)*theSinA/theMw/theCosB - theSinApB*brac1) +
			theQb2R *(   sqr(mb)*theSinA/theMw/theCosB + theSinApB*brac3)+trilinear));

	couplings[2] = make_pair(coup, coup);
	// heavier stop
	trilinear = theIncludeTriLinear ? 
	  -theQt2LR*0.5*mt/theMw*(Trit*theCosA + theMu*theSinA)/theSinB : Energy();
	coup = Complex(0.5*UnitRemoval::InvE*
		       (theQt2L *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac2) +
			theQt2R *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac4)+trilinear));
	   		
	couplings[3] = make_pair(coup, coup);
	// top
	coup = -0.25*(mt*theCosA/theMw/theSinB);

	couplings[4] = make_pair(coup, coup);
	// bottom
	coup = +0.25*(mb*theSinA/theMw/theCosB);

	couplings[5] = make_pair(coup, coup);
      }
      else {
	// lightest sbottom
	complex<Energy> trilinear = theIncludeTriLinear ? 
	   theQb1LR*0.5*mb/theMw*(theMu*theSinA - Trib*theCosA)/theCosB: Energy();
	Complex coup = 0.5*UnitRemoval::InvE*
	  (theQb1L *( - sqr(mb)*theCosA/theMw/theCosB + theCosApB*brac1) +
	   theQb1R *( - sqr(mb)*theCosA/theMw/theCosB - theCosApB*brac3)+trilinear);

	couplings[0] = make_pair(coup, coup);
	// lightest stop
	trilinear = theIncludeTriLinear ? 
	   -theQt1LR*0.5*mt/theMw*(-theMu*theCosA + Trit*theSinA)/theSinB: Energy();
	coup = Complex(0.5*UnitRemoval::InvE*
		       (theQt1L *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac2) +
			theQt1R *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac4)+trilinear));

	couplings[1] = make_pair(coup, coup);
	// heavier sbottom
	trilinear = theIncludeTriLinear ? 
	   theQb2LR*0.5*mb/theMw*(theMu*theSinA - Trib*theCosA)/theCosB: Energy();
	coup =  Complex(0.5*UnitRemoval::InvE*
			(theQb2L *( - sqr(mb)*theCosA/theMw/theCosB + theCosApB*brac1) +
			 theQb2R *( - sqr(mb)*theCosA/theMw/theCosB - theCosApB*brac3)+trilinear)); 
	
	couplings[2] = make_pair(coup, coup);
	// heavier stop
	trilinear = theIncludeTriLinear ? 
	   -theQt2LR*0.5*mt/theMw*(-theMu*theCosA + Trit*theSinA)/theSinB: Energy();
	coup =  Complex(0.5*UnitRemoval::InvE*
			(theQt2L *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac2) +
			 theQt2R *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac4)+trilinear));
	
	couplings[3] = make_pair(coup, coup);
	// top
	coup = -0.25*mt*theSinA/theMw/theSinB;
	couplings[4] = make_pair(coup, coup);
	// bottom
	coup = -0.25*mb*theCosA/theMw/theCosB;
	couplings[5] = make_pair(coup, coup);
      }
    }
    else {
      setNParticles(2);
      masses.resize(2);
      couplings.resize(2);
      masses[0] = mt;
      masses[1] = mb;
      type.resize(2,PDT::Spin1Half);
      Complex coup = 0.25*Complex(0., 1.)*mt/theMw/theTanB;
	  	
      couplings[0] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
      coup = 0.25*Complex(0., 1.)*mb/theMw*theTanB;
	 
      couplings[1] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
    }
    theq2last = q2;
    theLastID = higgs;
    theHaveCoeff = false;
  }
  norm(theCouplast);
  //calculate tensor coefficients
  if( !theHaveCoeff ) {
    VVSLoopVertex::setCoupling(q2, particle2, particle3, particle1);
    theHaveCoeff = true;
  }
}

 
#line 1 "./SSHPPVertex.cc"
// -*- C++ -*-
//
// SSHPPVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSHPPVertex class.
//

#include "SSHPPVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>
#include "Herwig/Looptools/clooptools.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSHPPVertex::SSHPPVertex() : theIncludeTriLinear(true),
			     thePseudoScalarTreatment(false),
			     theSw(0.), theMw(), theZfact(),
			     theQt1L(0.), theQt1R(0.), theQt1LR(0.),
			     theQt2L(0.), theQt2R(0.), theQt2LR(0.),
			     theQb1L(0.), theQb1R(0.), theQb1LR(0.),
			     theQb2L(0.), theQb2R(0.), theQb2LR(0.),
			     theLt1L(0.), theLt1R(0.), theLt1LR(0.),
			     theLt2L(0.), theLt2R(0.), theLt2LR(0.),
			     theSfmass(6,ZERO),
			     theTanB(0.),theSinA(0.), 
			     theCosA(0.), theSinB(0.), theCosB(0.), 
			     theSinApB(0.), theCosApB(0.), 
			     theSinBmA(0.), theCosBmA(0.), 
			     theCouplast(0.), 
			     theq2last(), theHaveCoeff(false), theLastID(0) {
  orderInGs(0);
  orderInGem(3);
  colourStructure(ColourStructure::SINGLET);
}

void SSHPPVertex::persistentOutput(PersistentOStream & os) const {
  os << theMSSM << theSw << ounit(theMw,GeV) << ounit(theZfact,GeV) 
     << theQt1L << theQt1R << theQt1LR << theQt2L << theQt2R << theQt2LR
     << theQb1L << theQb1R << theQb1LR << theQb2L << theQb2R << theQb2LR
     << theLt1L << theLt1R << theLt1LR << theLt2L << theLt2R << theLt2LR
     << thetop << thebot << thetau << theTanB << theIncludeTriLinear
     << theSinA << theCosA << theSinB << theCosB << theSinApB << theCosApB
     << theSinBmA << theCosBmA << thePseudoScalarTreatment
     << ounit(theSfmass, GeV) << theU << theV;
}

void SSHPPVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMSSM >> theSw >> iunit(theMw,GeV) >> iunit(theZfact,GeV) 
     >> theQt1L >> theQt1R >> theQt1LR >> theQt2L >> theQt2R >> theQt2LR 
     >> theQb1L >> theQb1R >> theQb1LR >> theQb2L >> theQb2R >> theQb2LR
     >> theLt1L >> theLt1R >> theLt1LR >> theLt2L >> theLt2R >> theLt2LR
     >> thetop >> thebot >> thetau >> theTanB >> theIncludeTriLinear
     >> theSinA >> theCosA >> theSinB >> theCosB >> theSinApB >> theCosApB
     >> theSinBmA >> theCosBmA >> thePseudoScalarTreatment
     >> iunit(theSfmass, GeV) >> theU >> theV;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSHPPVertex,VVSLoopVertex>
describeHerwigSSHPPVertex("Herwig::SSHPPVertex", "HwSusy.so");

void SSHPPVertex::Init() {
  
  static ClassDocumentation<SSHPPVertex> documentation
    ("This class implements the higgs-gluon-gluon effective "
     "vertex in the MSSM including stop, sbottom and top quarks " 
     "loops.");

  static Switch<SSHPPVertex,bool> interfaceIncludeTriLinear
    ("IncludeTriLinear",
     "Whether or not to include the A term squark trilinear couplings",
     &SSHPPVertex::theIncludeTriLinear, true, false, false);
  static SwitchOption interfaceIncludeTriLinearYes
    (interfaceIncludeTriLinear,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncludeTriLinearNo
    (interfaceIncludeTriLinear,
     "No",
     "Don't include them",
     false);

  static Switch<SSHPPVertex,bool> interfacePseudoScalarTreatment
    ("PseudoScalarTreatment",
     "Whether to treat the pseudoscalar as pseudoscalar or scalar, for testing only",
     &SSHPPVertex::thePseudoScalarTreatment, false, false, false);
  static SwitchOption interfacePseudoScalarTreatmentPseudoScalar
    (interfacePseudoScalarTreatment,
     "PseudoScalar",
     "Treat as a pseudoscalar, the right physics",
     false);
  static SwitchOption interfacePseudoScalarTreatmentScalar
    (interfacePseudoScalarTreatment,
     "Scalar",
     "Treat as a scalar for testing",
     true);

}
  
void SSHPPVertex::setCoupling(Energy2 q2, tcPDPtr particle2,
			      tcPDPtr particle3, tcPDPtr particle1) {
  long higgs(abs(particle1->id()));
  // check allowed
  assert ( higgs == ParticleID::h0 || higgs == ParticleID::H0 ||
	   higgs == ParticleID::A0 );
  assert(particle2->id() == ParticleID::gamma && 
	 particle3->id() == ParticleID::gamma );
  // couplings
  if( q2 != theq2last || theCouplast == 0. || higgs != theLastID ) {
    Looptools::clearcache();
    theCouplast = sqr(electroMagneticCoupling(q2))*weakCoupling(q2);
    Energy mt   = theMSSM->mass(q2, thetop);    
    Energy mb   = theMSSM->mass(q2, thebot);
    Energy mtau = theMSSM->mass(q2, thetau);
    masses.clear();
    type.clear();
    if( higgs == ParticleID::h0 || higgs == ParticleID::H0 ) {
      setNParticles(13);
      masses.insert(masses.begin(), theSfmass.begin(), theSfmass.end());
      masses.push_back(mt);
      masses.push_back(mb);
      masses.push_back(mtau);
      masses.push_back(getParticleData(ParticleID::Hplus)->mass());
      masses.push_back(theMw);
      masses.push_back(getParticleData(ParticleID::SUSY_chi_1plus)->mass());
      masses.push_back(getParticleData(ParticleID::SUSY_chi_2plus)->mass());
      type.resize(13, PDT::Spin0);
      type[6] = PDT::Spin1Half;
      type[7] = PDT::Spin1Half;
      type[8] = PDT::Spin1Half;
      type[9] = PDT::Spin0;
      type[10] = PDT::Spin1;
      type[11] = PDT::Spin1Half;
      type[12] = PDT::Spin1Half;
      couplings.resize(13, make_pair(0., 0.));
      complex<Energy> brac1 = theZfact*(0.5 + theMSSM->ed()*sqr(theSw));
      complex<Energy> brac2 = theZfact*(0.5 - theMSSM->eu()*sqr(theSw));
      complex<Energy> brac3 = theZfact*theMSSM->ed()*sqr(theSw);
      complex<Energy> brac4 = theZfact*theMSSM->eu()*sqr(theSw);
      complex<Energy> brac5 = theZfact*(0.5 + theMSSM->ee()*sqr(theSw));
      complex<Energy> brac6 = theZfact*theMSSM->ee()*sqr(theSw);
      Energy Trib=theMSSM->bottomTrilinear().real();
      Energy Trit=theMSSM->topTrilinear().real();
      Energy Trita=theMSSM->tauTrilinear().real();
      Energy theMu = theMSSM->muParameter();
      if( higgs == ParticleID::h0 ) {
	// lightest sbottom
	complex<Energy> trilinear = theIncludeTriLinear ? 
	  theQb1LR*0.5*mb/theMw*(Trib*theSinA + theMu*theCosA)/theCosB : Energy();
	Complex coup = 3.*UnitRemoval::InvE*sqr(theMSSM->ed())*
	  (theQb1L *(   sqr(mb)*theSinA/theMw/theCosB - theSinApB*brac1) +
	   theQb1R *(   sqr(mb)*theSinA/theMw/theCosB + theSinApB*brac3) +
	   trilinear);
	couplings[0] = make_pair(coup, coup);
	// lightest stop
	trilinear = theIncludeTriLinear ?
	  theQt1LR*0.5*mt/theMw*(Trit*theCosA + theMu*theSinA)/theSinB : Energy();
	coup = Complex(3.*UnitRemoval::InvE*sqr(theMSSM->eu())*
		       (theQt1L *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac2) +
			theQt1R *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac4) -
			trilinear));
	couplings[1] = make_pair(coup, coup);
	// lightest stau
	trilinear = theIncludeTriLinear ?
	  theLt1LR*0.5*mtau/theMw*(Trita*theSinA + theMu*theCosA)/theCosB : Energy();
	coup = Complex(UnitRemoval::InvE*sqr(theMSSM->ee())*
		       (theLt1L *(   sqr(mtau)*theSinA/theMw/theCosB - theSinApB*brac5) +
			theLt1R *(   sqr(mtau)*theSinA/theMw/theCosB + theSinApB*brac6) +
			trilinear));
	couplings[2] = make_pair(coup, coup);
	// heavier sbottom
	trilinear = theIncludeTriLinear ? 
	   theQb2LR*0.5*mb/theMw*(Trib*theSinA + theMu*theCosA)/theCosB : Energy();
	coup = Complex(3.*UnitRemoval::InvE*sqr(theMSSM->ed())*
		       (theQb2L *(   sqr(mb)*theSinA/theMw/theCosB - theSinApB*brac1) +
			theQb2R *(   sqr(mb)*theSinA/theMw/theCosB + theSinApB*brac3) +
			trilinear));
	couplings[3] = make_pair(coup, coup);
	// heavier stop
	trilinear = theIncludeTriLinear ? 
	  theQt2LR*0.5*mt/theMw*(Trit*theCosA + theMu*theSinA)/theSinB : Energy();
	coup = Complex(3.*UnitRemoval::InvE*sqr(theMSSM->eu())*
		       (theQt2L*( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac2) +
			theQt2R*( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac4) -
			trilinear));
	couplings[4] = make_pair(coup, coup);
	// heavier stau
	trilinear = theIncludeTriLinear ? 
	  theLt2LR*0.5*mtau/theMw*(Trita*theSinA + theMu*theCosA)/theCosB : Energy();
	coup = Complex(UnitRemoval::InvE*sqr(theMSSM->ee())*
		       (theLt2L *(   sqr(mtau)*theSinA/theMw/theCosB - theSinApB*brac5) +
			theLt2R *(   sqr(mtau)*theSinA/theMw/theCosB + theSinApB*brac6)+
			trilinear));
	couplings[5] = make_pair(coup, coup);
	// top
	coup = - 3.*mt*sqr(theMSSM->eu())*theCosA/2./theMw/theSinB;
	couplings[6] = make_pair(coup, coup);
	// bottom
	coup =   3.*mb*sqr(theMSSM->ed())*theSinA/2./theMw/theCosB;
	couplings[7] = make_pair(coup, coup);
	// tau
	coup =    mtau*sqr(theMSSM->ee())*theSinA/2./theMw/theCosB;
	couplings[8] = make_pair(coup, coup);
	// charged higgs
	coup = - UnitRemoval::InvE*theMw*(theSinBmA + 0.5/(1.-sqr(theSw))*
		       (sqr(theCosB)-sqr(theSinB))*theSinApB);
	couplings[9] = make_pair(coup, coup);
	// W boson
	coup = UnitRemoval::InvE*theMw*theSinBmA;
	couplings[10] = make_pair(coup, coup);
	// charginos
	for(unsigned int ix=0;ix<2;++ix) {
	  Complex Q = sqrt(0.5)*(*theV)(ix,0)*(*theU)(ix,1);
	  Complex S = sqrt(0.5)*(*theV)(ix,1)*(*theU)(ix,0);
	  coup = Q*theSinA-S*theCosA;
	  couplings[11+ix] = make_pair(conj(coup), coup);
	}
      }
      else {
	// lightest sbottom
	complex<Energy> trilinear = theIncludeTriLinear ? 
	  theQb1LR*0.5*mb/theMw*(theMu*theSinA - Trib*theCosA)/theCosB: Energy();
	Complex coup = 3.*UnitRemoval::InvE*sqr(theMSSM->ed())*
	  (theQb1L *( - sqr(mb)*theCosA/theMw/theCosB + theCosApB*brac1) +
	   theQb1R *( - sqr(mb)*theCosA/theMw/theCosB - theCosApB*brac3)+trilinear);
	couplings[0] = make_pair(coup, coup);
	// lightest stop
	trilinear = theIncludeTriLinear ? 
	   -theQt1LR*0.5*mt/theMw*(-theMu*theCosA + Trit*theSinA)/theSinB: Energy();
	coup = Complex(3.*UnitRemoval::InvE*sqr(theMSSM->eu())*
		       (theQt1L *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac2) +
			theQt1R *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac4)+trilinear));
	couplings[1] = make_pair(coup, coup);
	// lightest stau
	trilinear = theIncludeTriLinear ? 
	   theLt1LR*0.5*mtau/theMw*(theMu*theSinA - Trita*theCosA)/theCosB: Energy();
	coup = Complex(UnitRemoval::InvE*sqr(theMSSM->ee())*
		       (theLt1L *( - sqr(mtau)*theCosA/theMw/theCosB + theCosApB*brac5) +
			theLt1R *( - sqr(mtau)*theCosA/theMw/theCosB - theCosApB*brac6)+trilinear));
	couplings[2] = make_pair(coup, coup);
	// heavier sbottom
	trilinear = theIncludeTriLinear ? 
	   theQb2LR*0.5*mb/theMw*(theMu*theSinA - Trib*theCosA)/theCosB: Energy();
	coup = Complex(3.*UnitRemoval::InvE*sqr(theMSSM->ed())*
		       (theQb2L *( - sqr(mb)*theCosA/theMw/theCosB + theCosApB*brac1) +
			theQb2R *( - sqr(mb)*theCosA/theMw/theCosB - theCosApB*brac3)+trilinear)); 
	couplings[3] = make_pair(coup, coup);
	// heavier stop
	trilinear = theIncludeTriLinear ? 
	  -theQt2LR*0.5*mt/theMw*(-theMu*theCosA + Trit*theSinA)/theSinB: Energy();
	coup = Complex(3.*UnitRemoval::InvE*sqr(theMSSM->eu())*
		       (theQt2L *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac2) +
			theQt2R *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac4)+trilinear));
	couplings[4] = make_pair(coup, coup);
	// heavier stau
	trilinear = theIncludeTriLinear ? 
	   theLt2LR*0.5*mtau/theMw*(theMu*theSinA - Trita*theCosA)/theCosB: Energy();
	coup = Complex(UnitRemoval::InvE*sqr(theMSSM->ee())*
		       (theLt2L *( - sqr(mtau)*theCosA/theMw/theCosB + theCosApB*brac5) +
			theLt2R *( - sqr(mtau)*theCosA/theMw/theCosB - theCosApB*brac6)+trilinear));
	couplings[5] = make_pair(coup, coup);
	// top
	coup = -3.*mt*sqr(theMSSM->eu())*theSinA/2./theMw/theSinB;
	couplings[6] = make_pair(coup, coup);
	// bottom
	coup = -3.*mb*sqr(theMSSM->ed())*theCosA/2./theMw/theCosB;
	couplings[7] = make_pair(coup, coup);
	// tau
	coup = -mtau*sqr(theMSSM->ee())*theCosA/2./theMw/theCosB;
	couplings[8] = make_pair(coup, coup);
	// charged higgs
	coup = - UnitRemoval::InvE*theMw*(theCosBmA - 0.5/(1.-sqr(theSw))*
					  (sqr(theCosB)-sqr(theSinB))*theCosApB);
	couplings[9] = make_pair(coup, coup);
	// W boson
	coup = UnitRemoval::InvE*theMw*theCosBmA;
	couplings[10] = make_pair(coup, coup);
	// charginos
	for(unsigned int ix=0;ix<2;++ix) {
	  Complex Q = sqrt(0.5)*(*theV)(ix,0)*(*theU)(ix,1);
	  Complex S = sqrt(0.5)*(*theV)(ix,1)*(*theU)(ix,0);
	  coup = -Q*theCosA-S*theSinA;
	  couplings[11+ix] = make_pair(conj(coup), coup);
	}
      }
    }
    else {
      setNParticles(5);
      masses.resize(5);
      couplings.resize(5);
      masses[0] = mt;
      masses[1] = mb;
      masses[2] = mtau;
      masses[3] = getParticleData(ParticleID::SUSY_chi_1plus)->mass();
      masses[4] = getParticleData(ParticleID::SUSY_chi_2plus)->mass();
      type.resize(5,PDT::Spin1Half);
      // top
      Complex coup = 3.*Complex(0., 1.)*sqr(theMSSM->eu())*mt/2./theMw/theTanB;
      couplings[0] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
      // bottom
      coup = 3.*Complex(0., 1.)*sqr(theMSSM->ed())*mb/2./theMw*theTanB;
      couplings[1] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
      // tau
      coup = Complex(0., 1.)*sqr(theMSSM->ee())*mtau/2./theMw*theTanB;
      couplings[2] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
      // charginos
      for(unsigned int ix=0;ix<2;++ix) {
	Complex Q = sqrt(0.5)*(*theV)(ix,0)*(*theU)(ix,1);
	Complex S = sqrt(0.5)*(*theV)(ix,1)*(*theU)(ix,0);
	coup = - Complex(0., 1.)*(Q*theSinB+S*theCosB);
	couplings[3+ix] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
      }
    }
    theq2last = q2;
    theLastID = higgs;
    theHaveCoeff = false;
  }
  norm(theCouplast);
  //calculate tensor coefficients
  if( !theHaveCoeff ) {
    VVSLoopVertex::setCoupling(q2, particle2, particle3, particle1);
    theHaveCoeff = true;
  }
}

// functions for loops for testing
// namespace {

// Complex F0(double tau) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return tau*(1.-tau*ft);
// }

// Complex FHalf(double tau,double eta) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return -2.*tau*(eta+(1.-tau*eta)*ft);
// }

// Complex F1(double tau) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return 2.+3.*tau+3.*tau*(2.-tau)*ft;
// }
// }

void SSHPPVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(22,22,25);
  addToList(22,22,35);
  addToList(22,22,36);
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM ) 
    throw InitException()
      << "SSHPPVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  theMw = getParticleData(ParticleID::Wplus)->mass();
  thetop = getParticleData(ParticleID::t);
  thebot = getParticleData(ParticleID::b);
  thetau = getParticleData(ParticleID::tauminus);
  theSw = sqrt(sin2ThetaW());
  theZfact = getParticleData(ParticleID::Wplus)->mass()/(1. - sqr(theSw));
  
  theSinA = sin(theMSSM->higgsMixingAngle());
  theCosA = sqrt(1. - sqr(theSinA));
  theTanB = theMSSM->tanBeta();
  theSinB = theTanB/sqrt(1. + sqr(theTanB));
  theCosB = sqrt( 1. - sqr(theSinB) );
  theSinApB = theSinA*theCosB + theCosA*theSinB;
  theCosApB = theCosA*theCosB - theSinA*theSinB;
  theSinBmA =-theSinA*theCosB + theCosA*theSinB;
  theCosBmA = theCosA*theCosB + theSinA*theSinB;
  
  MixingMatrix stop = *theMSSM->stopMix();
  MixingMatrix sbot = *theMSSM->sbottomMix();
  MixingMatrix stau = *theMSSM->stauMix();
  theQt1L  = stop(0,0)*stop(0,0);
  theQt1R  = stop(0,1)*stop(0,1);
  theQt1LR = stop(0,1)*stop(0,0) + stop(0,1)*stop(0,0);
  theQt2L  = stop(1,0)*stop(1,0);
  theQt2R  = stop(1,1)*stop(1,1);
  theQt2LR = stop(1,1)*stop(1,0) + stop(1,0)*stop(1,1);
  theQb1L  = sbot(0,0)*sbot(0,0);
  theQb1R  = sbot(0,1)*sbot(0,1);
  theQb1LR = sbot(0,1)*sbot(0,0) + sbot(0,1)*sbot(0,0);
  theQb2L  = sbot(1,0)*sbot(1,0);
  theQb2R  = sbot(1,1)*sbot(1,1);
  theQb2LR = sbot(1,1)*sbot(1,0) + sbot(1,0)*sbot(1,1);
  theLt1L  = stau(0,0)*stau(0,0);
  theLt1R  = stau(0,1)*stau(0,1);
  theLt1LR = stau(0,1)*stau(0,0) + stau(0,1)*stau(0,0);
  theLt2L  = stau(1,0)*stau(1,0);
  theLt2R  = stau(1,1)*stau(1,1);
  theLt2LR = stau(1,1)*stau(1,0) + stau(1,0)*stau(1,1);
  theU = theMSSM->charginoUMix();
  theV = theMSSM->charginoVMix();

  assert( theSfmass.size() == 6 );
  theSfmass[0] = getParticleData(ParticleID::SUSY_b_1)->mass();
  theSfmass[1] = getParticleData(ParticleID::SUSY_t_1)->mass();
  theSfmass[2] = getParticleData(ParticleID::SUSY_tau_1minus)->mass();
  theSfmass[3] = getParticleData(ParticleID::SUSY_b_2)->mass();
  theSfmass[4] = getParticleData(ParticleID::SUSY_t_2)->mass();
  theSfmass[5] = getParticleData(ParticleID::SUSY_tau_2minus)->mass();

  VVSLoopVertex::doinit();
  // test calc of the width
//   for(unsigned int ix=0;ix<2;++ix) {
//     Energy mh   = getParticleData(25+long(ix)*10)->mass();
//     Energy mt   = theMSSM->mass(sqr(mh  ), thetop);    
//     Energy mb   = theMSSM->mass(sqr(mh  ), thebot);    
//     Energy mtau = theMSSM->mass(sqr(mh  ), thetau);
//     Energy mhp  = getParticleData(ParticleID::Hplus)->mass();
//     Energy mc[2] = {getParticleData(ParticleID::SUSY_chi_1plus)->mass(),
// 		    getParticleData(ParticleID::SUSY_chi_2plus)->mass()};
//     // sbottom
//     Complex rsb1,rsb2;
//     if(ix==0) {
//       rsb1 = 
// 	+theQb1L*(-sqr(mb/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  -(-0.5+sqr(theSw)/3.)*theSinApB)
// 	+theQb1R*(-sqr(mb/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  +sqr(theSw)/3.*theSinApB);
//       rsb2 = 
// 	+theQb2L*(-sqr(mb/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  -(-0.5+sqr(theSw)/3.)*theSinApB)
// 	+theQb2R*(-sqr(mb/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  +sqr(theSw)/3.*theSinApB);
//     }
//     else {
//       rsb1 = 
// 	+theQb1L*(+sqr(mb/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  +(-0.5+sqr(theSw)/3.)*theCosApB)
// 	+theQb1R*(+sqr(mb/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  -sqr(theSw)/3.*theCosApB);
//       rsb2 = 
// 	+theQb2L*(+sqr(mb/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  +(-0.5+sqr(theSw)/3.)*theCosApB)
// 	+theQb2R*(+sqr(mb/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  -sqr(theSw)/3.*theCosApB);
//     }
//     Complex Isb1 = 3.*sqr(1./3.)*rsb1*sqr(theMw/theSfmass[0])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[0]/mh));
//     Complex Isb2 = 3.*sqr(1./3.)*rsb2*sqr(theMw/theSfmass[3])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[3]/mh));
//     // stop
//     Complex rst1,rst2;
//     if(ix==0) {
//       rst1 = 
// 	+theQt1L*(+sqr(mt/theMw)*(1.-sqr(theSw))*theCosA/theSinB
// 		  -(+0.5-2.*sqr(theSw)/3.)*theSinApB)
// 	+theQt1R*(+sqr(mt/theMw)*(1.-sqr(theSw))*theCosA/theSinB
// 		  -2.*sqr(theSw)/3.*theSinApB);
//       rst2 = 
// 	+theQt2L*(+sqr(mt/theMw)*(1.-sqr(theSw))*theCosA/theSinB
// 		  -(+0.5-2.*sqr(theSw)/3.)*theSinApB)
// 	+theQt2R*(+sqr(mt/theMw)*(1.-sqr(theSw))*theCosA/theSinB
// 		  -2.*sqr(theSw)/3.*theSinApB);
//     }
//     else {
//       rst1 = 
// 	+theQt1L*(+sqr(mt/theMw)*(1.-sqr(theSw))*theSinA/theSinB
// 		  +(+0.5-2.*sqr(theSw)/3.)*theCosApB)
// 	+theQt1R*(+sqr(mt/theMw)*(1.-sqr(theSw))*theSinA/theSinB
// 		  +2.*sqr(theSw)/3.*theCosApB);
//       rst2 = 
// 	+theQt2L*(+sqr(mt/theMw)*(1.-sqr(theSw))*theSinA/theSinB
// 		  +(+0.5-2.*sqr(theSw)/3.)*theCosApB)
// 	+theQt2R*(+sqr(mt/theMw)*(1.-sqr(theSw))*theSinA/theSinB
// 		  +2.*sqr(theSw)/3.*theCosApB);
//     }
//     Complex Ist1 = 3.*sqr(2./3.)*rst1*sqr(theMw/theSfmass[1])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[1]/mh));
//     Complex Ist2 = 3.*sqr(2./3.)*rst2*sqr(theMw/theSfmass[4])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[4]/mh));

//     // stau
//     Complex rstau1,rstau2;
//     if(ix==0) {
//       rstau1 = 
// 	+theLt1L*(-sqr(mtau/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  -(-0.5+sqr(theSw))*theSinApB)
// 	+theLt1R*(-sqr(mtau/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  +sqr(theSw)*theSinApB);
//       rstau2 = 
// 	+theLt2L*(-sqr(mtau/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  -(-0.5+sqr(theSw))*theSinApB)
// 	+theLt2R*(-sqr(mtau/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  +sqr(theSw)*theSinApB);
//     }
//     else {
//       rstau1 = 
// 	+theLt1L*(+sqr(mtau/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  +(-0.5+sqr(theSw))*theCosApB)
// 	+theLt1R*(+sqr(mtau/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  -sqr(theSw)*theCosApB);
//       rstau2 = 
// 	+theLt2L*(+sqr(mtau/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  +(-0.5+sqr(theSw))*theCosApB)
// 	+theLt2R*(+sqr(mtau/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  -sqr(theSw)*theCosApB);
//     }
//     Complex Istau1 = rstau1*sqr(theMw/theSfmass[2])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[2]/mh));
//     Complex Istau2 = rstau2*sqr(theMw/theSfmass[5])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[5]/mh));
//     // charged higgs
//     Complex rh;
//     if(ix==0) {
//       rh = theSinBmA+0.5*(sqr(theCosB)-sqr(theSinB))*theSinApB/(1.-sqr(theSw));
//     }
//     else {
//       rh = theCosBmA-0.5*(sqr(theCosB)-sqr(theSinB))*theCosApB/(1.-sqr(theSw));
//     }
//     Complex Ih = rh*sqr(theMw/mhp)*F0(sqr(2.*mhp/mh));
//     // W
//     Complex rw;
//     if(ix==0) {
//       rw = theSinBmA;
//     }
//     else {
//       rw = theCosBmA;
//     }
//     Complex IW = rw*F1(sqr(2.*theMw/mh));
//     // top
//     Complex rt;
//     if(ix==0) {
//       rt = theCosA/theSinB;
//     }
//     else {
//       rt = theSinA/theSinB;
//     }
//     Complex Itop = 3.*sqr(2./3.)*rt*FHalf(sqr(2.*mt/mh),1.);
//     // bottom
//     Complex rb;
//     if(ix==0) {
//       rb =-theSinA/theCosB;
//     }
//     else {
//       rb = theCosA/theCosB;
//     }
//     Complex Ibot = 3.*sqr(1./3.)*rb*FHalf(sqr(2.*mb/mh),1.);
//     // tau
//     Complex rtau;
//     if(ix==0) {
//       rtau =-theSinA/theCosB;
//     }
//     else {
//       rtau = theCosA/theCosB;
//     }
//     Complex Itau = rtau*FHalf(sqr(2.*mtau/mh),1.);
//     // charginos
//     Complex rc[2],IC[2];
//     for(unsigned int ic=0;ic<2;++ic) {
//       Complex Q = sqrt(0.5)*(*theV)(ic,0)*(*theU)(ic,1);
//       Complex S = sqrt(0.5)*(*theV)(ic,1)*(*theU)(ic,0);
//       if(ix==0) {
// 	rc[ic] = 2.*(S*theCosA-Q*theSinA);
//       }
//       else {
// 	rc[ic] = 2.*(S*theSinA+Q*theCosA);
//       }
//       IC[ic] = rc[ic]*FHalf(sqr(2.*mc[ic]/mh),1.)*theMw/mc[ic];
//     }
//     Energy pre = sqr(mh/theMw)*mh/1024./pow(Constants::pi,3)
//       *sqr(weakCoupling(sqr(mh))*sqr(electroMagneticCoupling(sqr(mh)))/4./Constants::pi);
//     cerr << "testing lighter sbottom" << ix << " " 
// 	 << pre*std::norm(Isb1)/GeV << "\n";
//     cerr << "testing heavier sbottom" << ix << " " 
// 	 << pre*std::norm(Istau2)/GeV << "\n";
//     cerr << "testing lighter stop" << ix << " " 
// 	 << pre*std::norm(Ist1)/GeV << "\n";
//     cerr << "testing heavier stop" << ix << " " 
// 	 << pre*std::norm(Ist2)/GeV << "\n";
//     cerr << "testing lighter stau" << ix << " " 
// 	 << pre*std::norm(Istau1)/GeV << "\n";
//     cerr << "testing heavier stau" << ix << " " 
// 	 << pre*std::norm(Isb2)/GeV << "\n";
//     cerr << "testing top " << ix << " " 
// 	 << pre*std::norm(Itop)/GeV << "\n";
//     cerr << "testing bottom " << ix << " " 
// 	 << pre*std::norm(Ibot)/GeV << "\n";
//     cerr << "testing tau " << ix << " " 
// 	 << pre*std::norm(Itau)/GeV << "\n";
//     cerr << "testing higgs " << ix << " " 
// 	 << pre*std::norm(Ih)/GeV << "\n";
//     cerr << "testing W " << ix << " " 
// 	 << pre*std::norm(IW)/GeV << "\n";
//     cerr << "testing chi1 " << ix << " " 
// 	 << pre*std::norm(IC[0])/GeV << "\n";
//     cerr << "testing chi2 " << ix << " " 
// 	 << pre*std::norm(IC[1])/GeV << "\n";
//     cerr << "testing chi " << ix << " " 
// 	 << pre*std::norm(IC[0]+IC[1])/GeV << "\n";
//     cerr << "testing higgs width " << ix << " " 
// 	 << pre*std::norm(Isb1+Isb2+Ist1+Ist2+Istau1+Istau2+
// 			  Itop+Ibot+Itau+Ih+IW+IC[0]+IC[1])/GeV << "\n";
//   }
  if(loopToolsInitialized()) Looptools::ltexi();
}
 
#line 1 "./SSNNPVertex.cc"
// -*- C++ -*-
//
// SSNNPVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSNNPVertex class.
//

#include "SSNNPVertex.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/Susy/MixingMatrix.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Looptools/clooptools.h"
#include "Herwig/Utilities/Maths.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSNNPVertex::SSNNPVertex() : _includeOnShell(false), _realIntegral(false), 
			     _sw(0.), _cw(0.), _id1last(0), 
			     _id2last(0), _q2last(ZERO), _couplast(0.),
			     _leftlast(ZERO), _rightlast(ZERO) {
  orderInGem(3);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void SSNNPVertex::doinit() {
  Looptools::ltini();
  int ineu[5] = {1000022,1000023,1000025,1000035,1000045};
  for(unsigned int i = 0; i < 5; ++i) {
    for(unsigned int j = 0; j < 5; ++j) {
      addToList(ineu[i], ineu[j], 22);
    }
  }
  GeneralFFVVertex::doinit();
  tMSSMPtr theSS = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if(!theSS)
    throw InitException() << "SSNNPVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  
  _theN  = theSS->neutralinoMix();
  _theU = theSS->charginoUMix();
  _theV = theSS->charginoVMix();
  if(!_theN || !_theU || ! _theV)
    throw InitException() << "SSNNPVertex::doinit - The neutralino "
			  << "mixing matrix pointer is null." 
			  << Exception::abortnow;
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt(1 - _sw*_sw);
  _mw = getParticleData(ParticleID::Wplus)->mass();
  double tb = theSS->tanBeta();
  _sb = tb/sqrt(1 + sqr(tb));
  _cb = sqrt(1 - sqr(_sb));
  _stop = theSS->stopMix();
  _sbot = theSS->sbottomMix();
  _stau = theSS->stauMix();
  Looptools::ltexi();
}

void SSNNPVertex::dofinish() {
  Looptools::ltexi();
  GeneralFFVVertex::dofinish();
}

void SSNNPVertex::doinitrun() {
  Looptools::ltini();
  GeneralFFVVertex::doinitrun();
}

void SSNNPVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _cw << _theN << ounit(_mw,GeV) << _sb << _cb
     << _stop << _sbot << _stau << _theU << _theV << _includeOnShell
     << _realIntegral;
}

void SSNNPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _theN >> iunit(_mw,GeV) >> _sb >> _cb
     >> _stop >> _sbot >> _stau >> _theU >> _theV >> _includeOnShell
     >> _realIntegral;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSNNPVertex,Helicity::GeneralFFVVertex>
describeSSNNPVertex("Herwig::SSNNPVertex", "HwSusy.so");

void SSNNPVertex::Init() {

  static ClassDocumentation<SSNNPVertex> documentation
    ("The loop-mediated coupling of the photon to a pair of neutralinos");

  static Switch<SSNNPVertex,bool> interfaceIncludeOnShellIntermediates
    ("IncludeOnShellIntermediates",
     "Whether or not to include on-shell intermediate states",
     &SSNNPVertex::_includeOnShell, false, false, false);
  static SwitchOption interfaceIncludeOnShellIntermediatesYes
    (interfaceIncludeOnShellIntermediates,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncludeOnShellIntermediatesNo
    (interfaceIncludeOnShellIntermediates,
     "No",
     "Don't incldue them",
     false);

  static Switch<SSNNPVertex,bool> interfaceRealIntegral
    ("RealIntegral",
     "Only include the real parts of the integrals",
     &SSNNPVertex::_realIntegral, false, false, false);
  static SwitchOption interfaceRealIntegralYes
    (interfaceRealIntegral,
     "Yes",
     "Only include the real part",
     true);
  static SwitchOption interfaceRealIntegralNo
    (interfaceRealIntegral,
     "No",
     "Don't include the real part",
     false);

}

void SSNNPVertex::setCoupling(Energy2 q2, tcPDPtr part1,
#ifndef NDEBUG
			      tcPDPtr part2,tcPDPtr part3) {
#else
			      tcPDPtr part2,tcPDPtr) {
#endif
  int o[2]={1,0};
  long in1 = part1->id();
  long in2 = part2->id();
  Energy Mj = part1->mass();
  Energy Mi = part2->mass();
  // checks of the particle ids
  assert(part3->id()==ParticleID::gamma);
  assert(in1 == ParticleID::SUSY_chi_10 || in1 == ParticleID::SUSY_chi_20 ||
	 in1 == ParticleID::SUSY_chi_30 || in1 == ParticleID::SUSY_chi_40 || 
	 in1 == 1000045                 );
  assert(in2 == ParticleID::SUSY_chi_10 || in2 == ParticleID::SUSY_chi_20 ||
	 in2 == ParticleID::SUSY_chi_30 || in2 == ParticleID::SUSY_chi_40 ||
	 in2 == 1000045              );
  // normal couplings are zero
  setLeft (0.);
  setRight(0.);
  if(in1==in2) {
    _leftlast  = ZERO;
    _rightlast = ZERO;
    setLeftSigma (_leftlast );
    setRightSigma(_rightlast);
    return;
  }
  if(q2 != _q2last || _couplast==0.) {
    _q2last = q2;
    _couplast = sqr(weakCoupling(q2))*
      electroMagneticCoupling(q2)/32./sqr(Constants::pi);
  }
  if(in1 != _id1last || in2 != _id2last) {
    _leftlast  = ZERO;
    _rightlast = ZERO;
    _id1last = in1;
    _id2last = in2;
    unsigned int neu1(in1 - 1000022), neu2(in2 - 1000022);
    if(neu1 > 1) neu1 = (in1-1000005)/10;
    if(neu2 > 1) neu2 = (in2-1000005)/10;
    Complex n1prime[2] = { (*_theN)(neu2,0)*_cw + (*_theN)(neu2,1)*_sw ,
			   (*_theN)(neu1,0)*_cw + (*_theN)(neu1,1)*_sw };
    Complex n2prime[2] = { (*_theN)(neu2,1)*_cw - (*_theN)(neu2,0)*_sw ,
			   (*_theN)(neu1,1)*_cw - (*_theN)(neu1,0)*_sw };
    // sfermion/fermion loops
    for(long iferm=1;iferm<16;++iferm) {
      if(iferm==7) iferm=11;
      if(iferm%2==0&&iferm>11) ++iferm;
      tcPDPtr smf = getParticleData(iferm);
      Energy mf = smf->mass();
      double qf = smf->charge()/eplus;
      double y = 0.5*mf/_mw;
      Complex bracketl[2] = { qf*_sw*( conj(n1prime[0]) - _sw*conj(n2prime[0])/_cw ) ,
			      qf*_sw*( conj(n1prime[1]) - _sw*conj(n2prime[1])/_cw ) };
      double lambda(0.);
      //neutralino mixing element
      Complex nlf[2]={0.,0.};
      if( iferm % 2 == 0 ) {
 	y /= _sb;
 	lambda = -0.5 + qf*sqr(_sw);
 	nlf[0] = (*_theN)(neu2,3);
 	nlf[1] = (*_theN)(neu1,3);
      }
      else { 
	y /= _cb;
	lambda = 0.5 + qf*sqr(_sw);
	nlf[0] = (*_theN)(neu2,2);
	nlf[1] = (*_theN)(neu1,2);
      }
      Complex bracketr[2] = { _sw*qf*n1prime[0] - n2prime[0]*lambda/_cw ,
			      _sw*qf*n1prime[1] - n2prime[1]*lambda/_cw };
      for(long iy=0;iy<2;++iy) {
	long isf = 1000000*(1+iy)+iferm;
	Energy msf = getParticleData(isf)->mass();
	if(!_includeOnShell&&(mf+msf<Mj||mf+msf<Mi)) continue;
	Complex g[2][2];
	Complex ma1(0.), ma2(0.);
	// heavy fermions
	if( iferm == 5 || iferm == 6 || iferm == 15 ) {
	  if( iferm == 5 ) {
	    ma1 = (*_sbot)(iy,0);
	    ma2 = (*_sbot)(iy,1); 
	  } 
	  else if( iferm == 6 ) {
	    ma1 = (*_stop)(iy,0);
	    ma2 = (*_stop)(iy,1);
	  } 
	  else {
	    ma1 = (*_stau)(iy,0);
	    ma2 = (*_stau)(iy,1);
	  }
	}
	else if(iy==0) {
	  ma1 = 1.;
	}
	else {
	  ma2 = 1.;
	}
	for(unsigned int ix=0;ix<2;++ix) {
	  g[ix][0] = y*conj(nlf[ix])*ma1 - ma2*bracketl[ix];
	  g[ix][1] = y*nlf[ix]*ma2 + ma1*bracketr[ix];
	}
	swap(g[0][0],g[0][1]);
	complex<InvEnergy2> I,J,K,I2;
	loopIntegrals(Mi,Mj,msf,mf,I,J,K,I2);
	complex<InvEnergy> coup[2];
	for(unsigned int ix=0;ix<2;++ix) {
	  coup[ix] = Mj*(I2-K)*(g[0][ix]*g[1][o[ix]]-conj(g[0][o[ix]]*g[1][ix]))
	    +Mi*K*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]))
	    +mf*I*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]));
	}
	double fact = 4.*qf;
	if(iferm<=6) fact *=3.;
	_leftlast  += fact*coup[0];
	_rightlast += fact*coup[1];
      }
    }
    // the chargino W contribution
    for(unsigned int ic=0;ic<2;++ic) {
      long id = ic==0 ? 
    	ParticleID::SUSY_chi_1plus : ParticleID::SUSY_chi_2plus;
      Energy Mk = getParticleData(id)->mass();
      if(!_includeOnShell&&(Mk+_mw<Mj||Mk+_mw<Mi)) continue;
      complex<InvEnergy2> I,J,K,I2;
      loopIntegrals(Mi,Mj,_mw,Mk,I,J,K,I2);
      Complex g[2][2];
      for(unsigned int ix=0;ix<2;++ix) {
    	unsigned int in = ix==0 ? neu2 : neu1;
    	g[ix][0] = 
    	  conj((*_theN)(in, 1))*(*_theV)(ic, 0) - 
    	  conj((*_theN)(in, 3))*(*_theV)(ic, 1)/sqrt(2);
    	g[ix][1] = 
    	  (*_theN)(in, 1)*conj((*_theU)(ic, 0)) +
    	  (*_theN)(in, 2)*conj((*_theU)(ic, 1))/sqrt(2);
      }
      complex<InvEnergy> coup[2];
      for(unsigned int ix=0;ix<2;++ix) {
    	coup[ix] = 
    	  Mj*(I2-J-K)*(g[0][o[ix]]*g[1][o[ix]]-conj(g[0][ix]*g[1][ix]))-
    	  Mi*(J-K)*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]))+
    	  2.*Mk*J*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]));
      }
      _leftlast  += 4.*coup[0];
      _rightlast += 4.*coup[1];
    }
    // the chargino charged higgs contribution
    Energy mh = getParticleData(ParticleID::Hplus)->mass();
    for(unsigned int ic=0;ic<2;++ic) {
      long id = ic==0 ? 
    	ParticleID::SUSY_chi_1plus : ParticleID::SUSY_chi_2plus;
      Energy Mk = getParticleData(id)->mass();
      if(!_includeOnShell&&(Mk+mh<Mj||Mk+mh<Mi)) continue;
      complex<InvEnergy2> I,J,K,I2;
      loopIntegrals(Mi,Mj,mh,Mk,I,J,K,I2);
      Complex g[2][2];
      for(unsigned int ix=0;ix<2;++ix) {
    	unsigned int in = ix==0 ? neu2 : neu1;
    	g[ix][0] =  (*_theN)(in, 3)*(*_theV)(ic,0) 
    	  +               ((*_theN)(in,1) + (*_theN)(in,0)*_sw/_cw)*
    	  (*_theV)(ic,1)/sqrt(2);
    	g[ix][0] *= _cb;
    	g[ix][1] = conj((*_theN)(in, 2)*(*_theU)(ic,0) 
    			- ((*_theN)(in,1) + (*_theN)(in,0)*_sw/_cw)*
    			(*_theU)(ic,1)/sqrt(2));
    	g[ix][1] *= _sb;
      }
      swap(g[1][0],g[1][1]);
      complex<InvEnergy> coup[2];
      for(unsigned int ix=0;ix<2;++ix) {
    	coup[ix] = Mj*(I2-K)*(g[0][ix]*g[1][o[ix]]-conj(g[0][o[ix]]*g[1][ix]))
    	  +Mi*K*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]))
    	  +Mk*I*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]));
      }
      _leftlast  += 2.*coup[0];
      _rightlast += 2.*coup[1];
    }
    // the chargino goldstone contribution
    for(unsigned int ic=0;ic<2;++ic) {
      long id = ic==0 ? 
    	ParticleID::SUSY_chi_1plus : ParticleID::SUSY_chi_2plus;
      Energy Mk = getParticleData(id)->mass();
      if(!_includeOnShell&&(Mk+_mw<Mj||Mk+_mw<Mi)) continue;
      complex<InvEnergy2> I,J,K,I2;
      loopIntegrals(Mi,Mj,_mw,Mk,I,J,K,I2);
      Complex g[2][2];
      for(unsigned int ix=0;ix<2;++ix) {
    	unsigned int in = ix==0 ? neu2 : neu1;
    	g[ix][0] = (*_theN)(in, 3)*(*_theV)(ic,0) 
    	  + ((*_theN)(in,1) + (*_theN)(in,0)*_sw/_cw)*
    	  (*_theV)(ic,1)/sqrt(2);
    	g[ix][0] *=-_sb;
    	g[ix][1] = conj((*_theN)(in, 2)*(*_theU)(ic,0) 
    			- ((*_theN)(in,1) + (*_theN)(in,0)*_sw/_cw)*
    			(*_theU)(ic,1)/sqrt(2));
    	g[ix][1] *= _cb;
      }
      swap(g[1][0],g[1][1]);
      complex<InvEnergy> coup[2];
      for(unsigned int ix=0;ix<2;++ix) {
    	coup[ix] = Mj*(I2-K)*(g[0][ix]*g[1][o[ix]]-conj(g[0][o[ix]]*g[1][ix]))
    	  +Mi*K*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]))
    	  +Mk*I*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]));
      }
      _leftlast  += 2.*coup[0];
      _rightlast += 2.*coup[1];
    }
  }
  norm(_couplast);
  setLeftSigma ( _leftlast);
  setRightSigma(_rightlast);  
}

void SSNNPVertex::loopIntegrals(Energy Mi, Energy Mj, Energy M, Energy m,
				complex<InvEnergy2> & I, complex<InvEnergy2> & J,
				complex<InvEnergy2> & K, complex<InvEnergy2> & I2) {
  Energy2 m2(sqr(m)),M2(sqr(M)),Mi2(sqr(Mi)),Mj2(sqr(Mj));
  double min2  = Mj2*UnitRemoval::InvE2;
  double mout2 = Mi2*UnitRemoval::InvE2;
  double mf2   = m2 *UnitRemoval::InvE2;
  double ms2   = M2 *UnitRemoval::InvE2;
  I  = Looptools::C0i(Looptools::cc0,min2,mout2,0.,mf2,ms2,mf2)*UnitRemoval::InvE2;
  J  = Looptools::C0i(Looptools::cc0,min2,mout2,0.,ms2,mf2,ms2)*UnitRemoval::InvE2;
  I2 =-Looptools::C0i(Looptools::cc1,min2,mout2,0.,mf2,ms2,mf2)*UnitRemoval::InvE2;
  K  = (1.+Complex(m2*I+M2*J-Mj2*I2))/(Mi2-Mj2);
  if(_realIntegral) {
    I  = I .real();
    J  = J .real();
    I2 = I2.real();
    K  = K .real();
  }
}
#line 1 "./SSGNGVertex.cc"
// -*- C++ -*-
//
// SSGNGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGNGVertex class.
//

#include "SSGNGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/Susy/MixingMatrix.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Looptools/clooptools.h"

using namespace ThePEG::Helicity;
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
}

SSGNGVertex::SSGNGVertex() : _includeOnShell(false), _realIntegral(false), 
			     _omitLightQuarkYukawas(false),
			     _sw(0.), _cw(0.), _idlast(0), 
			     _q2last(ZERO), _couplast(0.),
			     _leftlast(ZERO), _rightlast(ZERO),
			     _initLoops(false) {
  orderInGem(1);
  orderInGs(2);
  colourStructure(ColourStructure::DELTA);
}

void SSGNGVertex::doinit() {
  if(!_initLoops) {
    Looptools::ltini();
    _initLoops = true;
  }
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "SSGNGVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  _theN  = model->neutralinoMix();
  if(!_theN )
    throw InitException() << "SSGNGVertex::doinit - The neutralino "
			  << "mixing matrix pointer is null." 
			  << Exception::abortnow;
  vector<long> ineu(4);
  ineu[0] = 1000022; ineu[1] = 1000023; 
  ineu[2] = 1000025; ineu[3] = 1000035;
  if(_theN->size().first==5)
    ineu.push_back(1000045);
  else if(_theN->size().first==7) {
    if(model->majoranaNeutrinos()) {
      ineu.push_back(17);
      ineu.push_back(18);
      ineu.push_back(19);
    }
    else {
      ineu.push_back(12);
      ineu.push_back(14);
      ineu.push_back(16);
    }
  }
  for(unsigned int i = 0; i < ineu.size(); ++i) {
    addToList(1000021, ineu[i], 21);
  }
  GeneralFFVVertex::doinit();
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt(1 - _sw*_sw);
  _mw = getParticleData(ParticleID::Wplus)->mass();
  double tb = model->tanBeta();
  _sb = tb/sqrt(1 + sqr(tb));
  _cb = sqrt(1 - sqr(_sb));
  _stop = model->stopMix();
  _sbot = model->sbottomMix();
  Looptools::ltexi();
}

void SSGNGVertex::dofinish() {
  Looptools::ltexi();
  GeneralFFVVertex::dofinish();
}

void SSGNGVertex::doinitrun() {
  if(!_initLoops) {
    Looptools::ltini();
    _initLoops = true;
  }
  GeneralFFVVertex::doinitrun();
}

void SSGNGVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _cw << _theN << ounit(_mw,GeV) << _sb << _cb
     << _stop << _sbot << _includeOnShell << _omitLightQuarkYukawas
     << _realIntegral;
}

void SSGNGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _theN >> iunit(_mw,GeV) >> _sb >> _cb
     >> _stop >> _sbot >> _includeOnShell >> _omitLightQuarkYukawas
     >> _realIntegral;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSGNGVertex,GeneralFFVVertex>
describeHerwigSSGNGVertex("Herwig::SSGNGVertex", "HwSusy.so");

void SSGNGVertex::Init() {

  static ClassDocumentation<SSGNGVertex> documentation
    ("The loop-mediated coupling of the gluino to a gluon and a neutralino.");

  static Switch<SSGNGVertex,bool> interfaceIncludeOnShellIntermediates
    ("IncludeOnShellIntermediates",
     "Whether or not to include on-shell intermediate states",
     &SSGNGVertex::_includeOnShell, false, false, false);
  static SwitchOption interfaceIncludeOnShellIntermediatesYes
    (interfaceIncludeOnShellIntermediates,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncludeOnShellIntermediatesNo
    (interfaceIncludeOnShellIntermediates,
     "No",
     "Don't incldue them",
     false);

  static Switch<SSGNGVertex,bool> interfaceOmitLightQuarkYukawas
    ("OmitLightQuarkYukawas",
     "Omit the yukawa type couplings for down, up, strange"
     " and charm quarks, mainly for testing vs ISAJET",
     &SSGNGVertex::_omitLightQuarkYukawas, false, false, false);
  static SwitchOption interfaceOmitLightQuarkYukawasNo
    (interfaceOmitLightQuarkYukawas,
     "No",
     "Include the Yukawas",
     false);
  static SwitchOption interfaceOmitLightQuarkYukawasYes
    (interfaceOmitLightQuarkYukawas,
     "Yes",
     "Omit Yukawas",
     true);

  static Switch<SSGNGVertex,bool> interfaceRealIntegral
    ("RealIntegral",
     "Only include the real parts of the integrals",
     &SSGNGVertex::_realIntegral, false, false, false);
  static SwitchOption interfaceRealIntegralYes
    (interfaceRealIntegral,
     "Yes",
     "Only include the real part",
     true);
  static SwitchOption interfaceRealIntegralNo
    (interfaceRealIntegral,
     "No",
     "Don't include the real part",
     false);
}

void SSGNGVertex::setCoupling(Energy2 q2, tcPDPtr part1,
#ifndef NDEBUG
			      tcPDPtr part2,tcPDPtr part3) {
#else
                              tcPDPtr part2,tcPDPtr) {
#endif
  int o[2]={1,0};
  long in1 = part1->id();
  long in2 = part2->id();
  Energy Mj = part1->mass();
  Energy Mi = part2->mass();
  if(in1!=ParticleID::SUSY_g) {
    swap(in1,in2);
    swap(Mj,Mi);
  }
  // checks of the particle ids
  assert(part3->id()==ParticleID::g);
  assert(in1 == ParticleID::SUSY_g);
  assert(in2 == ParticleID::SUSY_chi_10 || in2 == ParticleID::SUSY_chi_20 ||
	 in2 == ParticleID::SUSY_chi_30 || in2 == ParticleID::SUSY_chi_40 ||
	 in2 == 1000045 || in2==12 || in2==14 || in2==16 || in2==17 || in2==18 || in2==19 );
  // normal couplings are zero
  setLeft (0.);
  setRight(0.);
  if(in2 != _idlast || q2 !=_q2last) {
    if(!_initLoops) {
      Looptools::ltini();
      _initLoops = true;
    }
    Looptools::clearcache();
    _leftlast  = ZERO;
    _rightlast = ZERO;
    _idlast = in2;
    unsigned int neu = neutralinoIndex(in2);
    Complex n1prime = (*_theN)(neu,0)*_cw + (*_theN)(neu,1)*_sw;
    Complex n2prime = (*_theN)(neu,1)*_cw - (*_theN)(neu,0)*_sw;
    // squark/quark loops
    for(long iferm=1;iferm<7;++iferm) {
      tcPDPtr smf = getParticleData(iferm);
      Energy mf = smf->mass();
      double qf = smf->charge()/eplus;
      double y = (!(iferm<=4&&_omitLightQuarkYukawas)) ?
	0.5*double(dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel())->mass(q2,smf)/_mw) : 0.;
      Complex bracketl = qf*_sw*( conj(n1prime) - _sw*conj(n2prime)/_cw );
      double lambda(0.);
      // neutralino mixing element
      Complex nlf(0.);
      if( iferm % 2 == 0 ) {
  	y /= _sb;
  	lambda = -0.5 + qf*sqr(_sw);
   	nlf = (*_theN)(neu,3);
      }
      else { 
  	y /= _cb;
  	lambda = 0.5 + qf*sqr(_sw);
  	nlf = (*_theN)(neu,2);
      }
      Complex bracketr = _sw*qf*n1prime - n2prime*lambda/_cw;
      for(long iy=0;iy<2;++iy) {
	long isf = 1000000*(1+iy)+iferm;
	Energy msf = getParticleData(isf)->mass();
	if(!_includeOnShell&&(mf+msf<Mj||mf+msf<Mi)) continue;
	Complex g[2][2];
	Complex ma1(0.), ma2(0.);
	// heavy fermions
	if( iferm == 5 || iferm == 6 ) {
	  if( iferm == 5 ) {
	    ma1 = (*_sbot)(iy,0);
	    ma2 = (*_sbot)(iy,1); 
	  } 
	  else if( iferm == 6 ) {
	    ma1 = (*_stop)(iy,0);
	    ma2 = (*_stop)(iy,1);
	  }
	}
	else if(iy==0) {
	  ma1 = 1.;
	}
	else {
	  ma2 = 1.;
	}
	g[0][0] = y*conj(nlf)*ma1 - ma2*bracketl;
	g[0][1] = y*nlf      *ma2 + ma1*bracketr;
	g[1][0] =  - ma2;
	g[1][1] =  + ma1;
 	swap(g[0][0],g[0][1]);
	complex<InvEnergy2> I,J,K,I2;
	loopIntegrals(Mi,Mj,msf,mf,I,J,K,I2);
	complex<InvEnergy> coup[2];
	for(unsigned int ix=0;ix<2;++ix) {
	  coup[ix] = Mj*(I2-K)*(g[0][ix]*g[1][o[ix]]-conj(g[0][o[ix]]*g[1][ix]))
	    +Mi*K*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]))
	    +mf*I*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]));
	}
	_leftlast  += 2.*coup[0];
	_rightlast += 2.*coup[1];
      }
    }
  }
  if(q2 != _q2last || _couplast==0.) {
    _q2last = q2;
    _couplast = weakCoupling(q2)*sqr(strongCoupling(q2))/
      32./sqr(Constants::pi);
  }
  norm(_couplast);
  setLeftSigma ( _leftlast);
  setRightSigma(_rightlast);
}

void SSGNGVertex::loopIntegrals(Energy Mi, Energy Mj, Energy M, Energy m,
				complex<InvEnergy2> & I, complex<InvEnergy2> & J,
				complex<InvEnergy2> & K, complex<InvEnergy2> & I2) {
  Energy2 m2(sqr(m)),M2(sqr(M)),Mi2(sqr(Mi)),Mj2(sqr(Mj));
  double min2  = Mj2*UnitRemoval::InvE2;
  double mout2 = Mi2*UnitRemoval::InvE2;
  double mf2   = m2 *UnitRemoval::InvE2;
  double ms2   = M2 *UnitRemoval::InvE2;
  I  = Looptools::C0i(Looptools::cc0,min2,mout2,0.,mf2,ms2,mf2)*UnitRemoval::InvE2;
  J  = Looptools::C0i(Looptools::cc0,min2,mout2,0.,ms2,mf2,ms2)*UnitRemoval::InvE2;
  I2 =-Looptools::C0i(Looptools::cc1,min2,mout2,0.,mf2,ms2,mf2)*UnitRemoval::InvE2;
  K  = (1.+Complex(m2*I+M2*J-Mj2*I2))/(Mi2-Mj2);
  if(_realIntegral) {
    I  = I .real();
    J  = J .real();
    I2 = I2.real();
    K  = K .real();
  }
}
#line 1 "./SSNCTVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSNCTVertex class.
//

#include "SSNCTVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"

using namespace Herwig;

SSNCTVertex::SSNCTVertex() : MX_(2.e16*GeV), 
			     sw_(0.), cw_(0.), mw_(ZERO), 
			     sb_(0.), cb_(0.), q2last_(), couplast_(0.),
			     leftlast_(0.), rightlast_(0.), idlast_(0),
			     epsilon_(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

IBPtr SSNCTVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SSNCTVertex::fullclone() const {
  return new_ptr(*this);
}

void SSNCTVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(MX_,GeV) << nmix_ << sw_ << cw_ << ounit(mw_,GeV)
     << sb_ << cb_ << ounit(q2last_,GeV2) << couplast_
     << leftlast_ << rightlast_ << idlast_ << epsilon_;
}

void SSNCTVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(MX_,GeV) >> nmix_ >> sw_ >> cw_ >> iunit(mw_,GeV)
     >> sb_ >> cb_ >> iunit(q2last_,GeV2) >> couplast_
     >> leftlast_ >> rightlast_ >> idlast_ >> epsilon_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSNCTVertex,Helicity::FFSVertex>
describeHerwigSSNCTVertex("Herwig::SSNCTVertex", "HwSusy.so");

void SSNCTVertex::Init() {

  static ClassDocumentation<SSNCTVertex> documentation
    ("The SSNCTVertex class implements the flavour violating"
     " coupling of the top quark to a charm quark and a neutralino");

  static Parameter<SSNCTVertex,Energy> interfaceMX
    ("MX",
     "Unification scale for the loop",
     &SSNCTVertex::MX_, GeV, 2.e16*GeV, 1.e14*GeV, 1.e20*GeV,
     false, false, Interface::limited);

}

void SSNCTVertex::doinit() {
  long neut[5] = {1000022, 1000023, 1000025, 1000035, 1000045};
  for(unsigned int nl = 0; nl < 5; ++nl) {
    addToList( neut[nl],  4, -1000006 );
    addToList( neut[nl], -4,  1000006 );
  }
  FFSVertex::doinit();
  // get the MSSM
  MSSMPtr model = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "SSNCTVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  // standard SUSY couplings
  // neutralino mixing
  nmix_ = model->neutralinoMix();
  if(!nmix_) throw InitException() << "SSNCTVertex::doinit() "
				   << "The neutralino mixing matrix pointer is null."
				   << Exception::abortnow;
  sw_ = sqrt(sin2ThetaW());
  mw_ = getParticleData(24)->mass();
  double tb = model->tanBeta();
  cw_ = sqrt(1. - sqr(sw_));
  sb_ = tb/sqrt(1 + sqr(tb));
  cb_ = sqrt(1 - sqr(sb_));
  // susy breaking scale
  Energy mSUSY = 
    sqrt(max(sqr(getParticleData(ParticleID::Z0)->mass()),
	     model->Mq3L()*model->MtR()));
  // CKM factor
  ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer 
    CKMptr = ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardCKM>::
    transient_const_pointer>(model->CKM());
  if(!CKMptr)
    throw Exception() << "Must have access to the Herwig::StandardCKM object"
		      << "for the CKM matrix in SSNCTVertex::doinit()"
		      << Exception::runerror;
  vector< vector<Complex > > CKM;
  CKM = CKMptr->getUnsquaredMatrix(generator()->standardModel()->families());
  // SM masses
  Energy mb = getParticleData(ParticleID::b)->mass();
  Energy mt = getParticleData(ParticleID::t)->mass();
  // squark masses
  Energy mt1 = getParticleData(1000006)->mass();
  Energy mcL = getParticleData(1000004)->mass();
  // mixing parameter
  Complex pre = sqr(weakCoupling(sqr(mSUSY)))/16./sqr(Constants::pi)*
    log(MX_/mSUSY)*sqr(double(mb/mw_))/sqr(cb_)*conj(CKM[2][2])*CKM[1][2];
  complex<Energy2> deltaL = -0.5*pre*(sqr(model->Mq2L())+sqr(model->Mq3L()) +
				      2.*model->Mh12()+2.*sqr(model->MbR()) +
				      2.*real(     model->bottomTrilinear()*
					      conj(model->bottomTrilinear())));
  complex<Energy2> deltaR =  pre*mt*conj(model->bottomTrilinear());
  if(abs(mt1-mcL)/abs(mt1+mcL)<1e-10) {
    epsilon_ = 0.;
  }
  else {
    epsilon_ = (deltaL*(*model->stopMix())(0,0)-
		deltaR*(*model->stopMix())(0,1))/(sqr(mt1)-sqr(mcL));
  }
}


void SSNCTVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr ) {
  long ism(abs(part1->id())), ineut(abs(part2->id()));
  tcPDPtr smfermion = part1;
  if( ism / 1000000 == 1 )  {
    swap( ism, ineut);
    smfermion = part2;
  }
  if(q2!= q2last_ || couplast_==0.) {
    couplast_ = -sqrt(2)*weakCoupling(q2);
    q2last_=q2;
  }
  norm(couplast_);
  if( ineut != idlast_) {
    idlast_ = ineut;
    // determine neutralino
    unsigned nl(0);
    switch( ineut ) {
    case 1000022 : nl = 0;
      break;
    case 1000023 : nl = 1;
      break;
    case 1000025 : nl = 2;
      break;
    case 1000035 : nl = 3;
      break;
    case 1000045 : nl = 4;
      break;
    default : assert(false);
    }
    // common primed neutralino matrices
    Complex n2prime = (*nmix_)(nl,1)*cw_ - (*nmix_)(nl,0)*sw_;
    Complex n1prime = (*nmix_)(nl,0)*cw_ + (*nmix_)(nl,1)*sw_;
    tcPDPtr smf = getParticleData(ism);
    double qf = smf->charge()/eplus;
    //Complex bracketl = qf*sw_*( conj(n1prime) - sw_*conj(n2prime)/cw_ );
    double lambda(0.);
    //neutralino mixing element
    Complex nlf(0.);
    lambda = -0.5 + qf*sqr(sw_);
    nlf = (*nmix_)(nl,3);
    Complex bracketr = sw_*qf*n1prime - n2prime*lambda/cw_;
    leftlast_  = 0.;
    rightlast_ = epsilon_*bracketr;
  }
  //determine the helicity order of the vertex
  if( smfermion->id() < 0 ) {
    left(conj(rightlast_));
    right(conj(leftlast_));
  }
  else {
    left(leftlast_);
    right(rightlast_);
  }
}
#line 1 "./SSGVNHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGVNHVertex class.
//

#include "SSGVNHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "MSSM.h"

using namespace Herwig;

SSGVNHVertex::SSGVNHVertex() : sa_(0.), sb_(0.), ca_(0.), cb_(0.),
			       MPlanck_(2.4e18*GeV) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr SSGVNHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SSGVNHVertex::fullclone() const {
  return new_ptr(*this);
}

void SSGVNHVertex::persistentOutput(PersistentOStream & os) const {
  os << sa_ << sb_ << ca_ << cb_ << nmix_ << ounit(MPlanck_,GeV);
}

void SSGVNHVertex::persistentInput(PersistentIStream & is, int) {
  is >> sa_ >> sb_ >> ca_ >> cb_ >> nmix_ >> iunit(MPlanck_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSGVNHVertex,Helicity::RFSVertex>
describeHerwigSSGVNHVertex("Herwig::SSGVNHVertex", "HwSusy.so");

void SSGVNHVertex::Init() {

  static ClassDocumentation<SSGVNHVertex> documentation
    ("The SSGVNHVertex class implments the coupling of the Higgs"
     " bosons to a gravitino and a neutralino");

}

void SSGVNHVertex::doinit() {
  long neu[4] = {ParticleID::SUSY_chi_10, ParticleID::SUSY_chi_20,
		 ParticleID::SUSY_chi_30, ParticleID::SUSY_chi_40};
  long higgs[3] =  {25, 35, 36};
  for(unsigned int i = 0; i < 3; ++i) {
    for(unsigned int j = 0; j < 4; ++j) {
      addToList(ParticleID::SUSY_Gravitino, neu[j], higgs[i]);
    }
  }
  RFSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !model )
    throw InitException() 
      << "SSGVNHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  double tanb = model->tanBeta();
  sb_ = tanb/sqrt(1. + sqr(tanb));
  cb_ = sqrt( 1. - sqr(sb_) );
  sa_ = sin(model->higgsMixingAngle());
  ca_ = sqrt(1. - sqr(sa_));
  nmix_ = model->neutralinoMix();
  MPlanck_ = model->MPlanck();
}

void SSGVNHVertex::setCoupling(Energy2 ,
#ifndef NDEBUG
			       tcPDPtr part1,
#else
			       tcPDPtr,
#endif
			       tcPDPtr part2,tcPDPtr part3) {
  assert(part1->id()==ParticleID::SUSY_Gravitino);
  assert(part3->iSpin()==PDT::Spin0);
  unsigned int neut = part2->id() - ParticleID::SUSY_chi_10;
  if(neut>1) neut = ( neut == 13 ) ? 3 : 2;
  int hid = part3->id();
  Complex coup;
  switch(hid) {
  case ParticleID::h0 :
    left (1.);
    right(1.);
    coup = -(*nmix_)(neut,2)*sa_+(*nmix_)(neut,3)*ca_;
    break;
  case ParticleID::H0 :
    left (1.);
    right(1.);
    coup =  (*nmix_)(neut,2)*ca_+(*nmix_)(neut,3)*sa_;
    break;
  case ParticleID::A0 :
    left (Complex(0.,-1.));
    right(Complex(0., 1.));
    coup =  (*nmix_)(neut,2)*sb_+(*nmix_)(neut,3)*cb_;
    break;
  default :
    assert(false);
  }
  norm(coup/MPlanck_*UnitRemoval::E);
}
#line 1 "./SSGVNVVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGVNVVertex class.
//

#include "SSGVNVVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "MSSM.h"

using namespace Herwig;

SSGVNVVertex::SSGVNVVertex() : sw_(0.), cw_(0.), sb_(0.), cb_(0.),
			       mz_(91.1876*GeV), MPlanck_(2.4e18*GeV) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr SSGVNVVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SSGVNVVertex::fullclone() const {
  return new_ptr(*this);
}

void SSGVNVVertex::persistentOutput(PersistentOStream & os) const {
  os << sw_ << cw_ << sb_ << cb_ << ounit(mz_,GeV) << nmix_ << ounit(MPlanck_,GeV);
}

void SSGVNVVertex::persistentInput(PersistentIStream & is, int) {
  is >> sw_ >> cw_ >> sb_ >> cb_ >> iunit(mz_,GeV) >> nmix_ >> iunit(MPlanck_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSGVNVVertex,Helicity::RFVVertex>
describeHerwigSSGVNVVertex("Herwig::SSGVNVVertex", "SSGVNVVertex.so");

void SSGVNVVertex::Init() {

  static ClassDocumentation<SSGVNVVertex> documentation
    ("The SSGVNVVertex class implements the coupling of the gravitino"
     " to the neutralino and a photon or Z boson, or the gluino and gluon.");

}

void SSGVNVVertex::doinit() {
  long neu[4] = {ParticleID::SUSY_chi_10, ParticleID::SUSY_chi_20,
		 ParticleID::SUSY_chi_30, ParticleID::SUSY_chi_40};
  for(unsigned int j = 0; j < 4; ++j) {
    addToList(ParticleID::SUSY_Gravitino, neu[j], ParticleID::gamma);
    addToList(ParticleID::SUSY_Gravitino, neu[j], ParticleID::Z0);
  }
  addToList(ParticleID::SUSY_Gravitino, ParticleID::SUSY_g, ParticleID::g);
  RFVVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !model )
    throw InitException() 
      << "SSGVNVVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  double tanb = model->tanBeta();
  sb_ = tanb/sqrt(1. + sqr(tanb));
  cb_ = sqrt( 1. - sqr(sb_) );
  sw_ = sqrt(sin2ThetaW());
  cw_ = sqrt(1. - sin2ThetaW());
  nmix_ = model->neutralinoMix();
  MPlanck_ = model->MPlanck();
}

void SSGVNVVertex::setCoupling(Energy2 ,
#ifndef NDEBUG
			       tcPDPtr part1,
#else
			       tcPDPtr,
#endif
			       tcPDPtr part2,tcPDPtr part3) {
  assert(part1->id()==ParticleID::SUSY_Gravitino);
  assert(part3->iSpin()==PDT::Spin1);
  unsigned int neut = part2->id() - ParticleID::SUSY_chi_10;
  if(neut>1) neut = ( neut == 13 ) ? 3 : 2;
  int bid = part3->id();
  Complex coup[2];
  vector<Complex> lV,rV;
  switch(bid) {
  case ParticleID::gamma :
    coup[0] = (*nmix_)(neut,0)*cw_+(*nmix_)(neut,1)*sw_;
    lV.push_back(-coup[0]*part2->mass()*UnitRemoval::InvE);
    lV.push_back( coup[0]);
    lV.push_back(0.);
    rV=lV;
    break;
  case ParticleID::Z0 :
    coup[0] = -(*nmix_)(neut,0)*sw_+(*nmix_)(neut,1)*cw_;
    coup[1] = -(*nmix_)(neut,2)*cb_+(*nmix_)(neut,3)*sb_;
    lV.push_back((-coup[0]*part2->mass()-mz_*coup[1])*UnitRemoval::InvE);
    lV.push_back( coup[0]);
    lV.push_back(0.);
    rV=lV;
    rV[0] = Complex((-coup[0]*part2->mass()+mz_*coup[1])*UnitRemoval::InvE);
    break;
  case ParticleID::g :
    lV.push_back(-double(part2->mass()*UnitRemoval::InvE));
    lV.push_back(1.);
    lV.push_back(0.);
    rV=lV;
    break;
  default :
    assert(false);
  }
  left (lV);
  right(rV);
  norm(double(1./MPlanck_*UnitRemoval::E));
}
#line 1 "./SSGVFSVertex.cc"
// -*- C++ -*-
//
// SSGVFSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGVFSVertex class.
//

#include "SSGVFSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGVFSVertex::SSGVFSVertex() : MPlanck_(2.4e18*GeV) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SSGVFSVertex::persistentOutput(PersistentOStream & os) const {
  os << stop_ << sbot_ << stau_ << ounit(MPlanck_,GeV);
}

void SSGVFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> stop_ >> sbot_ >> stau_ >> iunit(MPlanck_,GeV);
}

void SSGVFSVertex::doinit() {
  //quarks
  for(long ix=1;ix<7;++ix){
    addToList( ParticleID::SUSY_Gravitino,  ix, -(1000000+ix) );
    addToList( ParticleID::SUSY_Gravitino,  ix, -(2000000+ix) );
    addToList( ParticleID::SUSY_Gravitino, -ix,  (1000000+ix) );
    addToList( ParticleID::SUSY_Gravitino, -ix,  (2000000+ix) );
  }
  //leptons
  for(long ix=11;ix<17;++ix) {
    addToList( ParticleID::SUSY_Gravitino,  ix, -(1000000+ix) );
    addToList( ParticleID::SUSY_Gravitino, -ix,  (1000000+ix) );
    
    if( ix % 2 != 0 ) {
      addToList( ParticleID::SUSY_Gravitino,  ix, -(2000000+ix) );
      addToList( ParticleID::SUSY_Gravitino, -ix,  (2000000+ix) );
    }
  }
  RFSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !model )
    throw InitException() << "SSGVFSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  MPlanck_ = model->MPlanck();
  stop_ = model->stopMix();
  sbot_ = model->sbottomMix();
  stau_ = model->stauMix();
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSGVFSVertex,RFSVertex>
describeHerwigSSGVFSVertex("Herwig::SSGVFSVertex", "libHwSusy.so");

void SSGVFSVertex::Init() {

  static ClassDocumentation<SSGVFSVertex> documentation
    ("The SSGVFSVertex implements the coupling of  the gravitino to "
     "a fermion-sfermion");
}

void SSGVFSVertex::setCoupling(Energy2 ,
#ifndef NDEBUG
			       tcPDPtr part1,
#else
			       tcPDPtr,
#endif
			       tcPDPtr part2,tcPDPtr part3) {
  assert(part1->id()==ParticleID::SUSY_Gravitino);
  assert(part3->iSpin()==PDT::Spin0);
  norm(double(sqrt(2.)/MPlanck_*UnitRemoval::E));
  // sfermion mass eigenstate
  unsigned int alpha(abs(part3->id())/1000000 - 1);
  unsigned int ism(abs(part2->id()));
  Complex lc,rc;
  //heavy quarks/sleptons
  if( ism == 5 || ism == 6 || ism == 15 ) {
    Complex ma1(0.), ma2(0.);
    if( ism == 5 ) {
      ma1 = (*sbot_)(alpha,0);
      ma2 = (*sbot_)(alpha,1); 
    } 
    else if( ism == 6 ) {
      ma1 = (*stop_)(alpha,0);
      ma2 = (*stop_)(alpha,1);
    } 
    else {
      ma1 = (*stau_)(alpha,0);
      ma2 = (*stau_)(alpha,1);
    }
    lc = - ma2;
    rc = + ma1;
  }
  else {
    if( alpha == 0 ) {
      lc =  0.;
      rc =  1.;
    } 
    else {
      lc = -1.;
      rc =  0.;
    }
  }
  // determine the helicity order of the vertex
  if( part2->id() < 0 ) {
    left (conj(rc));
    right(conj(lc));
  }
  else {
    left (lc);
    right(rc);
  }
}
