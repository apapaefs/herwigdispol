// Validation-only plugin. Never enabled by production cards.
// Compile with -DHARD_SPIN_NEW for the new API and optionally
// -DHARD_SPIN_FIXTURES to exercise detached trees on fixed native hard inputs.
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Handlers/StepHandler.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/SpinInfo.h"
#include "ThePEG/EventRecord/StandardSelectors.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDF/PDFBase.h"
#include "ThePEG/PDF/LHAPDF6.h"
#include "Herwig/Shower/QTilde/QTildeShowerHandler.h"
#include "Herwig/Shower/QTilde/Kinematics/KinematicsReconstructor.h"
#include "Herwig/Shower/QTilde/Base/PartnerFinder.h"
#include "Herwig/Shower/QTilde/Base/HardTree.h"
#include "Herwig/Shower/QTilde/Base/ShowerTree.h"
#include "Herwig/Shower/QTilde/Base/ShowerVertex.h"
#include "Herwig/Shower/QTilde/SplittingFunctions/Sudakov1to2FormFactor.h"
#include "Herwig/Shower/QTilde/SplittingFunctions/HalfHalfOneSplitFn.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "fastjet/ClusterSequence.hh"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <execinfo.h>

using namespace ThePEG;
using namespace Herwig;
namespace {
std::string pendingHard, pendingLHE;
std::vector<std::string> lheEvents;
std::ofstream audit;
long fixtureEvents = 0;
std::map<std::string,std::pair<unsigned long,unsigned long>> polarizedPDFCalls;
void require(bool condition, const char* message) {
  if (!condition) throw Exception() << "Hard-spin validation: " << message << Exception::runerror;
}
bool physicalAverage(const RhoDMatrix& rho, long id) {
  const int n = rho.iSpin();
  for (int i=0; i<n; ++i) for (int j=0; j<n; ++j) {
    const double target = i==j && !(id==21 && i==1) ? 0.5 : 0.;
    if (std::abs(rho(i,j)-target)>1.e-10) return false;
  }
  return true;
}
#ifdef HARD_SPIN_FIXTURES
class FlatPDF: public PDFBase {
public:
  mutable unsigned calls = 0;
  bool canHandleParticle(tcPDPtr) const { return true; }
  cPDVector partons(tcPDPtr p) const { return cPDVector(1,p); }
  double xfl(tcPDPtr, tcPDPtr, Energy2, double, Energy2 = ZERO) const { ++calls; return 1.; }
  double xfx(tcPDPtr, tcPDPtr, Energy2, double, double = 0., Energy2 = ZERO) const { ++calls; return 1.; }
  IBPtr clone() const { return new_ptr(*this); }
  IBPtr fullclone() const { return new_ptr(*this); }
};
class ProbeSudakov: public HalfHalfOneSplitFn {
public:
  using SudakovFormFactor::setPDF;
  using SudakovFormFactor::PDFVetoRatio;
};
class ProbeShower: public QTildeShowerHandler {
public:
  ProbeShower() { setCurrentHandler(); }
  ~ProbeShower() { unSetCurrentHandler(); }
};
void fixedInputTests(tSubProPtr sub, const tPPair& beams) {
  PerturbativeProcessPtr process = new_ptr(PerturbativeProcess());
  tPVector originals;
  originals.push_back(sub->incoming().first);
  originals.push_back(sub->incoming().second);
  originals.insert(originals.end(),sub->outgoing().begin(),sub->outgoing().end());
  std::vector<SpinPtr> saved;
  for (unsigned i=0; i<originals.size(); ++i) {
    saved.push_back(originals[i]->spinInfo());
    if (i<2) process->incoming().push_back(std::make_pair(originals[i],PerturbativeProcessPtr()));
    else process->outgoing().push_back(std::make_pair(originals[i],PerturbativeProcessPtr()));
  }
  ShowerTreePtr tree = new_ptr(ShowerTree(process,false));
  const auto particles = tree->extractProgenitorParticles();
  for (auto part: particles) {
    require(!part->spinInfo(),"detached progenitor inherited SpinInfo");
    ShowerBasisPtr basis = new_ptr(ShowerBasis());
    Lorentz5Momentum p = part->momentum();
    Lorentz5Momentum n(p);
    n.setX(-p.x()); n.setY(-p.y()); n.setZ(-p.z());
    if (!part->isFinalState()) {
      p = part->momentum().z()>ZERO ? beams.first->momentum() : beams.second->momentum();
      n = part->momentum().z()>ZERO ? beams.second->momentum() : beams.first->momentum();
    }
    basis->setBasis(p,n,ShowerBasis::BackToBack);
    part->showerBasis(basis,false);
    const RhoDMatrix rho = part->extractRhoMatrix(part->isFinalState());
    require(physicalAverage(rho,part->id()),"seed not a physical spin average");
    require(!part->spinInfo()->productionVertex(),"seed linked to hard tensor");
  }
  tree->clear();
  for (auto part: tree->extractProgenitorParticles())
    require(!part->spinInfo(),"retry retained a rejected spin vertex");
  for (unsigned i=0; i<originals.size(); ++i)
    require(saved[i]==originals[i]->spinInfo(),"hard original SpinInfo changed");
  auto ordinary = new_ptr(FlatPDF());
  auto difference = new_ptr(FlatPDF());
  require(!ShowerHandler::currentHandlerIsSet(),"fixture entered during a live shower");
  ProbeShower context; // PDFVetoRatio uses the current handler's scale multiplier.
  ProbeSudakov sudakov;
  for (auto beam: {beams.first,beams.second}) {
    const tcBeamPtr beamdata = dynamic_ptr_cast<tcBeamPtr>(beam->dataPtr());
    require(bool(beamdata),"missing beam particle data");
    for (auto part: particles) for (double q2: {0.1,1.,100.}) {
      // Install nonzero polarization and nonnull difference PDFs deliberately:
      // isPol=false must be sufficient, including after a prior true setting.
      sudakov.setPDF(ordinary,2.5*GeV,true,ThreeVector<double>(0.2,0.1,1.),difference,difference);
      sudakov.setPDF(ordinary,2.5*GeV,false,ThreeVector<double>(0.2,0.1,1.),difference,difference);
      const auto h = sudakov.calculateHMatrix(beamdata,part->dataPtr(),0.1,q2*GeV2);
      require(physicalAverage(h,part->id()),"terminal H is polarized with isPol=false");
      const double ratio = sudakov.PDFVetoRatio(q2*GeV2,0.1,0.5,part->dataPtr(),
                                               part->dataPtr(),beamdata,h,1.);
      require(std::isfinite(ratio) && ratio>=0.,"invalid unpolarized PDF veto ratio");
    }
  }
  require(difference->calls==0,"difference PDF used by unpolarized Sudakov");
  ++fixtureEvents;
}
#endif
}

namespace Herwig {
class AuditPolarizedPDF: public ThePEG::LHAPDF {
public:
  double xfx(tcPDPtr beam,tcPDPtr parton,Energy2 scale,double x,
             double eps=0.,Energy2 virtuality=ZERO) const {
    auto& calls=polarizedPDFCalls[fullName()];
    bool inShower=ShowerHandler::currentHandlerIsSet();
#ifdef HARD_SPIN_NEW
    if (inShower) {
      auto shower=dynamic_ptr_cast<Ptr<QTildeShowerHandler>::transient_pointer>(ShowerHandler::currentHandler());
      if (shower && !shower->hardProcessSpin()) {
        // Herwig's global current-handler pointer can survive a rejected event.
        // Distinguish the next hard ME's PDF evaluation from a shower PDF call.
        void* frames[32];
        const int size=backtrace(frames,32);
        char** symbols=backtrace_symbols(frames,size);
        bool hardExtractor=false, showerKernel=false;
        for (int i=0; i<size && symbols; ++i) {
          const std::string symbol(symbols[i]);
          hardExtractor |= symbol.find("PolarizedPartonExtractor") != std::string::npos;
          showerKernel |= symbol.find("Sudakov") != std::string::npos ||
                          symbol.find("SplittingGenerator") != std::string::npos;
        }
        free(symbols);
        require(hardExtractor && !showerKernel,"off shower called a polarized PDF");
        inShower=false;
      }
    }
#endif
    if (inShower) ++calls.second;
    else ++calls.first;
    return ThePEG::LHAPDF::xfx(beam,parton,scale,x,eps,virtuality);
  }
  IBPtr clone() const { return new_ptr(*this); }
  IBPtr fullclone() const { return new_ptr(*this); }
  static void Init() { static ClassDocumentation<AuditPolarizedPDF> doc("Validation-only polarized-PDF call counter."); }
};
class HardSpinCapture: public StepHandler {
public:
  void handle(EventHandler& eh, const tPVector&, const Hint&) {
    auto sub = eh.currentCollision()->primarySubProcess();
    require(sub->outgoing().size()==2,"export requires a 2-to-2 hard event");
    tPVector legs{sub->incoming().first,sub->incoming().second};
    legs.insert(legs.end(),sub->outgoing().begin(),sub->outgoing().end());
    std::map<tcColinePtr,int> colours;
    auto colour = [&](tcColinePtr line) {
      if (!line) return 0;
      if (!colours.count(line)) colours[line]=501+colours.size();
      return colours[line];
    };
    std::ostringstream hard, lhe;
    hard << std::setprecision(17) << "[";
    lhe << std::setprecision(17);
    const double scale = sqrt(eh.lastXCombPtr()->lastShowerScale())/GeV;
    lhe << "<event>\n4 1 1 " << scale << " -1 -1\n";
    for (unsigned i=0; i<legs.size(); ++i) {
      const auto p = legs[i]->momentum();
      const int c = colour(legs[i]->colourLine()), ac = colour(legs[i]->antiColourLine());
      if(i) hard << ",";
      hard << "[" << legs[i]->id() << "," << c << "," << ac << ","
           << p.x()/GeV << "," << p.y()/GeV << "," << p.z()/GeV << ","
           << p.e()/GeV << "," << p.mass()/GeV << "]";
      lhe << legs[i]->id() << " " << (i<2 ? -1 : 1) << " "
          << (i<2 ? 0 : 1) << " " << (i<2 ? 0 : 2) << " " << c << " " << ac << " "
          << p.x()/GeV << " " << p.y()/GeV << " " << p.z()/GeV << " "
          << p.e()/GeV << " " << p.mass()/GeV << " 0 9\n";
    }
    hard << "]";
    lhe << "</event>\n";
    pendingHard=hard.str(); pendingLHE=lhe.str();
#ifdef HARD_SPIN_FIXTURES
    fixedInputTests(sub,eh.currentEvent()->incoming());
#endif
  }
  IBPtr clone() const { return new_ptr(*this); }
  IBPtr fullclone() const { return new_ptr(*this); }
  static void Init() { static ClassDocumentation<HardSpinCapture> doc("Validation-only pre-shower hard-input capture."); }
};

class HardSpinAudit: public AnalysisHandler {
public:
  void doinitrun() {
    AnalysisHandler::doinitrun();
    lheEvents.clear(); pendingHard.clear(); fixtureEvents=0; polarizedPDFCalls.clear();
    audit.open(generator()->filename()+".hard-spin.jsonl");
    audit << std::setprecision(17);
  }
  void analyze(tEventPtr event,long,int loop,int state) {
    if (!event || loop>0 || state!=0) return;
    require(!pendingHard.empty(),"missing pre-shower capture");
    require(std::abs(event->weight()-1.)<1.e-10,"export only supports unit positive hard weights");
    lheEvents.push_back(pendingLHE);
    tPVector all;
    event->select(std::back_inserter(all),SelectAll());
    unsigned hardLinks=0, showerLinks=0;
    for(auto p: all) {
      auto sp=dynamic_ptr_cast<ShowerParticlePtr>(p);
      if(!sp || !sp->spinInfo()) continue;
      auto spin=sp->spinInfo();
      if(dynamic_ptr_cast<tcHardVertexPtr>(spin->productionVertex())) ++hardLinks;
      if(dynamic_ptr_cast<tcHardVertexPtr>(spin->decayVertex())) ++hardLinks;
      if(dynamic_ptr_cast<Ptr<ShowerVertex>::transient_const_pointer>(spin->productionVertex())) ++showerLinks;
      if(dynamic_ptr_cast<Ptr<ShowerVertex>::transient_const_pointer>(spin->decayVertex())) ++showerLinks;
    }
    auto shower=dynamic_ptr_cast<Ptr<QTildeShowerHandler>::transient_pointer>(generator()->eventHandler()->cascadeHandler());
    bool hardSpin=true;
#ifdef HARD_SPIN_NEW
    hardSpin=shower->hardProcessSpin();
    if(!hardSpin) require(hardLinks==0,"shower connected to hard spin tensor in off mode");
#endif
    std::vector<fastjet::PseudoJet> inputs;
    std::vector<std::string> finalState;
    for(auto p: event->getFinalState()) {
#ifndef HARD_SPIN_COMPACT
      std::ostringstream record; record << std::setprecision(17) << p->id() << ":"
        << p->momentum().x()/GeV << ":" << p->momentum().y()/GeV << ":"
        << p->momentum().z()/GeV << ":" << p->momentum().e()/GeV;
      finalState.push_back(record.str());
#endif
      const long id=std::abs(p->id());
      if(id==12 || id==14 || id==16 || id==18) continue;
      inputs.emplace_back(p->momentum().x()/GeV,p->momentum().y()/GeV,
                          p->momentum().z()/GeV,p->momentum().e()/GeV);
    }
    std::sort(finalState.begin(),finalState.end());
    fastjet::ClusterSequence cluster(inputs,fastjet::JetDefinition(fastjet::antikt_algorithm,0.5));
    std::vector<fastjet::PseudoJet> jets;
    for(auto jet: fastjet::sorted_by_pt(cluster.inclusive_jets(2.)))
      if(std::abs(jet.eta())<1.5) jets.push_back(jet);
    const bool baseline=jets.size()>=2 && jets[0].pt()>5. && jets[1].pt()>4.;
    std::vector<double> features{double(baseline)};
    for(double threshold: {2.,3.,4.,5.,6.,8.,10.}) {
      features.push_back(baseline && jets.size()>=3 && jets[2].pt()>threshold ? 1. : 0.);
      features.push_back(baseline && jets.size()>=4 && jets[3].pt()>threshold ? 1. : 0.);
    }
    for(unsigned j: {2u,3u}) for(double edge: {2.,3.,4.,5.,6.,8.,10.,15.,20.,30.})
      features.push_back(baseline && jets.size()>j && jets[j].pt()>edge ? 1. : 0.);
    audit << "{\"sequence\":" << lheEvents.size() << ",\"hard\":" << pendingHard
          << ",\"weight\":" << event->weight() << ",\"hard_process_spin\":" << (hardSpin?"true":"false")
          << ",\"spin_correlations\":" << shower->spinCorrelations()
          << ",\"hard_links\":" << hardLinks << ",\"shower_links\":" << showerLinks << ",\"features\":[";
    for(unsigned i=0;i<features.size();++i) { if(i) audit << ","; audit << features[i]; }
    audit << "],\"final_state\":[";
    for(unsigned i=0;i<finalState.size();++i) { if(i) audit << ","; audit << "\"" << finalState[i] << "\""; }
    audit << "]}\n";
    pendingHard.clear();
  }
  void dofinish() {
    AnalysisHandler::dofinish();
    audit.close();
    const double sigma=generator()->integratedXSec()/picobarn;
#ifndef HARD_SPIN_DIS
    std::ofstream lhe(generator()->filename()+".hard.lhe");
    lhe << std::setprecision(17) << "<LesHouchesEvents version=\"1.0\">\n<init>\n"
        << "2212 2212 255 255 0 0 0 0 3 1\n" << sigma << " 0 1 1\n</init>\n";
    for(auto record: lheEvents) lhe << record;
    lhe << "</LesHouchesEvents>\n";
#endif
    std::ofstream summary(generator()->filename()+".hard-spin-summary.json");
    summary << std::setprecision(17) << "{\"events\":" << lheEvents.size()
            << ",\"sigma_pb\":" << sigma << ",\"fixed_input_tests\":" << fixtureEvents
            << ",\"polarized_pdf_calls\":{";
    bool first=true;
    for (const auto& item: polarizedPDFCalls) {
      if(!first) summary << ",";
      first=false;
      summary << "\"" << item.first << "\":{\"hard\":" << item.second.first
              << ",\"shower\":" << item.second.second << "}";
    }
    summary << "},\"next_random\":" << UseRandom::rnd() << "}\n";
  }
  IBPtr clone() const { return new_ptr(*this); }
  IBPtr fullclone() const { return new_ptr(*this); }
  static void Init() { static ClassDocumentation<HardSpinAudit> doc("Validation-only hard-spin audit and LHE export; 510 GeV unweighted pp."); }
};
DescribeNoPIOClass<HardSpinCapture,StepHandler> describeCapture("Herwig::HardSpinCapture","HardSpinValidation.so");
DescribeNoPIOClass<HardSpinAudit,AnalysisHandler> describeAudit("Herwig::HardSpinAudit","HardSpinValidation.so");
DescribeNoPIOClass<AuditPolarizedPDF,ThePEG::LHAPDF> describePDF("Herwig::AuditPolarizedPDF","HardSpinValidation.so");
}
