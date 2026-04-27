// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Tools/RivetHepMC.hh"
#include "HepMC3/Attribute.h"
#include "fastjet/ClusterSequence.hh"
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <valarray>

namespace Rivet {

  /// @brief DIS dijets in the Breit frame (cleaned per your specs)
  class MC_DIS_BREIT : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(MC_DIS_BREIT);

    struct ClusteredJet {
      FourMomentum breitMom;
      FourMomentum labMom;
    };

    struct JetInputParticle {
      FourMomentum breitMom;
      FourMomentum labMom;
      int pid;
    };

    struct MEPartonDebug {
      size_t taggedAny = 0;
      size_t taggedAnyRole2 = 0;
      size_t taggedAnyRole3 = 0;
      size_t taggedOutgoing = 0;
      size_t taggedRole2 = 0;
      size_t taggedRole3 = 0;
      size_t graphHardPartons = 0;
      size_t graphRootSeeds = 0;
      size_t graphTerminalPartons = 0;
      size_t terminalFallback = 0;
      size_t inputPartons = 0;
      int source = 0; // 0 none, 1 POWHEG tags, 2 hard-line fallback
    };

    struct MEPartonDebugTotals {
      size_t events = 0;
      size_t tagSource = 0;
      size_t fallbackSource = 0;
      size_t noSource = 0;
      size_t inputGE2 = 0;
      size_t jetsGE2 = 0;
      size_t passPt = 0;
      size_t passRapidity = 0;
      size_t taggedAny = 0;
      size_t taggedAnyRole2 = 0;
      size_t taggedAnyRole3 = 0;
      size_t graphHardPartons = 0;
      size_t graphRootSeeds = 0;
      size_t graphTerminalPartons = 0;
    };

    struct AnalysisWeights {
      double sigma = 1.0;
      double deltaLL = 0.0;
      bool haveRivetWeights = false;
    };

    enum class JetInputMode {
      FULL,
      HARDPARTONS,
      TOP2PARTONS,
      TOP3PARTONS,
      MEPARTONS
    };

    void init() {
      const string jetinputopt = toUpper(getOption("JETINPUT", "FULL"));
      if (jetinputopt == "FULL") {
        _jetInputMode = JetInputMode::FULL;
      } else if (jetinputopt == "HARDPARTONS") {
        _jetInputMode = JetInputMode::HARDPARTONS;
      } else if (jetinputopt == "TOP2PARTONS") {
        _jetInputMode = JetInputMode::TOP2PARTONS;
      } else if (jetinputopt == "TOP3PARTONS") {
        _jetInputMode = JetInputMode::TOP3PARTONS;
      } else if (jetinputopt == "MEPARTONS") {
        _jetInputMode = JetInputMode::MEPARTONS;
      } else {
        MSG_WARNING("Unknown JETINPUT option " + jetinputopt + ". Defaulting to FULL.");
        _jetInputMode = JetInputMode::FULL;
      }
      const string rivetweightsopt = toUpper(getOption("RIVETWEIGHTS", "NO"));
      _rivetWeightsMode = (rivetweightsopt == "YES" || rivetweightsopt == "ON" ||
                           rivetweightsopt == "TRUE" || rivetweightsopt == "1");
      if (_rivetWeightsMode) {
        MSG_INFO("RIVETWEIGHTS mode enabled: using HERWIG_DIS_SIGMA and HERWIG_DIS_DELTA_LL optional weights.");
      }

      const string dismodeopt = toUpper(getOption("DISMODE", "NC"));
      DISMode dismode = DISMode::L2L;
      if (dismodeopt == "NC") {
        dismode = DISMode::L2L;
      } else if (dismodeopt == "CC") {
        dismode = DISMode::L2NU;
      } else {
        MSG_WARNING("Unknown DISMODE option " + dismodeopt + ". Defaulting to NC.");
      }

      // DIS kinematics (gives Q2, x, y, and the Breit/Lab boosts)
      const DISLepton dislep(Cuts::OPEN, LeptonReco::ALL, ObjOrdering::ENERGY,
                             0.0, 0.0, 0.0, dismode);
      const DISKinematics diskin(dislep);
      declare(diskin, "Kinematics");
      declare(FinalState(), "LabFS");

      // Histograms (no Q2-binned sets)
     
      book(_h_Q2, "Q2", 100, 100.0, 2500.0);
      book(_h_Pt, "Pt", 15, 5.0, 30.0);
      book(_h_XBj, "XBj", 20, 0.0, 1.0);
      book(_h_Mjj, "Mjj", logspace(15, 10.0, 100.0));
      book(_h_Eta, "Eta", 15, 0.0, 2.5);
      book(_h_Zeta, "Zeta", 12, -1.75, -0.25);
      book(_h_pT1, "pT1", 15, 5.0, 30.0);
      book(_h_pT2, "pT2", 15, 5.0, 30.0);
      book(_h_pT2OverpT1, "pT2OverpT1", 15, 0.0, 1.0);
      book(_h_pTAsym, "pTAsym", 15, 0.0, 1.0);
      book(_h_Q2PreCut, "Q2PreCut", 100, 100.0, 2500.0);
      book(_h_XBjPreCut, "XBjPreCut", 20, 0.0, 1.0);
      book(_h_YPreCut, "YPreCut", 40, 0.2, 0.6);
      book(_h_pT1PreCut, "pT1PreCut", 30, 0.0, 30.0);
      book(_h_pT2PreCut, "pT2PreCut", 30, 0.0, 30.0);
      if (_rivetWeightsMode) {
        book(_h_DQ2, "DQ2", 100, 100.0, 2500.0);
        book(_h_DPt, "DPt", 15, 5.0, 30.0);
        book(_h_DXBj, "DXBj", 20, 0.0, 1.0);
        book(_h_DMjj, "DMjj", logspace(15, 10.0, 100.0));
        book(_h_DEta, "DEta", 15, 0.0, 2.5);
        book(_h_DZeta, "DZeta", 12, -1.75, -0.25);
        book(_h_DpT1, "DpT1", 15, 5.0, 30.0);
        book(_h_DpT2, "DpT2", 15, 5.0, 30.0);
        book(_h_DpT2OverpT1, "DpT2OverpT1", 15, 0.0, 1.0);
        book(_h_DpTAsym, "DpTAsym", 15, 0.0, 1.0);
        book(_h_DQ2PreCut, "DQ2PreCut", 100, 100.0, 2500.0);
        book(_h_DXBjPreCut, "DXBjPreCut", 20, 0.0, 1.0);
        book(_h_DYPreCut, "DYPreCut", 40, 0.2, 0.6);
        book(_h_DpT1PreCut, "DpT1PreCut", 30, 0.0, 30.0);
        book(_h_DpT2PreCut, "DpT2PreCut", 30, 0.0, 30.0);
      }
      book(_h_METaggedOutgoingCount, "METaggedOutgoingCount", 8, -0.5, 7.5);
      book(_h_METaggedRole2Count, "METaggedRole2Count", 5, -0.5, 4.5);
      book(_h_METaggedRole3Count, "METaggedRole3Count", 5, -0.5, 4.5);
      book(_h_METerminalFallbackCount, "METerminalFallbackCount", 8, -0.5, 7.5);
      book(_h_MEInputPartonCount, "MEInputPartonCount", 8, -0.5, 7.5);
      book(_h_MEJetCountPreCut, "MEJetCountPreCut", 8, -0.5, 7.5);
      book(_h_MEInputSource, "MEInputSource", 4, -0.5, 3.5);
      book(_h_MESelectionStage, "MESelectionStage", 5, -0.5, 4.5);
    }

    void analyze(const Event& event) {
      const DISKinematics& dis = apply<DISKinematics>(event, "Kinematics");
      const AnalysisWeights analysisWeights = eventAnalysisWeights(event);
      const double weight = analysisWeights.sigma;
      const double deltaWeight = analysisWeights.deltaLL;

      const double Q2  = dis.Q2();
      const double xbj = dis.x();
      const double y   = dis.y();

      // Kinematic window: 100 < Q^2 < 2500 GeV^2 and 0.2 < y < 0.6
      if (!inRange(Q2, 100*GeV2, 2500*GeV2)) vetoEvent;
      if (!inRange(y,  0.2,     0.6))       vetoEvent;

      vector<ClusteredJet> jets;
      if (_jetInputMode == JetInputMode::HARDPARTONS) {
        jets = clusterHardPartonJets(event);
      } else if (_jetInputMode == JetInputMode::TOP2PARTONS) {
        jets = clusterTopHardPartonJets(event, 2);
      } else if (_jetInputMode == JetInputMode::TOP3PARTONS) {
        jets = clusterTopHardPartonJets(event, 3);
      } else if (_jetInputMode == JetInputMode::MEPARTONS) {
        jets = clusterMEPartonJets(event);
      } else {
        jets = clusterFullFinalStateJets(event);
      }

      // Match the POLDIS "pre-cut" logic: fill after the DIS cuts even when
      // fewer than two jets survive the clustering, with missing jet pT = 0.
      const double jet1PreCutPt = jets.size() > 0 ? jets[0].breitMom.pT() : 0.0;
      const double jet2PreCutPt = jets.size() > 1 ? jets[1].breitMom.pT() : 0.0;
      _h_Q2PreCut->fill(Q2, weight);
      _h_XBjPreCut->fill(xbj, weight);
      _h_YPreCut->fill(y, weight);
      _h_pT1PreCut->fill(jet1PreCutPt, weight);
      _h_pT2PreCut->fill(jet2PreCutPt, weight);
      if (_rivetWeightsMode) {
        _h_DQ2PreCut->fill(Q2, deltaWeight);
        _h_DXBjPreCut->fill(xbj, deltaWeight);
        _h_DYPreCut->fill(y, deltaWeight);
        _h_DpT1PreCut->fill(jet1PreCutPt, deltaWeight);
        _h_DpT2PreCut->fill(jet2PreCutPt, deltaWeight);
      }
      if (_jetInputMode == JetInputMode::MEPARTONS) {
        fillMEPartonDiagnostics(jets, weight);
      }

      if (jets.size() < 2) vetoEvent;

      const FourMomentum& j1Bmom = jets[0].breitMom;
      const FourMomentum& j2Bmom = jets[1].breitMom;
      const FourMomentum& j1LabMom = jets[0].labMom;
      const FourMomentum& j2LabMom = jets[1].labMom;

      // Breit-frame pT thresholds: 5 GeV (lead), 4 GeV (sublead)
      if (j1Bmom.pT() < 5*GeV || j2Bmom.pT() < 4*GeV) vetoEvent;

      // POLDIS applies the lab-frame acceptance to the same two leading Breit jets.
      if (!inRange(j1LabMom.rapidity(), -3.5, 3.5)) vetoEvent;
      if (!inRange(j2LabMom.rapidity(), -3.5, 3.5)) vetoEvent;

      // Dijet mean transverse momentum (Breit)
      const double dijetPt = 0.5*(j1Bmom.pT() + j2Bmom.pT());
      const double pT2OverpT1 = j2Bmom.pT()/j1Bmom.pT();
      const double pTAsym = (j1Bmom.pT() - j2Bmom.pT())/(j1Bmom.pT() + j2Bmom.pT());

      // Invariant dijet mass in Breit.
      const double Mjj = (j1Bmom + j2Bmom).mass();

      // Eta* stays defined from the Breit-frame pseudorapidities.
      const int orientation = dis.orientation();
      const double eta1S = orientation * j1Bmom.eta();
      const double eta2S = orientation * j2Bmom.eta();
      const double etastar = 0.5 * std::abs(eta1S - eta2S);

      // log10 zeta = log10( x_Bj * (1 + Mjj^2 / Q^2) )
      const double logZeta = std::log10(xbj * (1.0 + sqr(Mjj)/Q2));

      // Fill histograms
      
      _h_Q2  ->fill(Q2, weight);
      _h_pT1->fill(j1Bmom.pT(), weight);
      _h_pT2->fill(j2Bmom.pT(), weight);
      _h_XBj ->fill(xbj, weight);
      _h_Pt  ->fill(dijetPt, weight);
      _h_Mjj ->fill(Mjj, weight);
      _h_Eta ->fill(etastar, weight);
      _h_Zeta->fill(logZeta, weight);
      _h_pT2OverpT1->fill(pT2OverpT1, weight);
      _h_pTAsym->fill(pTAsym, weight);
      if (_rivetWeightsMode) {
        _h_DQ2  ->fill(Q2, deltaWeight);
        _h_DpT1->fill(j1Bmom.pT(), deltaWeight);
        _h_DpT2->fill(j2Bmom.pT(), deltaWeight);
        _h_DXBj ->fill(xbj, deltaWeight);
        _h_DPt  ->fill(dijetPt, deltaWeight);
        _h_DMjj ->fill(Mjj, deltaWeight);
        _h_DEta ->fill(etastar, deltaWeight);
        _h_DZeta->fill(logZeta, deltaWeight);
        _h_DpT2OverpT1->fill(pT2OverpT1, deltaWeight);
        _h_DpTAsym->fill(pTAsym, deltaWeight);
      }
    }

    void finalize() {
      if (_jetInputMode == JetInputMode::MEPARTONS && _mePartonDebugTotals.events > 0) {
        MSG_INFO("JETINPUT=MEPARTONS diagnostic summary after DIS cuts: events=" +
                 std::to_string(_mePartonDebugTotals.events) +
                 ", source(tags/fallback/none)=" +
                 std::to_string(_mePartonDebugTotals.tagSource) + "/" +
                 std::to_string(_mePartonDebugTotals.fallbackSource) + "/" +
                 std::to_string(_mePartonDebugTotals.noSource) +
                 ", input>=2=" + std::to_string(_mePartonDebugTotals.inputGE2) +
                 ", jets>=2=" + std::to_string(_mePartonDebugTotals.jetsGE2) +
                 ", passPt=" + std::to_string(_mePartonDebugTotals.passPt) +
                 ", passRapidity=" + std::to_string(_mePartonDebugTotals.passRapidity) +
                 ", allTagged(role2/role3/any)=" +
                 std::to_string(_mePartonDebugTotals.taggedAnyRole2) + "/" +
                 std::to_string(_mePartonDebugTotals.taggedAnyRole3) + "/" +
                 std::to_string(_mePartonDebugTotals.taggedAny) +
                 ", graph(hard/root/terminal)=" +
                 std::to_string(_mePartonDebugTotals.graphHardPartons) + "/" +
                 std::to_string(_mePartonDebugTotals.graphRootSeeds) + "/" +
                 std::to_string(_mePartonDebugTotals.graphTerminalPartons));
      }
      if (_rivetWeightsMode) {
        MSG_INFO("RIVETWEIGHTS summary: events=" +
                 std::to_string(_rivetWeightsEvents) +
                 ", used=" + std::to_string(_rivetWeightsUsed) +
                 ", missing=" + std::to_string(_rivetWeightsMissing));
      }
      if (sumW() == 0.0) return;
      const double sf = crossSection()/picobarn/sumW();
      scale(_h_Q2,   sf);
      scale(_h_XBj,  sf);
      scale(_h_Pt,   sf);
      scale(_h_Mjj,  sf);
      scale(_h_Eta,  sf);
      scale(_h_Zeta, sf);
      scale(_h_pT1,  sf);
      scale(_h_pT2, sf);
      scale(_h_pT2OverpT1, sf);
      scale(_h_pTAsym, sf);
      scale(_h_Q2PreCut, sf);
      scale(_h_XBjPreCut, sf);
      scale(_h_YPreCut, sf);
      scale(_h_pT1PreCut, sf);
      scale(_h_pT2PreCut, sf);
      if (_rivetWeightsMode) {
        scale(_h_DQ2, sf);
        scale(_h_DXBj, sf);
        scale(_h_DPt, sf);
        scale(_h_DMjj, sf);
        scale(_h_DEta, sf);
        scale(_h_DZeta, sf);
        scale(_h_DpT1, sf);
        scale(_h_DpT2, sf);
        scale(_h_DpT2OverpT1, sf);
        scale(_h_DpTAsym, sf);
        scale(_h_DQ2PreCut, sf);
        scale(_h_DXBjPreCut, sf);
        scale(_h_DYPreCut, sf);
        scale(_h_DpT1PreCut, sf);
        scale(_h_DpT2PreCut, sf);
      }
      scale(_h_METaggedOutgoingCount, sf);
      scale(_h_METaggedRole2Count, sf);
      scale(_h_METaggedRole3Count, sf);
      scale(_h_METerminalFallbackCount, sf);
      scale(_h_MEInputPartonCount, sf);
      scale(_h_MEJetCountPreCut, sf);
      scale(_h_MEInputSource, sf);
      scale(_h_MESelectionStage, sf);
    }

    double defaultEventWeight(const Event& event) const {
      const auto weights = event.weights();
      const size_t defaultIndex = defaultWeightIndex();
      if (weights.size() > defaultIndex) return weights[defaultIndex];
      return weights.size() > 0 ? weights[0] : 1.0;
    }

    bool namedHepMCWeight(const Event& event, const string& name, double& value) const {
      const GenEvent* ge = event.genEvent();
      if (ge == nullptr) return false;

      const vector<string> names = HepMCUtils::weightNames(*ge);
      const std::valarray<double> weights = HepMCUtils::weights(*ge);
      const size_t n = std::min(names.size(), weights.size());
      for (size_t i = 0; i < n; ++i) {
        if (names[i] == name) {
          value = weights[i];
          return std::isfinite(value);
        }
      }
      return false;
    }

    bool namedEventWeight(const Event& event, const string& name, double& value) const {
      const vector<string>& names = handler().weightNames();
      const std::valarray<double> weights = event.weights();
      const size_t n = std::min(names.size(), weights.size());
      for (size_t i = 0; i < n; ++i) {
        if (names[i] == name) {
          value = weights[i];
          return std::isfinite(value);
        }
      }
      return namedHepMCWeight(event, name, value);
    }

    string availableWeightNamesSummary() const {
      const vector<string>& names = handler().weightNames();
      if (names.empty()) return "<none>";
      string out;
      const size_t n = std::min<size_t>(names.size(), 12);
      for (size_t i = 0; i < n; ++i) {
        if (!out.empty()) out += ", ";
        out += names[i].empty() ? "<empty>" : names[i];
      }
      if (names.size() > n) out += ", ...";
      return out;
    }

    AnalysisWeights eventAnalysisWeights(const Event& event) {
      AnalysisWeights out;
      const double defaultWeight = defaultEventWeight(event);
      out.sigma = defaultWeight;
      if (!_rivetWeightsMode) return out;

      ++_rivetWeightsEvents;
      double sigma = 0.0;
      double delta = 0.0;
      const bool haveSigma = namedEventWeight(event, "HERWIG_DIS_SIGMA", sigma);
      const bool haveDelta = namedEventWeight(event, "HERWIG_DIS_DELTA_LL", delta);
      if (haveSigma && haveDelta) {
        ++_rivetWeightsUsed;
        // Rivet multiplies each explicit fill weight by the selected event
        // weight at event collapse. The Herwig DIS optional weights are
        // correlated helicity multipliers for the Default event weight.
        out.sigma = sigma;
        out.deltaLL = delta;
        out.haveRivetWeights = true;
        return out;
      }

      ++_rivetWeightsMissing;
      if (!_warnedMissingRivetWeights) {
        MSG_WARNING("RIVETWEIGHTS=YES but HERWIG_DIS_SIGMA/HERWIG_DIS_DELTA_LL weights are missing. Available Rivet weight names: " + availableWeightNamesSummary() + ". Falling back to the default event weight for unpolarized histograms.");
        _warnedMissingRivetWeights = true;
      }
      return out;
    }

    vector<ClusteredJet> clusterFullFinalStateJets(const Event& event) const {
      return clusterJetInputs(collectBreitFinalStateInputs(event));
    }

    vector<ClusteredJet> clusterHardPartonJets(const Event& event) const {
      const vector<JetInputParticle> allInputs = collectBreitFinalStateInputs(event);
      vector<JetInputParticle> inputs;
      inputs.reserve(allInputs.size());
      for (const JetInputParticle& input : allInputs) {
        if (isHardParton(input.pid)) {
          inputs.push_back(input);
        }
      }
      return clusterJetInputs(inputs);
    }

    vector<ClusteredJet> clusterTopHardPartonJets(const Event& event, size_t maxPartons) const {
      vector<JetInputParticle> inputs;
      for (const JetInputParticle& input : collectBreitFinalStateInputs(event)) {
        if (isHardParton(input.pid)) {
          inputs.push_back(input);
        }
      }

      std::sort(inputs.begin(), inputs.end(),
                [](const JetInputParticle& a, const JetInputParticle& b) {
                  return a.breitMom.pT() > b.breitMom.pT();
                });

      if (inputs.size() > maxPartons) {
        inputs.resize(maxPartons);
      }
      return clusterJetInputs(inputs);
    }

    vector<ClusteredJet> clusterMEPartonJets(const Event& event) {
      _lastMEPartonDebug = MEPartonDebug();
      const vector<JetInputParticle> inputs = collectMEPartonInputs(event);
      _lastMEPartonDebug.inputPartons = inputs.size();
      return clusterJetInputs(inputs);
    }

    vector<JetInputParticle> jetInputsFromGenParticles(const DISKinematics& dis,
                                                       const vector<ConstGenParticlePtr>& particles) const {
      const LorentzTransform& breitBoost = dis.boostBreit();
      vector<JetInputParticle> inputs;
      inputs.reserve(particles.size());

      for (ConstGenParticlePtr parton : particles) {
        if (!parton) continue;
        const Particle pLab(parton);
        inputs.push_back({breitBoost.transform(pLab.momentum()), pLab.momentum(), pLab.pid()});
      }

      return inputs;
    }

    vector<JetInputParticle> collectBreitFinalStateInputs(const Event& event) const {
      const DISKinematics& dis = apply<DISKinematics>(event, "Kinematics");
      const Particles& labParticles = apply<FinalState>(event, "LabFS").particles();
      vector<JetInputParticle> inputs;
      inputs.reserve(labParticles.size());
      const auto scatteredLepton = dis.scatteredLepton().genParticle();
      const LorentzTransform& breitBoost = dis.boostBreit();

      for (const Particle& pLab : labParticles) {
        if (pLab.genParticle() == scatteredLepton) continue;
        inputs.push_back({breitBoost.transform(pLab.momentum()), pLab.momentum(), pLab.pid()});
      }

      return inputs;
    }

    vector<JetInputParticle> collectMEPartonInputs(const Event& event) {
      const DISKinematics& dis = apply<DISKinematics>(event, "Kinematics");
      const GenEvent* ge = event.genEvent();
      if (ge == nullptr) {
        warnMEPartonSelection("GenEvent pointer is null");
        return {};
      }

      recordAllMEPartonDebug(*ge);
      const vector<ConstGenParticlePtr> taggedPartons = collectTaggedPOWHEGMEPartons(*ge);
      recordTaggedMEPartonDebug(taggedPartons);
      if (taggedPartons.size() >= 2) {
        _lastMEPartonDebug.source = 1;
        return jetInputsFromGenParticles(dis, taggedPartons);
      }
      if (!taggedPartons.empty()) {
        warnMEPartonSelection("found only " + std::to_string(taggedPartons.size()) +
                              " tagged outgoing POWHEG ME parton; trying the signal-process hard-line fallback");
      }

      std::unordered_set<int> hardSeedIds;
      vector<ConstGenParticlePtr> hardSeeds;
      std::unordered_set<int> acceptedIds;
      vector<ConstGenParticlePtr> terminalPartons;
      const ConstGenVertexPtr signalVertex = resolveSignalProcessVertex(*ge, false);
      if (signalVertex) {
        // In Herwig POWHEG DIS the realized hard emission can sit on the
        // incoming colored line before the hard-scattering vertex, not only on
        // the outgoing Born leg. Keep the incoming hard-line ancestry as the
        // anchor, then select final-state partons attached to that hard system.
        for (ConstGenParticlePtr incoming : signalVertex->particles_in()) {
          if (!incoming || !isHardParton(incoming->pdg_id())) continue;
          collectHardLineAncestors(incoming, hardSeedIds, hardSeeds);
        }
        for (ConstGenParticlePtr outgoing : signalVertex->particles_out()) {
          if (!outgoing || !isHardParton(outgoing->pdg_id())) continue;
          const int uniqueId = HepMCUtils::uniqueId(outgoing);
          if (hardSeedIds.insert(uniqueId).second) {
            hardSeeds.push_back(outgoing);
          }
        }
        _lastMEPartonDebug.graphRootSeeds = hardSeeds.size();
        if (hardSeeds.empty()) {
          warnMEPartonSelection("signal-process fallback found no colored hard seeds");
          return {};
        }
        for (ConstGenParticlePtr seed : hardSeeds) {
          collectTerminalHardSystemPartons(seed, acceptedIds, terminalPartons, true);
        }
      } else {
        collectGraphTerminalHardPartons(*ge, acceptedIds, terminalPartons);
      }
      _lastMEPartonDebug.graphTerminalPartons = terminalPartons.size();

      if (terminalPartons.empty()) {
        warnMEPartonSelection("hard-line fallback found no terminal colored partons");
        return {};
      }
      _lastMEPartonDebug.terminalFallback = terminalPartons.size();
      if (terminalPartons.size() < 2) {
        warnMEPartonSelection("hard-line fallback found only " +
                              std::to_string(terminalPartons.size()) +
                              " terminal colored parton");
      }

      _lastMEPartonDebug.source = 2;
      return jetInputsFromGenParticles(dis, terminalPartons);
    }

    void recordAllMEPartonDebug(const GenEvent& ge) {
      for (ConstGenParticlePtr particle : ge.particles()) {
        if (!particle || !isHardParton(particle->pdg_id())) continue;
        const int role = powhegMEPartonRole(particle);
        if (!role) continue;
        ++_lastMEPartonDebug.taggedAny;
        if (role == 2) ++_lastMEPartonDebug.taggedAnyRole2;
        if (role == 3) ++_lastMEPartonDebug.taggedAnyRole3;
      }
    }

    void recordTaggedMEPartonDebug(const vector<ConstGenParticlePtr>& taggedPartons) {
      _lastMEPartonDebug.taggedOutgoing = taggedPartons.size();
      for (ConstGenParticlePtr parton : taggedPartons) {
        const int role = powhegMEPartonRole(parton);
        if (role == 2) ++_lastMEPartonDebug.taggedRole2;
        if (role == 3) ++_lastMEPartonDebug.taggedRole3;
      }
    }

    void fillMEPartonDiagnostics(const vector<ClusteredJet>& jets, double weight) {
      ++_mePartonDebugTotals.events;
      if (_lastMEPartonDebug.source == 1) {
        ++_mePartonDebugTotals.tagSource;
      } else if (_lastMEPartonDebug.source == 2) {
        ++_mePartonDebugTotals.fallbackSource;
      } else {
        ++_mePartonDebugTotals.noSource;
      }
      _mePartonDebugTotals.taggedAny += _lastMEPartonDebug.taggedAny;
      _mePartonDebugTotals.taggedAnyRole2 += _lastMEPartonDebug.taggedAnyRole2;
      _mePartonDebugTotals.taggedAnyRole3 += _lastMEPartonDebug.taggedAnyRole3;
      _mePartonDebugTotals.graphHardPartons += _lastMEPartonDebug.graphHardPartons;
      _mePartonDebugTotals.graphRootSeeds += _lastMEPartonDebug.graphRootSeeds;
      _mePartonDebugTotals.graphTerminalPartons += _lastMEPartonDebug.graphTerminalPartons;

      _h_METaggedOutgoingCount->fill(_lastMEPartonDebug.taggedOutgoing, weight);
      _h_METaggedRole2Count->fill(_lastMEPartonDebug.taggedRole2, weight);
      _h_METaggedRole3Count->fill(_lastMEPartonDebug.taggedRole3, weight);
      _h_METerminalFallbackCount->fill(_lastMEPartonDebug.terminalFallback, weight);
      _h_MEInputPartonCount->fill(_lastMEPartonDebug.inputPartons, weight);
      _h_MEJetCountPreCut->fill(jets.size(), weight);
      _h_MEInputSource->fill(_lastMEPartonDebug.source, weight);

      _h_MESelectionStage->fill(0.0, weight);
      if (_lastMEPartonDebug.inputPartons >= 2) {
        ++_mePartonDebugTotals.inputGE2;
        _h_MESelectionStage->fill(1.0, weight);
      }
      if (jets.size() < 2) return;

      ++_mePartonDebugTotals.jetsGE2;
      _h_MESelectionStage->fill(2.0, weight);
      const bool passPt = jets[0].breitMom.pT() >= 5*GeV && jets[1].breitMom.pT() >= 4*GeV;
      if (!passPt) return;

      ++_mePartonDebugTotals.passPt;
      _h_MESelectionStage->fill(3.0, weight);
      const bool passRapidity =
        inRange(jets[0].labMom.rapidity(), -3.5, 3.5) &&
        inRange(jets[1].labMom.rapidity(), -3.5, 3.5);
      if (passRapidity) {
        ++_mePartonDebugTotals.passRapidity;
        _h_MESelectionStage->fill(4.0, weight);
      }
    }

    vector<ClusteredJet> clusterJetInputs(const vector<JetInputParticle>& jetInputs) const {
      if (jetInputs.empty()) return {};

      vector<fastjet::PseudoJet> inputs;
      vector<FourMomentum> breitMomenta;
      vector<FourMomentum> labMomenta;
      inputs.reserve(jetInputs.size());
      breitMomenta.reserve(jetInputs.size());
      labMomenta.reserve(jetInputs.size());

      for (const JetInputParticle& jetInput : jetInputs) {
        fastjet::PseudoJet input(jetInput.breitMom.px(), jetInput.breitMom.py(),
                                 jetInput.breitMom.pz(), jetInput.breitMom.E());
        input.set_user_index(static_cast<int>(breitMomenta.size()));
        inputs.push_back(input);
        breitMomenta.push_back(jetInput.breitMom);
        labMomenta.push_back(jetInput.labMom);
      }

      fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, 1.0,
                                    fastjet::E_scheme, fastjet::Best);
      fastjet::ClusterSequence cs(inputs, jetDef);
      const vector<fastjet::PseudoJet> fjJets = fastjet::sorted_by_pt(cs.inclusive_jets());

      vector<ClusteredJet> jets;
      jets.reserve(fjJets.size());
      for (const fastjet::PseudoJet& fjJet : fjJets) {
        const vector<fastjet::PseudoJet> constituents = cs.constituents(fjJet);
        if (constituents.empty()) continue;

        const size_t first = static_cast<size_t>(constituents.front().user_index());
        FourMomentum breitMom = breitMomenta[first];
        FourMomentum labMom = labMomenta[first];
        for (size_t i = 1; i < constituents.size(); ++i) {
          const size_t idx = static_cast<size_t>(constituents[i].user_index());
          breitMom += breitMomenta[idx];
          labMom += labMomenta[idx];
        }
        jets.push_back({breitMom, labMom});
      }

      std::sort(jets.begin(), jets.end(),
                [](const ClusteredJet& a, const ClusteredJet& b) {
                  return a.breitMom.pT() > b.breitMom.pT();
                });
      return jets;
    }

  private:
    template <typename...>
    using VoidT = void;

    template <typename EventT, typename = void>
    struct HasSignalProcessVertexAccessor : std::false_type {};

    template <typename EventT>
    struct HasSignalProcessVertexAccessor<EventT, VoidT<decltype(std::declval<const EventT&>().signal_process_vertex())>>
      : std::true_type {};

    static bool isHardParton(int pid) {
      const int apid = std::abs(pid);
      return (apid >= 1 && apid <= 6) || pid == 21;
    }

    static bool isRemnant(int pid) {
      return pid == 82;
    }

    static int powhegMEPartonRoleFromStatus(int status) {
      if (status == 990071) return 1;
      if (status == 990072) return 2;
      if (status == 990073) return 3;
      return 0;
    }

    int powhegMEPartonRole(ConstGenParticlePtr particle) const {
      if (!particle) return 0;
      const auto attr = particle->attribute<RivetHepMC::IntAttribute>("herwig_powheg_me_parton");
      if (attr) return attr->value();
      return powhegMEPartonRoleFromStatus(particle->status());
    }

    bool isTaggedPOWHEGMEOutgoing(ConstGenParticlePtr particle) const {
      const int role = powhegMEPartonRole(particle);
      return role == 2 || role == 3;
    }

    bool hasTaggedPOWHEGMEContinuationChild(ConstGenParticlePtr particle) const {
      if (!particle) return false;
      const int role = powhegMEPartonRole(particle);
      const int pid = particle->pdg_id();
      for (ConstGenParticlePtr child : directChildren(particle)) {
        if (!child || !isHardParton(child->pdg_id())) continue;
        // Collapse only the copy chain of the same tagged ME leg.
        // The hard POWHEG branching itself can legitimately give a tagged
        // quark/gluon sibling pair. Treating any tagged child as a continuation
        // erroneously removes one of the two real-emission legs and leaves the
        // event with a single jet input.
        if (powhegMEPartonRole(child) == role && child->pdg_id() == pid) {
          return true;
        }
      }
      return false;
    }

    vector<ConstGenParticlePtr> collectTaggedPOWHEGMEPartons(const GenEvent& ge) const {
      std::unordered_set<int> acceptedIds;
      vector<ConstGenParticlePtr> taggedPartons;

      for (ConstGenParticlePtr particle : ge.particles()) {
        if (!particle || !isHardParton(particle->pdg_id())) continue;
        if (!isTaggedPOWHEGMEOutgoing(particle)) continue;
        // Copy chains inherit the generator-side tag. Keep the last tagged
        // partonic copy of each ME leg, but do not collapse across the actual
        // hard-emission sibling branch.
        if (hasTaggedPOWHEGMEContinuationChild(particle)) continue;

        const int uniqueId = HepMCUtils::uniqueId(particle);
        if (acceptedIds.insert(uniqueId).second) {
          taggedPartons.push_back(particle);
        }
      }

      return taggedPartons;
    }

    static ConstGenVertexPtr normalizeSignalProcessVertex(const ConstGenVertexPtr& vertex) {
      return vertex;
    }

    static ConstGenVertexPtr normalizeSignalProcessVertex(const std::shared_ptr<RivetHepMC::GenVertex>& vertex) {
      return vertex;
    }

    static ConstGenVertexPtr normalizeSignalProcessVertex(const RivetHepMC::GenVertex* vertex) {
      return vertex ? vertex->shared_from_this() : ConstGenVertexPtr();
    }

    template <typename EventT>
    static ConstGenVertexPtr directSignalProcessVertex(const EventT& ge, std::true_type) {
      return normalizeSignalProcessVertex(ge.signal_process_vertex());
    }

    template <typename EventT>
    static ConstGenVertexPtr directSignalProcessVertex(const EventT&, std::false_type) {
      return ConstGenVertexPtr();
    }

    ConstGenVertexPtr resolveSignalProcessVertex(const GenEvent& ge, bool warnIfMissing = true) {
      // Some HepMC compatibility layers expose a direct signal-process-vertex
      // handle. Rivet's HepMC3 event may not, so probe at compile time and
      // otherwise fall back to the ThePEG-written vertex-id attribute.
      const ConstGenVertexPtr directVertex =
        directSignalProcessVertex(ge, HasSignalProcessVertexAccessor<GenEvent>{});
      if (directVertex) return directVertex;

      const auto signalVertexAttr = ge.attribute<RivetHepMC::IntAttribute>("signal_process_vertex");
      if (!signalVertexAttr) {
        if (warnIfMissing) warnMEPartonSelection("missing signal_process_vertex metadata");
        return ConstGenVertexPtr();
      }

      const int signalVertexId = signalVertexAttr->value();
      for (ConstGenVertexPtr vertex : HepMCUtils::vertices(&ge)) {
        if (vertex && vertex->id() == signalVertexId) return vertex;
      }

      if (warnIfMissing) {
        warnMEPartonSelection("signal_process_vertex id " + std::to_string(signalVertexId) +
                              " did not resolve to a GenVertex");
      }
      return ConstGenVertexPtr();
    }

    void collectHardLineAncestors(ConstGenParticlePtr parton,
                                  std::unordered_set<int>& hardSeedIds,
                                  vector<ConstGenParticlePtr>& hardSeeds) const {
      if (!parton || !isHardParton(parton->pdg_id())) return;

      vector<ConstGenParticlePtr> stack{parton};
      while (!stack.empty()) {
        ConstGenParticlePtr current = stack.back();
        stack.pop_back();
        if (!current || !isHardParton(current->pdg_id())) continue;

        const int uniqueId = HepMCUtils::uniqueId(current);
        if (!hardSeedIds.insert(uniqueId).second) continue;
        hardSeeds.push_back(current);

        const ConstGenVertexPtr production = current->production_vertex();
        if (!production) continue;

        for (ConstGenParticlePtr parent : production->particles_in()) {
          if (!parent || !isHardParton(parent->pdg_id()) || isRemnant(parent->pdg_id())) continue;
          stack.push_back(parent);
        }
      }
    }

    bool hasHardPartonChild(ConstGenParticlePtr parton) const {
      if (!parton) return false;
      for (ConstGenParticlePtr child : directChildren(parton)) {
        if (child && isHardParton(child->pdg_id()) && !isRemnant(child->pdg_id())) return true;
      }
      return false;
    }

    vector<ConstGenParticlePtr> directChildren(ConstGenParticlePtr parton) const {
      vector<ConstGenParticlePtr> children;
      if (!parton) return children;

      const ConstGenVertexPtr decay = parton->end_vertex();
      if (!decay) return children;

      const int parentId = HepMCUtils::uniqueId(parton);
      for (ConstGenParticlePtr child : decay->particles_out()) {
        if (!child) continue;
        if (HepMCUtils::uniqueId(child) == parentId) continue;
        children.push_back(child);
      }
      return children;
    }

    void collectGraphTerminalHardPartons(const GenEvent& ge,
                                         std::unordered_set<int>& acceptedIds,
                                         vector<ConstGenParticlePtr>& terminals) {
      for (ConstGenParticlePtr particle : ge.particles()) {
        if (!particle || !isHardParton(particle->pdg_id()) || isRemnant(particle->pdg_id())) continue;
        ++_lastMEPartonDebug.graphHardPartons;
        if (hasHardPartonChild(particle)) continue;

        const int uniqueId = HepMCUtils::uniqueId(particle);
        if (acceptedIds.insert(uniqueId).second) {
          terminals.push_back(particle);
        }
      }
      _lastMEPartonDebug.graphTerminalPartons = terminals.size();
    }

    bool hasRemnantAncestor(ConstGenParticlePtr parton) const {
      if (!parton) return false;

      std::unordered_set<int> visited;
      vector<ConstGenParticlePtr> stack{parton};
      while (!stack.empty()) {
        ConstGenParticlePtr current = stack.back();
        stack.pop_back();
        if (!current) continue;

        const ConstGenVertexPtr production = current->production_vertex();
        if (!production) continue;

        for (ConstGenParticlePtr parent : production->particles_in()) {
          if (!parent) continue;
          const int uniqueId = HepMCUtils::uniqueId(parent);
          if (!visited.insert(uniqueId).second) continue;
          if (isRemnant(parent->pdg_id())) return true;
          stack.push_back(parent);
        }
      }

      return false;
    }

    void collectTerminalHardSystemPartons(ConstGenParticlePtr parton,
                                          std::unordered_set<int>& acceptedIds,
                                          vector<ConstGenParticlePtr>& terminals,
                                          bool rejectRemnantAncestors) const {
      if (!parton || !isHardParton(parton->pdg_id())) return;
      if (rejectRemnantAncestors && hasRemnantAncestor(parton)) return;

      vector<ConstGenParticlePtr> partonChildren;
      // Work with the HepMC copy graph directly: in the RIVETFO record the
      // physically relevant hard-system partons may survive only as terminal
      // partonic copies rather than Rivet FinalState particles.
      for (ConstGenParticlePtr child : directChildren(parton)) {
        if (!child || !isHardParton(child->pdg_id())) continue;
        if (rejectRemnantAncestors && hasRemnantAncestor(child)) continue;
        partonChildren.push_back(child);
      }

      if (partonChildren.empty()) {
        const int uniqueId = HepMCUtils::uniqueId(parton);
        if (acceptedIds.insert(uniqueId).second) {
          terminals.push_back(parton);
        }
        return;
      }

      for (ConstGenParticlePtr child : partonChildren) {
        collectTerminalHardSystemPartons(child, acceptedIds, terminals, rejectRemnantAncestors);
      }
    }

    void warnMEPartonSelection(const string& reason) {
      if (_warnedMEPartonSelection) return;
      MSG_WARNING("JETINPUT=MEPARTONS selection note: " + reason + ".");
      _warnedMEPartonSelection = true;
    }

    Histo1DPtr _h_Q2, _h_XBj, _h_Pt, _h_Mjj, _h_Eta, _h_Zeta, _h_pT1, _h_pT2,
      _h_pT2OverpT1, _h_pTAsym;
    Histo1DPtr _h_Q2PreCut, _h_XBjPreCut, _h_YPreCut, _h_pT1PreCut, _h_pT2PreCut;
    Histo1DPtr _h_DQ2, _h_DXBj, _h_DPt, _h_DMjj, _h_DEta, _h_DZeta,
      _h_DpT1, _h_DpT2, _h_DpT2OverpT1, _h_DpTAsym;
    Histo1DPtr _h_DQ2PreCut, _h_DXBjPreCut, _h_DYPreCut,
      _h_DpT1PreCut, _h_DpT2PreCut;
    Histo1DPtr _h_METaggedOutgoingCount, _h_METaggedRole2Count, _h_METaggedRole3Count,
      _h_METerminalFallbackCount, _h_MEInputPartonCount, _h_MEJetCountPreCut,
      _h_MEInputSource, _h_MESelectionStage;
    MEPartonDebug _lastMEPartonDebug;
    MEPartonDebugTotals _mePartonDebugTotals;
    JetInputMode _jetInputMode = JetInputMode::FULL;
    bool _rivetWeightsMode = false;
    bool _warnedMEPartonSelection = false;
    bool _warnedMissingRivetWeights = false;
    size_t _rivetWeightsEvents = 0;
    size_t _rivetWeightsUsed = 0;
    size_t _rivetWeightsMissing = 0;
  };

  RIVET_DECLARE_PLUGIN(MC_DIS_BREIT);

}
