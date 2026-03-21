// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "fastjet/ClusterSequence.hh"
#include <algorithm>
#include <cmath>

namespace Rivet {

  /// @brief DIS dijets with shower-sensitive observables in the Breit frame
  class MC_DIS_PS : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(MC_DIS_PS);

    struct ClusteredJet {
      FourMomentum breitMom;
      FourMomentum labMom;
    };

    struct JetInputParticle {
      FourMomentum breitMom;
      FourMomentum labMom;
      int pid;
    };

    enum class JetInputMode {
      FULL,
      HARDPARTONS,
      TOP2PARTONS,
      TOP3PARTONS
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
      } else {
        MSG_WARNING("Unknown JETINPUT option " + jetinputopt + ". Defaulting to FULL.");
        _jetInputMode = JetInputMode::FULL;
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

      const DISLepton dislep(Cuts::OPEN, LeptonReco::ALL, ObjOrdering::ENERGY,
                             0.0, 0.0, 0.0, dismode);
      const DISKinematics diskin(dislep);
      declare(diskin, "Kinematics");
      declare(FinalState(), "LabFS");
      declare(VisibleFinalState(), "VisibleFS");

      // Control observables inherited from MC_DIS_BREIT.
      book(_h_Q2, "Q2", 100, 49.0, 2500.0);
      book(_h_Pt, "Pt", 15, 5.0, 30.0);
      book(_h_XBj, "XBj", 20, 0.0, 1.0);
      book(_h_Mjj, "Mjj", logspace(15, 10.0, 100.0));
      book(_h_Eta, "Eta", 15, 0.0, 2.5);
      book(_h_Zeta, "Zeta", 12, -1.75, -0.25);
      book(_h_pT1, "pT1", 15, 5.0, 30.0);
      book(_h_pT2, "pT2", 15, 5.0, 30.0);
      book(_h_pT2OverpT1, "pT2OverpT1", 15, 0.0, 1.0);
      book(_h_pTAsym, "pTAsym", 15, 0.0, 1.0);
      book(_h_Q2PreCut, "Q2PreCut", 100, 49.0, 2500.0);
      book(_h_XBjPreCut, "XBjPreCut", 20, 0.0, 1.0);
      book(_h_YPreCut, "YPreCut", 40, 0.2, 0.6);
      book(_h_pT1PreCut, "pT1PreCut", 30, 0.0, 30.0);
      book(_h_pT2PreCut, "pT2PreCut", 30, 0.0, 30.0);

      // Shower-sensitive observables.
      book(_h_NJets, "NJets", 7, 1.5, 8.5);
      book(_h_pT3, "pT3", 15, 0.0, 30.0);
      book(_h_pT3OverpT1, "pT3OverpT1", 15, 0.0, 1.0);
      book(_h_SumPtExtra, "SumPtExtra", 20, 0.0, 40.0);
      book(_h_Phi3, "Phi3", 16, 0.0, M_PI);
      book(_h_PhiCurrentHemi, "PhiCurrentHemi", 20, 0.0, M_PI);
      book(_h_Broadening, "Broadening", 20, 0.0, 0.5);
      book(_h_Cos2PhiCurrentHemiNumQ2, "Cos2PhiCurrentHemiNumQ2", 10, 49.0, 2500.0);
      book(_h_Cos2PhiCurrentHemiDenQ2, "Cos2PhiCurrentHemiDenQ2", 10, 49.0, 2500.0);
    }

    void analyze(const Event& event) {
      const DISKinematics& dis = apply<DISKinematics>(event, "Kinematics");
      const auto weights = event.weights();
      const double weight = weights.size() > 0 ? weights[0] : 1.0;

      const double Q2 = dis.Q2();
      const double xbj = dis.x();
      const double y = dis.y();

      if (!inRange(Q2, 49*GeV2, 2500*GeV2)) vetoEvent;
      if (!inRange(y, 0.2, 0.6)) vetoEvent;

      vector<ClusteredJet> jets;
      if (_jetInputMode == JetInputMode::HARDPARTONS) {
        jets = clusterHardPartonJets(event);
      } else if (_jetInputMode == JetInputMode::TOP2PARTONS) {
        jets = clusterTopHardPartonJets(event, 2);
      } else if (_jetInputMode == JetInputMode::TOP3PARTONS) {
        jets = clusterTopHardPartonJets(event, 3);
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

      if (jets.size() < 2) vetoEvent;

      const ClusteredJet& jet1 = jets[0];
      const ClusteredJet& jet2 = jets[1];

      if (jet1.breitMom.pT() < 5*GeV || jet2.breitMom.pT() < 4*GeV) vetoEvent;
      if (!inRange(jet1.labMom.rapidity(), -3.5, 3.5)) vetoEvent;
      if (!inRange(jet2.labMom.rapidity(), -3.5, 3.5)) vetoEvent;

      const double dijetPt = 0.5 * (jet1.breitMom.pT() + jet2.breitMom.pT());
      const double pT2OverpT1 = jet2.breitMom.pT() / jet1.breitMom.pT();
      const double pTAsym = (jet1.breitMom.pT() - jet2.breitMom.pT()) /
                            (jet1.breitMom.pT() + jet2.breitMom.pT());
      const double Mjj = (jet1.breitMom + jet2.breitMom).mass();

      const int orientation = dis.orientation();
      const double eta1S = orientation * jet1.breitMom.eta();
      const double eta2S = orientation * jet2.breitMom.eta();
      const double etastar = 0.5 * std::abs(eta1S - eta2S);
      const double logZeta = std::log10(xbj * (1.0 + sqr(Mjj) / Q2));

      _h_Q2->fill(Q2, weight);
      _h_pT1->fill(jet1.breitMom.pT(), weight);
      _h_pT2->fill(jet2.breitMom.pT(), weight);
      _h_XBj->fill(xbj, weight);
      _h_Pt->fill(dijetPt, weight);
      _h_Mjj->fill(Mjj, weight);
      _h_Eta->fill(etastar, weight);
      _h_Zeta->fill(logZeta, weight);
      _h_pT2OverpT1->fill(pT2OverpT1, weight);
      _h_pTAsym->fill(pTAsym, weight);

      vector<ClusteredJet> acceptedJets;
      acceptedJets.reserve(jets.size());
      for (const ClusteredJet& jet : jets) {
        if (jet.breitMom.pT() <= 1.0*GeV) continue;
        if (!inRange(jet.labMom.rapidity(), -3.5, 3.5)) continue;
        acceptedJets.push_back(jet);
      }

      _h_NJets->fill(static_cast<double>(acceptedJets.size()), weight);

      double sumPtExtra = 0.0;
      for (size_t index = 2; index < acceptedJets.size(); ++index) {
        sumPtExtra += acceptedJets[index].breitMom.pT();
      }
      _h_SumPtExtra->fill(sumPtExtra, weight);

      const LorentzTransform& breitBoost = dis.boostBreit();
      const FourMomentum beamLeptonBreit = breitBoost.transform(dis.beamLepton().momentum());
      const double phiLepton = beamLeptonBreit.phi();
      const auto scatteredLepton = dis.scatteredLepton().genParticle();
      const Particles& visibleParticles = apply<VisibleFinalState>(event, "VisibleFS").particles();

      double broadeningPtSum = 0.0;
      double broadeningDenom = 0.0;
      for (const Particle& particle : visibleParticles) {
        if (particle.genParticle() == scatteredLepton) continue;
        const FourMomentum breitMom = breitBoost.transform(particle.momentum());
        if (orientation * breitMom.pz() >= 0.0) continue;

        const double dphi = foldToPi(breitMom.phi() - phiLepton);
        _h_PhiCurrentHemi->fill(dphi, weight);
        _h_Cos2PhiCurrentHemiNumQ2->fill(Q2, std::cos(2.0 * dphi) * weight);
        _h_Cos2PhiCurrentHemiDenQ2->fill(Q2, weight);

        broadeningPtSum += breitMom.pT();
        broadeningDenom += breitMom.vector3().mod();
      }

      if (broadeningDenom > 0.0) {
        _h_Broadening->fill(broadeningPtSum / (2.0 * broadeningDenom), weight);
      }

      if (acceptedJets.size() < 3) return;

      const ClusteredJet& jet3 = acceptedJets[2];
      _h_pT3->fill(jet3.breitMom.pT(), weight);
      _h_pT3OverpT1->fill(jet3.breitMom.pT() / jet1.breitMom.pT(), weight);
      const double phi3 = foldToPi(jet3.breitMom.phi() - phiLepton);
      _h_Phi3->fill(phi3, weight);
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double sf = crossSection() / picobarn / sumW();
      scale(_h_Q2, sf);
      scale(_h_XBj, sf);
      scale(_h_Pt, sf);
      scale(_h_Mjj, sf);
      scale(_h_Eta, sf);
      scale(_h_Zeta, sf);
      scale(_h_pT1, sf);
      scale(_h_pT2, sf);
      scale(_h_pT2OverpT1, sf);
      scale(_h_pTAsym, sf);
      scale(_h_Q2PreCut, sf);
      scale(_h_XBjPreCut, sf);
      scale(_h_YPreCut, sf);
      scale(_h_pT1PreCut, sf);
      scale(_h_pT2PreCut, sf);
      scale(_h_NJets, sf);
      scale(_h_pT3, sf);
      scale(_h_pT3OverpT1, sf);
      scale(_h_SumPtExtra, sf);
      scale(_h_Phi3, sf);
      scale(_h_PhiCurrentHemi, sf);
      scale(_h_Broadening, sf);
      scale(_h_Cos2PhiCurrentHemiNumQ2, sf);
      scale(_h_Cos2PhiCurrentHemiDenQ2, sf);
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
    static bool isHardParton(int pid) {
      return std::abs(pid) <= 6 || pid == 21;
    }

    static double foldToPi(double angle) {
      const double twoPi = 2.0 * M_PI;
      double folded = std::fmod(std::abs(angle), twoPi);
      if (folded > M_PI) folded = twoPi - folded;
      return folded;
    }

    Histo1DPtr _h_Q2, _h_XBj, _h_Pt, _h_Mjj, _h_Eta, _h_Zeta, _h_pT1, _h_pT2,
      _h_pT2OverpT1, _h_pTAsym;
    Histo1DPtr _h_Q2PreCut, _h_XBjPreCut, _h_YPreCut, _h_pT1PreCut, _h_pT2PreCut;
    Histo1DPtr _h_NJets, _h_pT3, _h_pT3OverpT1, _h_SumPtExtra, _h_Phi3,
      _h_PhiCurrentHemi, _h_Broadening, _h_Cos2PhiCurrentHemiNumQ2,
      _h_Cos2PhiCurrentHemiDenQ2;
    JetInputMode _jetInputMode = JetInputMode::FULL;
  };

  RIVET_DECLARE_PLUGIN(MC_DIS_PS);

}
