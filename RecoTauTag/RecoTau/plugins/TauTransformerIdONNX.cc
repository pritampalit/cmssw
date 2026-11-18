// TauTransformerIdONNX: computes the Tau Transformer (TauT) scores for HPS taus.
// The inputs are built in the same way as for the TauT training. Taus that are not HPS taus get -1.

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <memory>

#include "DataFormats/TauReco/interface/TauDiscriminatorContainer.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/AbstractServices/interface/ResourceInformation.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "RecoTauTag/RecoTau/interface/PFRecoTauClusterVariables.h"

namespace taut {

  struct Cache {
    std::unique_ptr<cms::Ort::ONNXRuntime> session;
  };

  // values used for missing quantities
  constexpr float kDefault = std::numeric_limits<float>::lowest();
  constexpr float kDefaultInt = float(std::numeric_limits<int>::lowest());
  constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();

  // true for a valid value
  inline bool isFiniteF(float x) {
    return std::isfinite(x) && x != std::numeric_limits<float>::lowest() && x != std::numeric_limits<float>::max();
  }
  inline bool isFiniteI(float x) { return x != kDefaultInt; }
  // replace invalid values
  inline float nanToNum(float x) {
    if (std::isnan(x))
      return 0.f;
    if (std::isinf(x))
      return x > 0 ? std::numeric_limits<float>::max() : std::numeric_limits<float>::lowest();
    return x;
  }
  // difference in phi
  inline float dPhi(float phi, float tauPhi) {
    float d = phi - tauPhi;
    if (!(d <= float(M_PI)))
      d -= float(2 * M_PI);
    if (!(d >= float(-M_PI)))
      d += float(2 * M_PI);
    return d;
  }
  // particle type code used in the training
  inline int pfParticleType(int pdgId) {
    switch (std::abs(pdgId)) {
      case 211:
        return 1;
      case 11:
        return 2;
      case 13:
        return 3;
      case 22:
        return 4;
      case 130:
        return 5;
      case 1:
        return 6;
      case 2:
        return 7;
      default:
        return 0;
    }
  }

  // fills the inputs of one particle
  class Row {
  public:
    Row(const std::map<std::string, size_t>& idx, float* out) : idx_(idx), out_(out) {}
    void set(const std::string& name, float value) {
      auto it = idx_.find(name);
      if (it != idx_.end())
        out_[it->second] = nanToNum(value);
    }

  private:
    const std::map<std::string, size_t>& idx_;
    float* out_;
  };

  struct MatchSetup {  // cone sizes
    double tau_boostedTau_dR = 0.2, tau_jet_dR = 0.4, tau_fatJet_dR = 0.8, jet_fatJet_dR = 0.8;
    double jet_maxAbsEta = 3.4, fatJet_maxAbsEta = 3.8;
    double tau_cone = 0.5, boostedTau_cone = 0.5, jet_cone = 0.8, fatJet_cone = 0.8;
  };

  // input names, in the order expected by the model
  const std::vector<std::string> pfCandFeatures = {"r", "theta", "rel_pt", "particle_type", "charge", "pvAssociationQuality", "fromPV", "tauLeadChargedHadrCand", "puppiWeight", "lostInnerHits", "nPixelHits", "hasTrackDetails", "nHits", "nPixelLayers", "nStripLayers", "vertex_dx", "vertex_dy", "vertex_z_valid", "vertex_dz", "vertex_dx_tauFL", "vertex_dy_tauFL", "vertex_dz_tauFL", "dxy", "dxy_sig_valid", "dxy_sig", "dz_valid", "dz", "dz_sig_valid", "dz_sig", "track_ndof_valid", "track_ndof", "chi2_ndof", "hcalFraction", "rawCaloFraction", "rawHcalFraction"};
  const std::vector<std::string> eleFeatures = {"r", "theta", "rel_pt", "particle_type", "cc_valid", "cc_ele_rel_energy", "cc_gamma_rel_energy", "cc_n_gamma", "rel_trackMomentumAtVtx", "rel_trackMomentumAtCalo", "rel_trackMomentumOut", "rel_trackMomentumAtEleClus", "rel_trackMomentumAtVtxWithConstraint", "rel_ecalEnergy", "ecalEnergy_sig", "eSuperClusterOverP", "eSeedClusterOverP", "eSeedClusterOverPout", "eEleClusterOverPout", "deltaEtaSuperClusterTrackAtVtx", "deltaEtaSeedClusterTrackAtCalo", "deltaEtaEleClusterTrackAtCalo", "deltaPhiEleClusterTrackAtCalo", "deltaPhiSuperClusterTrackAtVtx", "deltaPhiSeedClusterTrackAtCalo", "mva_valid", "mvaInput_earlyBrem", "mvaInput_lateBrem", "mvaInput_sigmaEtaEta", "mvaInput_hadEnergy", "mvaInput_deltaEta", "gsfTrack_normalizedChi2", "gsfTrack_numberOfValidHits", "rel_gsfTrack_pt", "gsfTrack_pt_sig", "has_closestCtfTrack", "closestCtfTrack_normalizedChi2", "closestCtfTrack_numberOfValidHits", "sigmaEtaEta", "sigmaIetaIeta", "sigmaIphiIphi", "sigmaIetaIphi", "e1x5", "e2x5Max", "e5x5", "r9", "hcalDepth1OverEcal", "hcalDepth2OverEcal", "hcalDepth1OverEcalBc", "hcalDepth2OverEcalBc", "eLeft", "eRight", "eBottom", "eTop", "full5x5_sigmaEtaEta", "full5x5_sigmaIetaIeta", "full5x5_sigmaIphiIphi", "full5x5_sigmaIetaIphi", "full5x5_e1x5", "full5x5_e2x5Max", "full5x5_e5x5", "full5x5_r9", "full5x5_hcalDepth1OverEcal", "full5x5_hcalDepth2OverEcal", "full5x5_hcalDepth1OverEcalBc", "full5x5_hcalDepth2OverEcalBc", "full5x5_eLeft", "full5x5_eRight", "full5x5_eBottom", "full5x5_eTop", "full5x5_e2x5Left", "full5x5_e2x5Right", "full5x5_e2x5Bottom", "full5x5_e2x5Top"};
  const std::vector<std::string> muonFeatures = {"r", "theta", "rel_pt", "particle_type", "dxy", "dxy_sig", "normalizedChi2_valid", "normalizedChi2", "numberOfValidHits", "segmentCompatibility", "caloCompatibility", "rel_pfEcalEnergy", "n_matches_DT_1", "n_matches_DT_2", "n_matches_DT_3", "n_matches_DT_4", "n_matches_CSC_1", "n_matches_CSC_2", "n_matches_CSC_3", "n_matches_CSC_4", "n_matches_RPC_1", "n_matches_RPC_2", "n_matches_RPC_3", "n_matches_RPC_4", "n_hits_DT_1", "n_hits_DT_2", "n_hits_DT_3", "n_hits_DT_4", "n_hits_CSC_1", "n_hits_CSC_2", "n_hits_CSC_3", "n_hits_CSC_4", "n_hits_RPC_1", "n_hits_RPC_2", "n_hits_RPC_3", "n_hits_RPC_4"};
  const std::vector<std::string> globalFeatures = {"particle_type", "rho", "tau_pt", "tau_eta", "tau_mass", "tau_E_over_pt", "tau_charge", "tau_n_charged_prongs", "tau_n_neutral_prongs", "tau_chargedIsoPtSum", "tau_chargedIsoPtSumdR03_over_dR05", "tau_footprintCorrection", "tau_neutralIsoPtSum", "tau_neutralIsoPtSumWeight_over_neutralIsoPtSum", "tau_neutralIsoPtSumWeightdR03_over_neutralIsoPtSum", "tau_neutralIsoPtSumdR03_over_dR05", "tau_photonPtSumOutsideSignalCone", "tau_puCorrPtSum", "tau_hasSecondaryVertex", "tau_dxy_valid", "tau_dxy", "tau_dxy_sig", "tau_ip3d", "tau_ip3d_sig", "tau_dz_sig_valid", "tau_dz", "tau_dz_sig", "tau_flightLength_x", "tau_flightLength_y", "tau_flightLength_z", "tau_flightLength_sig", "tau_pt_weighted_deta_strip", "tau_pt_weighted_dphi_strip", "tau_pt_weighted_dr_signal", "tau_pt_weighted_dr_iso", "tau_e_ratio", "tau_gj_angle_diff_valid", "tau_gj_angle_diff", "tau_n_photons", "tau_emFraction", "tau_inside_ecal_crack", "tau_leadingTrackNormChi2", "tau_leadChargedCand_etaAtEcalEntrance_minus_tau_eta"};

  // muon chamber matches and hits per station
  struct MuonStationCounts {
    unsigned matches[3][4] = {}, hits[3][4] = {};
  };
  inline MuonStationCounts countMuonStations(const pat::Muon& m) {
    MuonStationCounts c;
    for (const auto& seg : m.matches()) {
      if (seg.segmentMatches.empty() && seg.rpcMatches.empty() && seg.gemMatches.empty() && seg.me0Matches.empty())
        continue;
      const int det = seg.detector(), st = seg.station();
      if (det >= 1 && det <= 3 && st >= 1 && st <= 4)
        ++c.matches[det - 1][st - 1];
    }
    if (m.outerTrack().isNonnull()) {
      const auto& hp = m.outerTrack()->hitPattern();
      for (int i = 0; i < hp.numberOfAllHits(reco::HitPattern::TRACK_HITS); ++i) {
        const auto id = hp.getHitPattern(reco::HitPattern::TRACK_HITS, i);
        if (id == 0)
          break;
        const auto type = hp.getHitType(id);
        if (!hp.muonHitFilter(id) || !(type == TrackingRecHit::valid || type == TrackingRecHit::bad))
          continue;
        const int st = hp.getMuonStation(id);
        const int det = hp.muonDTHitFilter(id) ? 0 : hp.muonCSCHitFilter(id) ? 1 : hp.muonRPCHitFilter(id) ? 2 : -1;
        if (det >= 0 && st >= 1 && st <= 4)
          ++c.hits[det][st - 1];
      }
    }
    return c;
  }

  // outputs
  constexpr const char* kOutputs[7] = {"probTau", "probEle", "probMu", "probJet", "VSjet", "VSe", "VSmu"};

}  // namespace taut

class TauTransformerIdONNX : public edm::stream::EDProducer<edm::GlobalCache<taut::Cache>> {
public:
  explicit TauTransformerIdONNX(const edm::ParameterSet& cfg, const taut::Cache*)
      : tausToken_(consumes<pat::TauCollection>(cfg.getParameter<edm::InputTag>("taus"))),
        boostedTausToken_(consumes<pat::TauCollection>(cfg.getParameter<edm::InputTag>("boostedTaus"))),
        jetsToken_(consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("jets"))),
        fatJetsToken_(consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("fatJets"))),
        candsToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("pfCandidates"))),
        electronsToken_(consumes<pat::ElectronCollection>(cfg.getParameter<edm::InputTag>("electrons"))),
        muonsToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
        verticesToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
        svToken_(consumes<std::vector<reco::VertexCompositePtrCandidate>>(cfg.getParameter<edm::InputTag>("secondVertices"))),
        rhoToken_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
        useBoostedTaus_(cfg.getParameter<bool>("useBoostedTaus")),
        dumpFile_(cfg.getParameter<std::string>("debugDumpFile")) {
    for (const auto& [coll, key] : {std::pair<int, const char*>{0, "pfCandFeatures"},
                                    {1, "eleFeatures"},
                                    {2, "muonFeatures"},
                                    {3, "globalFeatures"}}) {
      const auto names = cfg.getParameter<std::vector<std::string>>(key);
      for (size_t i = 0; i < names.size(); ++i)
        featIdx_[coll][names[i]] = i;
      nFeat_[coll] = names.size();
    }
    for (const char* name : taut::kOutputs)
      produces<reco::TauDiscriminatorContainer>(name);
    if (!dumpFile_.empty())
      dump_ = std::make_unique<std::ofstream>(dumpFile_, std::ios::binary);
  }

  static std::unique_ptr<taut::Cache> initializeGlobalCache(const edm::ParameterSet& cfg) {
    cms::Ort::Backend backend = cms::Ort::Backend::cpu;
    edm::Service<edm::ResourceInformation> ri;
    if (ri.isAvailable() && ri->hasGpuNvidia())
      backend = cms::Ort::Backend::cuda;
    edm::LogPrint("TauTransformerIdONNX") << "running TauT on " << (backend == cms::Ort::Backend::cuda ? "CUDA" : "CPU");
    auto options = cms::Ort::ONNXRuntime::defaultSessionOptions(backend);
    auto cache = std::make_unique<taut::Cache>();
    cache->session = std::make_unique<cms::Ort::ONNXRuntime>(
        edm::FileInPath(cfg.getParameter<std::string>("model")).fullPath(), &options);
    return cache;
  }
  static void globalEndJob(const taut::Cache*) {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("taus", edm::InputTag("slimmedTaus"));
    desc.add<edm::InputTag>("boostedTaus", edm::InputTag("slimmedTausBoosted"));
    desc.add<bool>("useBoostedTaus", true);
    desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
    desc.add<edm::InputTag>("fatJets", edm::InputTag("slimmedJetsAK8"));
    desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
    desc.add<edm::InputTag>("electrons", edm::InputTag("slimmedElectrons"));
    desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
    desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
    desc.add<edm::InputTag>("secondVertices", edm::InputTag("slimmedSecondaryVertices"));
    desc.add<edm::InputTag>("rho", edm::InputTag("fixedGridRhoAll"));
    desc.add<std::string>("model", "RecoTauTag/RecoTau/data/taut_tim-ep59.onnx");
    desc.add<std::vector<std::string>>("pfCandFeatures", taut::pfCandFeatures);
    desc.add<std::vector<std::string>>("eleFeatures", taut::eleFeatures);
    desc.add<std::vector<std::string>>("muonFeatures", taut::muonFeatures);
    desc.add<std::vector<std::string>>("globalFeatures", taut::globalFeatures);
    desc.add<std::string>("debugDumpFile", "");
    descriptions.add("tauTransformerIdONNX", desc);
  }

  void produce(edm::Event& event, const edm::EventSetup&) override;

private:
  struct Tokens {  // particles selected for one tau
    std::vector<size_t> cands, electrons, muons;
  };
  Tokens selectObjects(const pat::Tau& tau,
                       const pat::TauCollection& boostedTaus,
                       const pat::JetCollection& jets,
                       const pat::JetCollection& fatJets,
                       const pat::PackedCandidateCollection& cands,
                       const pat::ElectronCollection& electrons,
                       const pat::MuonCollection& muons,
                       const std::vector<reco::VertexCompositePtrCandidate>& svs) const;
  void fillPfCand(float* out, const pat::PackedCandidate& c, const pat::Tau& tau, const reco::Vertex& pv) const;
  void fillElectron(float* out, const pat::Electron& e, const pat::Tau& tau) const;
  void fillMuon(float* out, const pat::Muon& m, const pat::Tau& tau) const;
  void fillGlobal(float* out, const pat::Tau& tau, float rho) const;

  edm::EDGetTokenT<pat::TauCollection> tausToken_, boostedTausToken_;
  edm::EDGetTokenT<pat::JetCollection> jetsToken_, fatJetsToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> candsToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate>> svToken_;
  edm::EDGetTokenT<double> rhoToken_;
  const bool useBoostedTaus_;
  const std::string dumpFile_;
  std::unique_ptr<std::ofstream> dump_;
  std::map<std::string, size_t> featIdx_[4];
  size_t nFeat_[4];
  taut::MatchSetup setup_;
};

TauTransformerIdONNX::Tokens TauTransformerIdONNX::selectObjects(const pat::Tau& tau,
                                                                 const pat::TauCollection& boostedTaus,
                                                                 const pat::JetCollection& jets,
                                                                 const pat::JetCollection& fatJets,
                                                                 const pat::PackedCandidateCollection& cands,
                                                                 const pat::ElectronCollection& electrons,
                                                                 const pat::MuonCollection& muons,
                                                                 const std::vector<reco::VertexCompositePtrCandidate>& svs) const {
  // jets matched to the tau
  const pat::Tau* boostedTau = nullptr;
  double best = std::numeric_limits<double>::infinity();
  if (useBoostedTaus_)
    for (const auto& bt : boostedTaus) {
      const double dr = reco::deltaR(bt.polarP4(), tau.polarP4());
      if (dr < setup_.tau_boostedTau_dR && dr < best) {
        best = dr;
        boostedTau = &bt;
      }
    }
  const pat::Jet* jet = nullptr;
  best = std::numeric_limits<double>::infinity();
  for (const auto& j : jets) {
    if (!(std::abs(j.polarP4().eta()) < setup_.jet_maxAbsEta))
      continue;
    const double dr = reco::deltaR(j.polarP4(), tau.polarP4());
    if (dr < setup_.tau_jet_dR && dr < best) {
      best = dr;
      jet = &j;
    }
  }
  const pat::Jet* fatJet = nullptr;
  double bestTau = std::numeric_limits<double>::infinity(), bestJet = std::numeric_limits<double>::infinity();
  for (const auto& fj : fatJets) {
    if (!(std::abs(fj.polarP4().eta()) < setup_.fatJet_maxAbsEta))
      continue;
    double drTau = reco::deltaR(fj.polarP4(), tau.polarP4());
    double drJet = jet ? reco::deltaR(fj.polarP4(), jet->polarP4()) : std::numeric_limits<double>::infinity();
    if (!(drTau < setup_.tau_fatJet_dR))
      drTau = std::numeric_limits<double>::infinity();
    if (!(drJet < setup_.jet_fatJet_dR))
      drJet = std::numeric_limits<double>::infinity();
    if (std::isinf(drTau) && std::isinf(drJet))
      continue;
    if (drTau < bestTau || (drTau == bestTau && drJet < bestJet)) {
      bestTau = drTau;
      bestJet = drJet;
      fatJet = &fj;
    }
  }
  const auto hasMatch = [&](const auto& p4) {
    if (reco::deltaR(p4, tau.polarP4()) < setup_.tau_cone)
      return true;
    if (boostedTau && reco::deltaR(p4, boostedTau->polarP4()) < setup_.boostedTau_cone)
      return true;
    if (jet && reco::deltaR(p4, jet->polarP4()) < setup_.jet_cone)
      return true;
    if (fatJet && reco::deltaR(p4, fatJet->polarP4()) < setup_.fatJet_cone)
      return true;
    return false;
  };
  const auto svHasCand = [](const reco::VertexCompositePtrCandidate& sv, const pat::PackedCandidate* c) {
    for (size_t k = 0; k < sv.numberOfSourceCandidatePtrs(); ++k)
      if (&(*sv.sourceCandidatePtr(k)) == c)
        return true;
    return false;
  };
  // particles around the tau
  std::map<size_t, bool> selCands, selSVs;
  size_t nCands = 0, nSVs = 0;
  while (true) {
    for (size_t i = 0; i < cands.size(); ++i) {
      if (selCands.count(i))
        continue;
      bool inSV = false;
      for (const auto& [k, _] : selSVs)
        if (svHasCand(svs.at(k), &cands[i])) {
          inSV = true;
          break;
        }
      if (inSV || hasMatch(cands[i].polarP4()))
        selCands[i] = true;
    }
    for (size_t k = 0; k < svs.size(); ++k) {
      if (selSVs.count(k))
        continue;
      bool hasSel = false;
      for (const auto& [i, _] : selCands)
        if (svHasCand(svs[k], &cands[i])) {
          hasSel = true;
          break;
        }
      if (hasSel || hasMatch(svs[k].polarP4()))
        selSVs[k] = true;
    }
    if (nCands == selCands.size() && nSVs == selSVs.size())
      break;
    nCands = selCands.size();
    nSVs = selSVs.size();
  }
  Tokens t;
  for (const auto& [i, _] : selCands)
    t.cands.push_back(i);
  for (size_t i = 0; i < electrons.size(); ++i)
    if (hasMatch(electrons[i].polarP4()))
      t.electrons.push_back(i);
  for (size_t i = 0; i < muons.size(); ++i)
    if (hasMatch(muons[i].polarP4()))
      t.muons.push_back(i);
  return t;
}

void TauTransformerIdONNX::fillPfCand(float* out, const pat::PackedCandidate& c, const pat::Tau& tau, const reco::Vertex& pv) const {
  using namespace taut;
  Row r(featIdx_[0], out);
  const float tauPt = tau.polarP4().pt(), tauEta = tau.polarP4().eta(), tauPhi = tau.polarP4().phi();
  const float dphi = dPhi(float(c.polarP4().phi()), tauPhi), deta = float(c.polarP4().eta()) - tauEta;
  r.set("r", std::sqrt(deta * deta + dphi * dphi));
  r.set("theta", std::atan2(dphi, deta));
  r.set("rel_pt", float(c.polarP4().pt()) / tauPt);
  r.set("particle_type", float(pfParticleType(c.pdgId()) - 1));
  r.set("charge", float(c.charge()));
  r.set("pvAssociationQuality", float(int(c.pvAssociationQuality())));
  r.set("fromPV", float(int(c.fromPV())));
  const auto* lead = dynamic_cast<const pat::PackedCandidate*>(tau.leadChargedHadrCand().get());
  r.set("tauLeadChargedHadrCand", float(lead == &c));
  r.set("puppiWeight", c.puppiWeight());
  r.set("lostInnerHits", float(int(c.lostInnerHits())));
  r.set("nPixelHits", float(c.numberOfPixelHits()));
  r.set("nHits", float(c.numberOfHits()));
  r.set("nPixelLayers", float(c.pixelLayersWithMeasurement()));
  r.set("nStripLayers", float(c.stripLayersWithMeasurement()));
  const float vx = c.vertex().x(), vy = c.vertex().y(), vz = c.vertex().z();
  const float pvx = pv.position().x(), pvy = pv.position().y(), pvz = pv.position().z();
  const float flx = tau.flightLength().x(), fly = tau.flightLength().y(), flz = tau.flightLength().z();
  const bool vzValid = isFiniteF(vz);
  r.set("vertex_dx", vx - pvx);
  r.set("vertex_dy", vy - pvy);
  r.set("vertex_z_valid", float(vzValid));
  r.set("vertex_dz", vzValid ? vz - pvz : kNaN);
  r.set("vertex_dx_tauFL", vx - pvx - flx);
  r.set("vertex_dy_tauFL", vy - pvy - fly);
  r.set("vertex_dz_tauFL", vzValid ? vz - pvz - flz : kNaN);
  const bool trk = c.hasTrackDetails();
  r.set("hasTrackDetails", float(trk));
  const float dxy = c.dxy(), dxyErr = trk ? c.dxyError() : kDefault, dz = c.dz(), dzErr = trk ? c.dzError() : kDefault;
  const float ndof = trk ? float(c.bestTrack()->ndof()) : kDefault, chi2 = trk ? float(c.bestTrack()->chi2()) : kDefault;
  const bool dxyValid = trk && isFiniteF(dxy) && isFiniteF(dxyErr), dzValid = trk && isFiniteF(dz) && isFiniteF(dzErr);
  const bool ndofValid = trk && ndof > 0;
  r.set("dxy", dxy);
  r.set("dxy_sig_valid", float(dxyValid));
  r.set("dxy_sig", dxyValid ? std::abs(dxy) / dxyErr : kNaN);
  r.set("dz_sig_valid", float(dzValid));
  r.set("dz_sig", dzValid ? std::abs(dz) / dzErr : kNaN);
  r.set("dz_valid", float(dzValid));
  r.set("dz", isFiniteF(dz) ? dz : kNaN);
  r.set("track_ndof_valid", float(ndofValid));
  r.set("track_ndof", ndofValid ? ndof : kNaN);
  r.set("chi2_ndof", ndofValid ? chi2 / ndof : kNaN);
  r.set("hcalFraction", c.hcalFraction());
  r.set("rawCaloFraction", c.rawCaloFraction());
  r.set("rawHcalFraction", c.rawHcalFraction());
}

void TauTransformerIdONNX::fillElectron(float* out, const pat::Electron& e, const pat::Tau& tau) const {
  using namespace taut;
  Row r(featIdx_[1], out);
  const float pt = e.polarP4().pt();
  const float dphi = dPhi(float(e.polarP4().phi()), float(tau.polarP4().phi())), deta = float(e.polarP4().eta()) - float(tau.polarP4().eta());
  r.set("r", std::sqrt(deta * deta + dphi * dphi));
  r.set("theta", std::atan2(dphi, deta));
  r.set("rel_pt", pt / float(tau.polarP4().pt()));
  r.set("particle_type", 7.f);
  // electron cluster quantities
  float ccEle = 0, ccGamma = 0, ccN = 0;
  const auto& sc = e.superCluster();
  if (sc.isNonnull() && sc.isAvailable() && sc->clusters().isNonnull() && sc->clusters().isAvailable()) {
    for (auto it = sc->clustersBegin(); it != sc->clustersEnd(); ++it) {
      const float en = (*it)->energy();
      if (it == sc->clustersBegin())
        ccEle += en;
      else {
        ccGamma += en;
        ++ccN;
      }
    }
  } else {
    ccEle = ccGamma = kDefault;
    ccN = kDefaultInt;
  }
  const bool ccValid = ccEle >= 0;
  r.set("cc_valid", float(ccValid));
  r.set("cc_ele_rel_energy", ccValid ? ccEle / pt : kNaN);
  r.set("cc_gamma_rel_energy", ccValid ? ccGamma / ccEle : kNaN);
  r.set("cc_n_gamma", ccValid ? ccN : kNaN);
  r.set("rel_trackMomentumAtVtx", float(e.trackMomentumAtVtx().R()) / pt);
  r.set("rel_trackMomentumAtCalo", float(e.trackMomentumAtCalo().R()) / pt);
  r.set("rel_trackMomentumOut", float(e.trackMomentumOut().R()) / pt);
  r.set("rel_trackMomentumAtEleClus", float(e.trackMomentumAtEleClus().R()) / pt);
  r.set("rel_trackMomentumAtVtxWithConstraint", float(e.trackMomentumAtVtxWithConstraint().R()) / pt);
  r.set("rel_ecalEnergy", float(e.ecalEnergy()) / pt);
  r.set("ecalEnergy_sig", float(e.ecalEnergy()) / float(e.ecalEnergyError()));
  r.set("eSuperClusterOverP", e.eSuperClusterOverP());
  r.set("eSeedClusterOverP", e.eSeedClusterOverP());
  r.set("eSeedClusterOverPout", e.eSeedClusterOverPout());
  r.set("eEleClusterOverPout", e.eEleClusterOverPout());
  r.set("deltaEtaSuperClusterTrackAtVtx", e.deltaEtaSuperClusterTrackAtVtx());
  r.set("deltaEtaSeedClusterTrackAtCalo", e.deltaEtaSeedClusterTrackAtCalo());
  r.set("deltaEtaEleClusterTrackAtCalo", e.deltaEtaEleClusterTrackAtCalo());
  r.set("deltaPhiEleClusterTrackAtCalo", e.deltaPhiEleClusterTrackAtCalo());
  r.set("deltaPhiSuperClusterTrackAtVtx", e.deltaPhiSuperClusterTrackAtVtx());
  r.set("deltaPhiSeedClusterTrackAtCalo", e.deltaPhiSeedClusterTrackAtCalo());
  const bool isHGCAL = e.hasUserFloat("hgcElectronID:sigmaUU");
  r.set("mvaInput_earlyBrem", !isHGCAL ? float(e.mvaInput().earlyBrem) : kDefaultInt);
  r.set("mvaInput_lateBrem", !isHGCAL ? float(e.mvaInput().lateBrem) : kDefaultInt);
  r.set("mvaInput_sigmaEtaEta", !isHGCAL ? e.mvaInput().sigmaEtaEta : kDefault);
  r.set("mvaInput_hadEnergy", !isHGCAL ? e.mvaInput().hadEnergy : kDefault);
  r.set("mvaInput_deltaEta", !isHGCAL ? e.mvaInput().deltaEta : kDefault);
  const auto& gsf = e.gsfTrack();
  const float gsfPt = gsf.isNonnull() ? float(gsf->pt()) : kDefault, gsfPtErr = gsf.isNonnull() ? float(gsf->ptError()) : kDefault;
  r.set("gsfTrack_normalizedChi2", gsf.isNonnull() ? float(gsf->normalizedChi2()) : kDefault);
  r.set("gsfTrack_numberOfValidHits", gsf.isNonnull() ? float(gsf->numberOfValidHits()) : kDefaultInt);
  r.set("rel_gsfTrack_pt", gsfPt / pt);
  r.set("gsfTrack_pt_sig", gsfPt / gsfPtErr);
  const auto& ctf = e.closestCtfTrackRef();
  const float ctfChi2 = ctf.isNonnull() ? float(ctf->normalizedChi2()) : kDefault;
  const bool hasCtf = ctfChi2 >= 0;
  r.set("has_closestCtfTrack", float(hasCtf));
  r.set("closestCtfTrack_normalizedChi2", hasCtf ? ctfChi2 : kNaN);
  r.set("closestCtfTrack_numberOfValidHits", hasCtf ? float(ctf->numberOfValidHits()) : kNaN);
  // electron shower shapes
  const bool shape = !isHGCAL && pt > 5;
  const float e5x5 = shape ? e.e5x5() : kDefault;
  const bool mvaValid = isFiniteF(e5x5);
  r.set("mva_valid", float(mvaValid));
  const auto S = [&](const char* name, float v) { r.set(name, mvaValid ? (shape ? v : kDefault) : kNaN); };
  S("sigmaEtaEta", e.sigmaEtaEta());
  S("sigmaIetaIeta", e.sigmaIetaIeta());
  S("sigmaIphiIphi", e.sigmaIphiIphi());
  S("sigmaIetaIphi", e.sigmaIetaIphi());
  S("e1x5", e.e1x5());
  S("e2x5Max", e.e2x5Max());
  S("e5x5", e.e5x5());
  S("r9", e.r9());
  S("hcalDepth1OverEcal", e.hcalOverEcal(1));
  S("hcalDepth2OverEcal", e.hcalOverEcal(2));
  S("hcalDepth1OverEcalBc", e.hcalOverEcalBc(1));
  S("hcalDepth2OverEcalBc", e.hcalOverEcalBc(2));
  S("eLeft", e.eLeft());
  S("eRight", e.eRight());
  S("eBottom", e.eBottom());
  S("eTop", e.eTop());
  S("full5x5_sigmaEtaEta", e.full5x5_sigmaEtaEta());
  S("full5x5_sigmaIetaIeta", e.full5x5_sigmaIetaIeta());
  S("full5x5_sigmaIphiIphi", e.full5x5_sigmaIphiIphi());
  S("full5x5_sigmaIetaIphi", e.full5x5_sigmaIetaIphi());
  S("full5x5_e1x5", e.full5x5_e1x5());
  S("full5x5_e2x5Max", e.full5x5_e2x5Max());
  S("full5x5_e5x5", e.full5x5_e5x5());
  S("full5x5_r9", e.full5x5_r9());
  S("full5x5_hcalDepth1OverEcal", e.full5x5_hcalOverEcal(1));
  S("full5x5_hcalDepth2OverEcal", e.full5x5_hcalOverEcal(2));
  S("full5x5_hcalDepth1OverEcalBc", e.full5x5_hcalOverEcalBc(1));
  S("full5x5_hcalDepth2OverEcalBc", e.full5x5_hcalOverEcalBc(2));
  S("full5x5_eLeft", e.full5x5_eLeft());
  S("full5x5_eRight", e.full5x5_eRight());
  S("full5x5_eBottom", e.full5x5_eBottom());
  S("full5x5_eTop", e.full5x5_eTop());
  S("full5x5_e2x5Left", e.full5x5_e2x5Left());
  S("full5x5_e2x5Right", e.full5x5_e2x5Right());
  S("full5x5_e2x5Bottom", e.full5x5_e2x5Bottom());
  S("full5x5_e2x5Top", e.full5x5_e2x5Top());
}

void TauTransformerIdONNX::fillMuon(float* out, const pat::Muon& m, const pat::Tau& tau) const {
  using namespace taut;
  Row r(featIdx_[2], out);
  const float pt = m.polarP4().pt();
  const float dphi = dPhi(float(m.polarP4().phi()), float(tau.polarP4().phi())), deta = float(m.polarP4().eta()) - float(tau.polarP4().eta());
  r.set("r", std::sqrt(deta * deta + dphi * dphi));
  r.set("theta", std::atan2(dphi, deta));
  r.set("rel_pt", pt / float(tau.polarP4().pt()));
  r.set("particle_type", 8.f);
  const float dxy = m.dB(pat::Muon::PV2D), dxyErr = m.edB(pat::Muon::PV2D);
  r.set("dxy", dxy);
  r.set("dxy_sig", std::abs(dxy) / dxyErr);
  const float chi2 = m.globalTrack().isNonnull() ? float(m.normChi2()) : kDefault;
  const bool chi2Valid = chi2 > 0 && isFiniteF(chi2);
  r.set("normalizedChi2_valid", float(chi2Valid));
  r.set("normalizedChi2", chi2Valid ? chi2 : kNaN);
  r.set("numberOfValidHits", chi2Valid ? (m.innerTrack().isNonnull() ? float(int(m.numberOfValidHits())) : kDefaultInt) : kNaN);
  r.set("segmentCompatibility", float(m.segmentCompatibility()));
  r.set("caloCompatibility", m.caloCompatibility());
  r.set("rel_pfEcalEnergy", m.pfEcalEnergy() / pt);
  const MuonStationCounts counts = countMuonStations(m);
  static const char* const subdets[3] = {"DT", "CSC", "RPC"};
  for (int d = 0; d < 3; ++d)
    for (int st = 1; st <= 4; ++st) {
      const std::string sfx = std::string(subdets[d]) + "_" + std::to_string(st);
      r.set("n_matches_" + sfx, float(counts.matches[d][st - 1]));
      r.set("n_hits_" + sfx, float(counts.hits[d][st - 1]));
    }
}

void TauTransformerIdONNX::fillGlobal(float* out, const pat::Tau& tau, float rho) const {
  using namespace taut;
  Row r(featIdx_[3], out);
  const float pt = tau.polarP4().pt(), eta = tau.polarP4().eta(), mass = tau.polarP4().mass();
  r.set("particle_type", 9.f);
  r.set("rho", rho);
  r.set("tau_pt", pt);
  r.set("tau_eta", eta);
  r.set("tau_mass", mass);
  r.set("tau_E_over_pt", std::sqrt((pt * std::cosh(eta)) * (pt * std::cosh(eta)) + mass * mass) / pt);
  r.set("tau_charge", float(tau.charge()));
  r.set("tau_n_charged_prongs", float(tau.decayMode() / 5 + 1));
  r.set("tau_n_neutral_prongs", float(tau.decayMode() % 5));
  const float chIso = tau.tauID("chargedIsoPtSum"), chIso03 = tau.tauID("chargedIsoPtSumdR03");
  const float nIso = tau.tauID("neutralIsoPtSum"), nIsoW = tau.tauID("neutralIsoPtSumWeight");
  const float nIsoW03 = tau.tauID("neutralIsoPtSumWeightdR03"), nIso03 = tau.tauID("neutralIsoPtSumdR03");
  r.set("tau_chargedIsoPtSum", chIso);
  r.set("tau_chargedIsoPtSumdR03_over_dR05", chIso != 0 ? chIso03 / chIso : kNaN);
  r.set("tau_footprintCorrection", tau.tauID("footprintCorrection"));
  r.set("tau_neutralIsoPtSum", nIso);
  r.set("tau_neutralIsoPtSumWeight_over_neutralIsoPtSum", nIso != 0 ? nIsoW / nIso : kNaN);
  r.set("tau_neutralIsoPtSumWeightdR03_over_neutralIsoPtSum", nIso != 0 ? nIsoW03 / nIso : kNaN);
  r.set("tau_neutralIsoPtSumdR03_over_dR05", nIso != 0 ? nIso03 / nIso : kNaN);
  r.set("tau_photonPtSumOutsideSignalCone", tau.tauID("photonPtSumOutsideSignalCone"));
  r.set("tau_puCorrPtSum", tau.tauID("puCorrPtSum"));
  r.set("tau_hasSecondaryVertex", float(tau.hasSecondaryVertex()));
  const float dxy = tau.dxy(), dxyErr = tau.dxy_error();
  const bool dxyValid = isFiniteF(dxy) && dxy > -10 && isFiniteF(dxyErr) && dxyErr > 0;
  r.set("tau_dxy_valid", float(dxyValid));
  r.set("tau_dxy", dxyValid ? dxy : kNaN);
  r.set("tau_dxy_sig", dxyValid ? std::abs(dxy) / dxyErr : kNaN);
  r.set("tau_ip3d", tau.ip3d());
  r.set("tau_ip3d_sig", std::abs(float(tau.ip3d())) / float(tau.ip3d_error()));
  const auto* lead = dynamic_cast<const pat::PackedCandidate*>(tau.leadChargedHadrCand().get());
  const float dz = lead ? lead->dz() : kDefault, dzErr = lead && lead->hasTrackDetails() ? lead->dzError() : kDefault;
  const bool dzValid = isFiniteF(dz) && isFiniteF(dzErr) && dzErr > 0;
  r.set("tau_dz_sig_valid", float(dzValid));
  r.set("tau_dz", dzValid ? dz : kNaN);
  r.set("tau_dz_sig", dzValid ? std::abs(dz) / dzErr : kNaN);
  r.set("tau_flightLength_x", tau.flightLength().x());
  r.set("tau_flightLength_y", tau.flightLength().y());
  r.set("tau_flightLength_z", tau.flightLength().z());
  r.set("tau_flightLength_sig", tau.flightLengthSig());
  r.set("tau_pt_weighted_deta_strip", reco::tau::pt_weighted_deta_strip(tau, tau.decayMode()));
  r.set("tau_pt_weighted_dphi_strip", reco::tau::pt_weighted_dphi_strip(tau, tau.decayMode()));
  r.set("tau_pt_weighted_dr_signal", reco::tau::pt_weighted_dr_signal(tau, tau.decayMode()));
  r.set("tau_pt_weighted_dr_iso", reco::tau::pt_weighted_dr_iso(tau, tau.decayMode()));
  r.set("tau_e_ratio", reco::tau::eratio(tau));
  // decay angle quantity
  float gj = kDefault;
  if (tau.hasSecondaryVertex()) {
    static constexpr double mTau = 1.77682;
    const double mA = tau.p4().M(), pA = tau.p();
    const double argMax = (mTau * mTau - mA * mA) / (2 * mTau * pA);
    const double argMeas = tau.p4().Vect().Dot(tau.flightLength()) / (pA * tau.flightLength().R());
    if (std::abs(argMax) <= 1. && std::abs(argMeas) <= 1.)
      gj = float(std::acos(argMeas) - std::asin(argMax));
  }
  const bool gjValid = gj >= 0;
  r.set("tau_gj_angle_diff_valid", float(gjValid));
  r.set("tau_gj_angle_diff", gjValid ? gj : kNaN);
  r.set("tau_n_photons", float(int(reco::tau::n_photons_total(tau))));
  r.set("tau_emFraction", tau.emFraction_MVA());
  const double aeta = std::abs(tau.polarP4().eta());
  r.set("tau_inside_ecal_crack", float(aeta > 1.46 && aeta < 1.558));
  r.set("tau_leadingTrackNormChi2", tau.leadingTrackNormChi2());
  r.set("tau_leadChargedCand_etaAtEcalEntrance_minus_tau_eta", tau.etaAtEcalEntranceLeadChargedCand() - eta);
}

void TauTransformerIdONNX::produce(edm::Event& event, const edm::EventSetup&) {
  const auto& taus = event.get(tausToken_);
  const auto& jets = event.get(jetsToken_);
  const auto& fatJets = event.get(fatJetsToken_);
  const auto& cands = event.get(candsToken_);
  const auto& electrons = event.get(electronsToken_);
  const auto& muons = event.get(muonsToken_);
  const auto& vertices = event.get(verticesToken_);
  const auto& svs = event.get(svToken_);
  const float rho = event.get(rhoToken_);
  static const pat::TauCollection noTaus;
  const pat::TauCollection* boosted = &noTaus;
  edm::Handle<pat::TauCollection> boostedHandle;
  if (useBoostedTaus_ && event.getByToken(boostedTausToken_, boostedHandle) && boostedHandle.isValid())
    boosted = boostedHandle.product();

  std::vector<size_t> sel;  // HPS taus
  for (size_t i = 0; i < taus.size(); ++i)
    if (taus[i].isTauIDAvailable("chargedIsoPtSum"))
      sel.push_back(i);
  const size_t nTau = sel.size();
  std::vector<std::vector<float>> scores(7, std::vector<float>(taus.size(), -1.f));
  if (nTau > 0 && !vertices.empty()) {
    const auto& pv = vertices.front();
    std::vector<Tokens> tok(nTau);
    size_t nPf = 1, nEle = 1, nMu = 1;
    for (size_t i = 0; i < nTau; ++i) {
      tok[i] = selectObjects(taus[sel[i]], *boosted, jets, fatJets, cands, electrons, muons, svs);
      nPf = std::max(nPf, tok[i].cands.size());
      nEle = std::max(nEle, tok[i].electrons.size());
      nMu = std::max(nMu, tok[i].muons.size());
    }
    cms::Ort::FloatArrays in(4);
    in[0].assign(nTau * nPf * nFeat_[0], 0.f);
    in[1].assign(nTau * nEle * nFeat_[1], 0.f);
    in[2].assign(nTau * nMu * nFeat_[2], 0.f);
    in[3].assign(nTau * nFeat_[3], 0.f);
    for (size_t i = 0; i < nTau; ++i) {
      for (size_t j = 0; j < tok[i].cands.size(); ++j)
        fillPfCand(&in[0][(i * nPf + j) * nFeat_[0]], cands[tok[i].cands[j]], taus[sel[i]], pv);
      for (size_t j = 0; j < tok[i].electrons.size(); ++j)
        fillElectron(&in[1][(i * nEle + j) * nFeat_[1]], electrons[tok[i].electrons[j]], taus[sel[i]]);
      for (size_t j = 0; j < tok[i].muons.size(); ++j)
        fillMuon(&in[2][(i * nMu + j) * nFeat_[2]], muons[tok[i].muons[j]], taus[sel[i]]);
      fillGlobal(&in[3][i * nFeat_[3]], taus[sel[i]], rho);
    }
    const int64_t B = nTau;
    const std::vector<std::vector<int64_t>> shapes = {{B, int64_t(nPf), int64_t(nFeat_[0])},
                                                      {B, int64_t(nEle), int64_t(nFeat_[1])},
                                                      {B, int64_t(nMu), int64_t(nFeat_[2])},
                                                      {B, int64_t(nFeat_[3])}};
    const auto out = globalCache()->session->run({"pfcand", "ele", "muon", "global"}, in, shapes, {"probs"}, B)[0];
    for (size_t i = 0; i < nTau; ++i) {
      const float* p = &out[i * 4];  // tau, e, mu, jet
      const size_t t = sel[i];
      for (int k = 0; k < 4; ++k)
        scores[k][t] = p[k];
      scores[4][t] = p[0] / (p[0] + p[3]);
      scores[5][t] = p[0] / (p[0] + p[1]);
      scores[6][t] = p[0] / (p[0] + p[2]);
    }
    if (dump_) {  // optional debug output
      for (size_t i = 0; i < nTau; ++i) {
        const uint64_t id[3] = {event.id().run(), event.id().luminosityBlock(), event.id().event()};
        const auto& tp = taus[sel[i]].polarP4();
        const float p4[3] = {float(tp.pt()), float(tp.eta()), float(tp.phi())};
        const int32_t n[4] = {int32_t(sel[i]), int32_t(tok[i].cands.size()), int32_t(tok[i].electrons.size()), int32_t(tok[i].muons.size())};
        dump_->write(reinterpret_cast<const char*>(id), sizeof(id));
        dump_->write(reinterpret_cast<const char*>(p4), sizeof(p4));
        dump_->write(reinterpret_cast<const char*>(n), sizeof(n));
        dump_->write(reinterpret_cast<const char*>(&in[0][i * nPf * nFeat_[0]]), sizeof(float) * n[1] * nFeat_[0]);
        dump_->write(reinterpret_cast<const char*>(&in[1][i * nEle * nFeat_[1]]), sizeof(float) * n[2] * nFeat_[1]);
        dump_->write(reinterpret_cast<const char*>(&in[2][i * nMu * nFeat_[2]]), sizeof(float) * n[3] * nFeat_[2]);
        dump_->write(reinterpret_cast<const char*>(&in[3][i * nFeat_[3]]), sizeof(float) * nFeat_[3]);
        dump_->write(reinterpret_cast<const char*>(&out[i * 4]), sizeof(float) * 4);
      }
    }
  }
  edm::Handle<pat::TauCollection> tauHandle;
  event.getByToken(tausToken_, tauHandle);
  // store the outputs
  for (int k = 0; k < 7; ++k) {
    std::vector<reco::SingleTauDiscriminatorContainer> buffer(taus.size());
    for (size_t t = 0; t < taus.size(); ++t)
      buffer[t].rawValues.push_back(scores[k][t]);
    auto out = std::make_unique<reco::TauDiscriminatorContainer>();
    reco::TauDiscriminatorContainer::Filler filler(*out);
    filler.insert(tauHandle, buffer.begin(), buffer.end());
    filler.fill();
    event.put(std::move(out), taut::kOutputs[k]);
  }
}

DEFINE_FWK_MODULE(TauTransformerIdONNX);
