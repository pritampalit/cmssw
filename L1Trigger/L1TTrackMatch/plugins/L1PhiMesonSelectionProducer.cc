// -*- C++ -*-
//
// Package:    L1Trigger/L1TTrackMatch
// Class:      L1PhiMesonSelectionProducer
//
/**\class L1PhiMesonSelectionProducer L1PhiMesonSelectionProducer.cc L1Trigger/L1TTrackMatch/plugins/L1PhiMesonSelectionProducer.cc

 Description: Selects two set of positively and negatively charged L1Tracks corresponding to Kaons which already passed the criteria for for Light Meson track selection

 Implementation:
     Inputs:
         std::vector<TTTrack> - Each floating point TTTrack inside this collection inherits from
                                a bit-accurate TTTrack_TrackWord, used for emulation purposes.
     Outputs:
         std::vector<TTTrack> - A collection of TTTracks selected from cuts on the TTTrack properties
         std::vector<TTTrack> - A collection of TTTracks selected from cuts on the TTTrack_TrackWord properties
*/
//
// Original Author:  Alexx Perloff
//         Created:  Thu, 16 Dec 2021 19:02:50 GMT
//
//

// system include files
#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <TMath.h>
#include <cmath>

// Xilinx HLS includes
#include <ap_fixed.h>
#include <ap_int.h>

// user include files
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/L1TCorrelator/interface/TkPhiCandidate.h"
#include "DataFormats/L1TCorrelator/interface/TkPhiCandidateFwd.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1Trigger/interface/Vertex.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "CommonTools/Utils/interface/AndSelector.h"
#include "CommonTools/Utils/interface/EtaRangeSelector.h"
#include "CommonTools/Utils/interface/MinSelector.h"
#include "CommonTools/Utils/interface/MinFunctionSelector.h"
#include "CommonTools/Utils/interface/MinNumberSelector.h"
#include "CommonTools/Utils/interface/PtMinSelector.h"
#include "CommonTools/Utils/interface/Selection.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/Math/interface/LorentzVector.h"

//
// class declaration

//
using namespace std;
using namespace edm;
using namespace l1t;

class L1PhiMesonSelectionProducer : public edm::global::EDProducer<> {
public:
  explicit L1PhiMesonSelectionProducer(const edm::ParameterSet&);
  ~L1PhiMesonSelectionProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static constexpr double kmass = 0.493;

private:
  // ----------constants, enums and typedefs ---------
  // Relevant constants for the converted track word

  typedef TTTrack<Ref_Phase2TrackerDigi_> L1Track;
  typedef std::vector<L1Track> TTTrackCollection;
  typedef edm::Ref<TTTrackCollection> TTTrackRef;
  typedef edm::RefVector<TTTrackCollection> TTTrackRefCollection;
  typedef edm::Handle<TTTrackRefCollection> TTTrackCollectionHandle;
  typedef std::unique_ptr<TTTrackRefCollection> TTTrackRefCollectionUPtr;


  // ----------member functions ----------------------
  /*  void printDebugInfo(const TTTrackCollectionHandle& l1PosKaonTracksHandle,
		      const TTTrackCollectionHandle& l1NegKaonTracksHandle,
                      const TTTrackRefCollectionUPtr& vTTTrackOutput,
                      const TTTrackRefCollectionUPtr& vTTTrackEmulationOutput) const;
		      void printTrackInfo(edm::LogInfo& log, const TTTrackRef& track, bool printEmulation = false) const;*/
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  // ----------selectors -----------------------------
  // Based on recommendations from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGenericSelectors

  // ----------member data ---------------------------
  const edm::EDGetTokenT<TTTrackRefCollection> l1PosKaonTracksToken_;
  const edm::EDGetTokenT<TTTrackRefCollection> l1NegKaonTracksToken_;
  const std::string outputCollectionName_;
  const edm::ParameterSet cutSet_;
  const double dRmax_, dxymax_, dzmax_, tkpairMmin_, tkpairMmax_;
  //bool processSimulatedTracks_, processEmulatedTracks_;
  int debug_;

};

//
// constructors and destructor
//
L1PhiMesonSelectionProducer::L1PhiMesonSelectionProducer(const edm::ParameterSet& iConfig)
  : l1PosKaonTracksToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("l1PosKaonTracksInputTag"))),
    l1NegKaonTracksToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("l1NegKaonTracksInputTag"))),
      outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")),
    cutSet_(iConfig.getParameter<edm::ParameterSet>("cutSet")),
    dRmax_(cutSet_.getParameter<double>("dRmax")),
    dxymax_(cutSet_.getParameter<double>("dxymax")),
    dzmax_(cutSet_.getParameter<double>("dzmax")),
    tkpairMmin_(cutSet_.getParameter<double>("tkpairMmin")),
    tkpairMmax_(cutSet_.getParameter<double>("tkpairMmax")),
      debug_(iConfig.getParameter<int>("debug")) {
  // Confirm the the configuration makes sense
  produces<TkPhiCandidateCollection>(outputCollectionName_);
  //produces<TkPhiCandidateRefVector>(outputCollectionName_);
}

L1PhiMesonSelectionProducer::~L1PhiMesonSelectionProducer() {}

//
// member functions
//
/*
void L1PhiMesonSelectionProducer::printDebugInfo(const TTTrackCollectionHandle& l1PosKaonTracksHandle,
						 const TTTrackCollectionHandle& l1NegKaonTracksHandle,
                                              const TTTrackRefCollectionUPtr& vTTTrackOutput,
                                              const TTTrackRefCollectionUPtr& vTTTrackEmulationOutput) const {
  edm::LogInfo log("L1PhiMesonSelectionProducer");
  log << "The original Positive Kaon track collection (pt, eta, phi, nstub, bendchi2, chi2rz, chi2rphi, z0) values are ... \n";
  for (const auto& track : *l1PosKaonTracksHandle) {
    printTrackInfo(log, track, debug_ >= 4);
  }
  log << "\t---\n\tNumber of Positive Kaon tracks in this selection = " << l1PosKaonTracksHandle->size() << "\n\n";
  for (const auto& track : *l1NegKaonTracksHandle) {
    printTrackInfo(log, track, debug_ >= 4);
  }
  log << "\t---\n\tNumber of Negative tracks in this selection = " << l1NegKaonTracksHandle->size() << "\n\n";
  if (processSimulatedTracks_) {
    log << "The selected phi collection (pt, eta, phi, nstub, bendchi2, chi2rz, chi2rphi, z0) values are ... \n";
    for (const auto& track : *vTTTrackOutput) {
      printTrackInfo(log, track, debug_ >= 4);
    }
    log << "\t---\n\tNumber of tracks in this selection = " << vTTTrackOutput->size() << "\n\n";
  }
  if (processEmulatedTracks_) {
    log << "The emulation selected track collection (pt, eta, phi, nstub, bendchi2, chi2rz, chi2rphi, z0) values are "
           "... \n";
    for (const auto& track : *vTTTrackEmulationOutput) {
      printTrackInfo(log, track, debug_ >= 4);
    }
    log << "\t---\n\tNumber of tracks in this selection = " << vTTTrackEmulationOutput->size() << "\n\n";
  }
  if (processSimulatedTracks_ && processEmulatedTracks_) {
    TTTrackRefCollection inSimButNotEmu;
    TTTrackRefCollection inEmuButNotSim;
    std::set_difference(vTTTrackOutput->begin(),
                        vTTTrackOutput->end(),
                        vTTTrackEmulationOutput->begin(),
                        vTTTrackEmulationOutput->end(),
                        std::back_inserter(inSimButNotEmu));
    std::set_difference(vTTTrackEmulationOutput->begin(),
                        vTTTrackEmulationOutput->end(),
                        vTTTrackOutput->begin(),
                        vTTTrackOutput->end(),
                        std::back_inserter(inEmuButNotSim));
    log << "The set of tracks selected via cuts on the simulated values which are not in the set of tracks selected "
           "by cutting on the emulated values ... \n";
    for (const auto& track : inSimButNotEmu) {
      printTrackInfo(log, track, debug_ >= 3);
    }
    log << "\t---\n\tNumber of tracks in this selection = " << inSimButNotEmu.size() << "\n\n"
        << "The set of tracks selected via cuts on the emulated values which are not in the set of tracks selected "
           "by cutting on the simulated values ... \n";
    for (const auto& track : inEmuButNotSim) {
      printTrackInfo(log, track, debug_ >= 3);
    }
    log << "\t---\n\tNumber of tracks in this selection = " << inEmuButNotSim.size() << "\n\n";
  }
}

void L1PhiMesonSelectionProducer::printTrackInfo(edm::LogInfo& log, const TTTrackRef& track, bool printEmulation) const {
  log << "\t(" << track->momentum().perp() << ", " << track->momentum().eta() << ", " << track->momentum().phi() << ", "
      << track->getStubRefs().size() << ", " << track->stubPtConsistency() << ", " << track->chi2ZRed() << ", "
      << track->chi2XYRed() << ", " << track->z0() << ")\n";

  if (printEmulation) {
    ap_uint<TrackBitWidths::kPtSize> ptEmulationBits = track->getTrackWord()(
        TTTrack_TrackWord::TrackBitLocations::kRinvMSB - 1, TTTrack_TrackWord::TrackBitLocations::kRinvLSB);
    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize> ptEmulation;
    ptEmulation.V = ptEmulationBits.range();
    TTTrack_TrackWord::tanl_t etaEmulationBits = track->getTanlWord();
    ap_fixed<TrackBitWidths::kEtaSize, TrackBitWidths::kEtaMagSize> etaEmulation;
    etaEmulation.V = etaEmulationBits.range();
    log << "\t\t(" << ptEmulation.to_double() << ", " << etaEmulation.to_double() << ", " << track->getPhi() << ", "
        << track->getNStubs() << ", " << track->getBendChi2() << ", " << track->getChi2RZ() << ", " << track->getChi2RPhi()
        << ", " << track->getZ0() << ")\n";
  }
}
*/
// ------------ method called to produce the data  ------------
void L1PhiMesonSelectionProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  unique_ptr<TkPhiCandidateCollection> L1PhiMesonOutput(new TkPhiCandidateCollection);
  //unique_ptr<TkPhiCandidateRefVector> L1PhiMesonOutputRefVec(new TkPhiCandidateRefVector);

  TTTrackCollectionHandle l1PosKaonTracksHandle;
  TTTrackCollectionHandle l1NegKaonTracksHandle;

  iEvent.getByToken(l1PosKaonTracksToken_, l1PosKaonTracksHandle);
  iEvent.getByToken(l1NegKaonTracksToken_, l1NegKaonTracksHandle);
  size_t nPosKaonOutputApproximate = l1PosKaonTracksHandle->size();
  size_t nNegKaonOutputApproximate = l1NegKaonTracksHandle->size();
  size_t nPhiMesonOutputApproximate = nPosKaonOutputApproximate + nNegKaonOutputApproximate;

  L1PhiMesonOutput->reserve(nPhiMesonOutputApproximate);
  //L1PhiMesonOutputRefVec->reserve(nPhiMesonOutputApproximate);
  
  for (size_t i = 0; i < nPosKaonOutputApproximate; i++) {
    const auto& trackPosKaonRef = l1PosKaonTracksHandle->at(i);
    const auto& trackPosKaon = *trackPosKaonRef;

    for (size_t j = 0; j < nNegKaonOutputApproximate; j++) {
    const auto& trackNegKaonRef = l1NegKaonTracksHandle->at(j);
    const auto& trackNegKaon = *trackNegKaonRef;

    const edm::Ptr<L1Track>& trackPosKaonReftoPtr = edm::refToPtr(trackPosKaonRef);
    const edm::Ptr<L1Track>& trackNegKaonReftoPtr = edm::refToPtr(trackNegKaonRef);

    float l1postkpt = trackPosKaon.momentum().perp();
    float l1postketa = trackPosKaon.momentum().eta();
    float l1postkphi = trackPosKaon.momentum().phi();
    float l1postkpx = l1postkpt*cos(l1postkphi);
    float l1postkpy = l1postkpt*sin(l1postkphi);
    float l1postkpz = l1postkpt*sinh(l1postketa);
    float l1postke = l1postkpt*cosh(l1postketa);

    math::XYZTLorentzVector PosKaonP4(l1postkpx, l1postkpy, l1postkpz, l1postke);

    float l1negtkpt = trackNegKaon.momentum().perp();
    float l1negtketa = trackNegKaon.momentum().eta();
    float l1negtkphi = trackNegKaon.momentum().phi();
    float l1negtkpx = l1negtkpt*cos(l1negtkphi);
    float l1negtkpy = l1negtkpt*sin(l1negtkphi);
    float l1negtkpz = l1negtkpt*sinh(l1negtketa);
    float l1negtke = l1negtkpt*cosh(l1negtketa);

    math::XYZTLorentzVector NegKaonP4(l1negtkpx, l1negtkpy, l1negtkpz, l1negtke);

    TkPhiCandidate tkphi(PosKaonP4 + NegKaonP4, trackPosKaonReftoPtr, trackNegKaonReftoPtr);
    
    if (tkphi.dxyTrkPair() > dxymax_) continue;
    if (std::fabs(tkphi.dzTrkPair()) > dzmax_) continue;
    if (tkphi.dRTrkPair() > dRmax_) continue;
    // std::cout << "phi mass : " << tkphi.p4().M() << std::endl;
    if (tkphi.p4().M() < tkpairMmin_ || tkphi.p4().M() > tkpairMmax_) continue;

    L1PhiMesonOutput->push_back(tkphi);
    }
  }

  //for (size_t iphi = 0; iphi < L1PhiMesonOutput->size(); iphi++) {
    //std::cout << "tkphi M : " << L1PhiMesonOutput->at(iphi).p4().M() << std::endl;
    //L1PhiMesonOutputRefVec->push_back(TkPhiCandidateRef(L1PhiMesonOutput, *iphi));
  //}

  //for (std::vector<TkPhiCandidate>::iterator it = L1PhiMesonOutput->begin(); it != L1PhiMesonOutput->end(); ++it) {
  //std::cout << "tkphi M : " << it->p4().M() << std::endl;
  //L1PhiMesonOutputRefVec->push_back(TkPhiCandidateRef(L1PhiMesonOutput, *it));
  //}

  //TkPhiCandidateRef L1PhiMesonOutputRef = L1PhiMesonOutput;


    /*
  if (debug_ >= 2) {
    printDebugInfo(l1TracksHandle,
                   vTTPosTrackOutput,
                   vTTPosTrackEmulationOutput);

    printDebugInfo(l1TracksHandle,
                   vTTNegTrackOutput,
                   vTTNegTrackEmulationOutput);
		   }*/

  // Put the outputs into the event
  //L1PhiMesonOutputRefVec = TkPhiCandidateRefVector(L1PhiMesonOutput);
  iEvent.put(std::move(L1PhiMesonOutput), outputCollectionName_);

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1PhiMesonSelectionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //L1PhiMesonSelectionProducer
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1PosKaonTracksInputTag", edm::InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"));
  desc.add<edm::InputTag>("l1NegKaonTracksInputTag", edm::InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"));
  desc.add<std::string>("outputCollectionName", "Level1TTKaonTracksSelected");
  {
    edm::ParameterSetDescription descCutSet;
    descCutSet.add<double>("dRmax", 0.12)->setComment("dr must be less than this value, []");
    descCutSet.add<double>("dxymax", 1.0)->setComment("dxy must be less than this value, [cm]");
    descCutSet.add<double>("dzmax", 1.0)->setComment("dz must be less than this value, [cm]");
    descCutSet.add<double>("tkpairMmin", 1.0)->setComment("tkpair mass must be greater than this value, [GeV]");
    descCutSet.add<double>("tkpairMmax", 1.03)->setComment("tkpair mass must be less than this value, [GeV]");
    desc.add<edm::ParameterSetDescription>("cutSet", descCutSet);

  }
  desc.add<int>("debug", 0)->setComment("Verbosity levels: 0, 1, 2, 3");
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1PhiMesonSelectionProducer);
