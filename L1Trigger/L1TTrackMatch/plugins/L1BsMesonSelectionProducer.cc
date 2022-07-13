// -*- C++ -*-
//
// Package:    L1Trigger/L1TTrackMatch
// Class:      L1BsMesonSelectionProducer
//
/**\class L1BsMesonSelectionProducer L1BsMesonSelectionProducer.cc L1Trigger/L1TTrackMatch/plugins/L1BsMesonSelectionProducer.cc

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
#include "DataFormats/L1TCorrelator/interface/TkBsCandidate.h"
#include "DataFormats/L1TCorrelator/interface/TkPhiCandidateFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkBsCandidateFwd.h"
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

class L1BsMesonSelectionProducer : public edm::global::EDProducer<> {
public:
  explicit L1BsMesonSelectionProducer(const edm::ParameterSet&);
  ~L1BsMesonSelectionProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static constexpr double kmass = 0.493;

private:
  // ----------constants, enums and typedefs ---------
  // Relevant constants for the converted track word

  //typedef std::vector<TkPhiCandidate> TkPhiCandidateCollection;
  //typedef edm::Ref<TkPhiCandidateCollection> TkPhiCandidateRef;
  //typedef edm::RefVector<TkPhiCandidateCollection> TkPhiCandidateRefCollection;
  //typedef edm::Handle<TkPhiCandidateRefCollection> TkPhiCandidateCollectionHandle;
  //typedef std::unique_ptr<TkPhiCandidateRefCollection> TkPhiCandidateRefCollectionUPtr;

  typedef edm::Handle<TkPhiCandidateCollection> TkPhiCandidateCollectionHandle;
  //typedef edm::Handle<TkPhiCandidateRefVector> TkPhiCandidateCollectionHandle;
  

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
  const edm::EDGetTokenT<TkPhiCandidateCollection> l1PhiCandsToken_;
  const std::string outputCollectionName_;
  const edm::ParameterSet cutSet_;
  const double dRmax_, dRmin_, dxymax_, dzmax_, phipairMmin_, phipairMmax_;
  //bool processSimulatedTracks_, processEmulatedTracks_;
  int debug_;

};

//
// constructors and destructor
//
L1BsMesonSelectionProducer::L1BsMesonSelectionProducer(const edm::ParameterSet& iConfig)
  : l1PhiCandsToken_(consumes<TkPhiCandidateCollection>(iConfig.getParameter<edm::InputTag>("l1PhiCandsInputTag"))),
    outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")),
    cutSet_(iConfig.getParameter<edm::ParameterSet>("cutSet")),
    dRmax_(cutSet_.getParameter<double>("dRmax")),
    dRmin_(cutSet_.getParameter<double>("dRmin")),
    dxymax_(cutSet_.getParameter<double>("dxymax")),
    dzmax_(cutSet_.getParameter<double>("dzmax")),
    phipairMmin_(cutSet_.getParameter<double>("phipairMmin")),
    phipairMmax_(cutSet_.getParameter<double>("phipairMmax")),
  debug_(iConfig.getParameter<int>("debug")) {
  // Confirm the the configuration makes sense
  produces<TkBsCandidateCollection>(outputCollectionName_);
  //produces<TkBsCandidateRefVector>(outputCollectionName_);
}

L1BsMesonSelectionProducer::~L1BsMesonSelectionProducer() {}

//
// member functions
//
/*
void L1BsMesonSelectionProducer::printDebugInfo(const TTTrackCollectionHandle& l1PosKaonTracksHandle,
						 const TTTrackCollectionHandle& l1NegKaonTracksHandle,
                                              const TTTrackRefCollectionUPtr& vTTTrackOutput,
                                              const TTTrackRefCollectionUPtr& vTTTrackEmulationOutput) const {
  edm::LogInfo log("L1BsMesonSelectionProducer");
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

void L1BsMesonSelectionProducer::printTrackInfo(edm::LogInfo& log, const TTTrackRef& track, bool printEmulation) const {
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
void L1BsMesonSelectionProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  unique_ptr<TkBsCandidateCollection> L1BsMesonOutput(new TkBsCandidateCollection);
  //unique_ptr<TkBsCandidateRefVector> L1BsMesonOutputRef(new TkBsCandidateRefVector);


  TkPhiCandidateCollectionHandle l1PhiCandsHandle;

  iEvent.getByToken(l1PhiCandsToken_, l1PhiCandsHandle);
  size_t nPhiMesonOutputApproximate = l1PhiCandsHandle->size();
  size_t nBsMesonOutputApproximate = 2*nPhiMesonOutputApproximate;

  L1BsMesonOutput->reserve(nBsMesonOutputApproximate);
  
  for (size_t i = 0; i < nPhiMesonOutputApproximate; i++) {
    const auto& trackPhiCands1 = l1PhiCandsHandle->at(i);

    for (size_t j = i+1; j < nPhiMesonOutputApproximate; j++) {
      const auto& trackPhiCands2 = l1PhiCandsHandle->at(j);

      //    const edm::Ptr<L1Track>& trackPosKaonReftoPtr = edm::refToPtr(trackPosKaonRef);
      //const edm::Ptr<L1Track>& trackNegKaonReftoPtr = edm::refToPtr(trackNegKaonRef);

      float l1Phi1pt = trackPhiCands1.p4().Pt();
    float l1Phi1eta = trackPhiCands1.p4().Eta();
    float l1Phi1phi = trackPhiCands1.p4().Phi();
    float l1Phi1px = l1Phi1pt*cos(l1Phi1phi);
    float l1Phi1py = l1Phi1pt*sin(l1Phi1phi);
    float l1Phi1pz = l1Phi1pt*sinh(l1Phi1eta);
    float l1Phi1e = l1Phi1pt*cosh(l1Phi1eta);

    math::XYZTLorentzVector Phi1P4(l1Phi1px, l1Phi1py, l1Phi1pz, l1Phi1e);

    float l1Phi2pt = trackPhiCands2.p4().Pt();
    float l1Phi2eta = trackPhiCands2.p4().Eta();
    float l1Phi2phi = trackPhiCands2.p4().Phi();
    float l1Phi2px = l1Phi2pt*cos(l1Phi2phi);
    float l1Phi2py = l1Phi2pt*sin(l1Phi2phi);
    float l1Phi2pz = l1Phi2pt*sinh(l1Phi2eta);
    float l1Phi2e = l1Phi2pt*cosh(l1Phi2eta);

    math::XYZTLorentzVector Phi2P4(l1Phi2px, l1Phi2py, l1Phi2pz, l1Phi2e);

    

    TkBsCandidate tkBs(Phi1P4 + Phi2P4, trackPhiCands1, trackPhiCands2);
    
    if (tkBs.dxyPhiPair() > dxymax_) continue;
    if (std::fabs(tkBs.dzPhiPair()) > dzmax_) continue;
    if (tkBs.dRPhiPair() > dRmax_) continue;
    if (tkBs.dRPhiPair() < dRmin_) continue;
    //std::cout << "phi mass : " << tkBs.p4().M() << std::endl;
    if (tkBs.p4().M() < phipairMmin_ || tkBs.p4().M() > phipairMmax_) continue;

    L1BsMesonOutput->push_back(tkBs);
    }
  }

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
  iEvent.put(std::move(L1BsMesonOutput), outputCollectionName_);

  

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1BsMesonSelectionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //L1BsMesonSelectionProducer
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1PhiCandsInputTag", edm::InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"));
  desc.add<std::string>("outputCollectionName", "Level1TTKaonTracksSelected");
  {
    edm::ParameterSetDescription descCutSet;
    descCutSet.add<double>("dRmax", 0.12)->setComment("dr must be less than this value, []");
    descCutSet.add<double>("dRmin", 0.12)->setComment("dr must be greater than this value, []");
    descCutSet.add<double>("dxymax", 1.0)->setComment("dxy must be less than this value, [cm]");
    descCutSet.add<double>("dzmax", 1.0)->setComment("dz must be less than this value, [cm]");
    descCutSet.add<double>("phipairMmin", 1.0)->setComment("phipair mass must be greater than this value, [GeV]");
    descCutSet.add<double>("phipairMmax", 1.03)->setComment("phipair mass must be less than this value, [GeV]");
    desc.add<edm::ParameterSetDescription>("cutSet", descCutSet);

  }
  desc.add<int>("debug", 0)->setComment("Verbosity levels: 0, 1, 2, 3");
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1BsMesonSelectionProducer);
