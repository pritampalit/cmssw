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
#include <bitset>
// Xilinx HLS includes
#include <ap_fixed.h>
#include <ap_int.h>
#include <stdio.h>
#include <cassert>
#include <cstdlib>

// user include files
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/L1TCorrelator/interface/TkPhiCandidate.h"
#include "DataFormats/L1TCorrelator/interface/TkPhiCandidateFwd.h"
#include "DataFormats/L1Trigger/interface/TkLightMesonWord.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
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
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
//#include "hls_math.h"
//
// class declaration

//
using namespace std;
using namespace edm;
using namespace l1t;

class L1BsMesonSelectionEmulationProducer : public edm::global::EDProducer<> {
public:
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  explicit L1BsMesonSelectionEmulationProducer(const edm::ParameterSet&);
  ~L1BsMesonSelectionEmulationProducer() override;

  static constexpr double kmass = 0.493;
  double ETAPHI_LSB = M_PI / (1 << 12);
  double Z0_LSB = 0.05;


private:
  // ----------constants, enums and typedefs ---------
  // Relevant constants for the converted track word

  enum TrackBitWidths {
    kPtSize = TTTrack_TrackWord::TrackBitWidths::kRinvSize - 1,  // Width of pt                                                                                         
    kPtMagSize = 9,                                              // Width of pt magnitude (unsigned)                                                                    
    kEtaSize = TTTrack_TrackWord::TrackBitWidths::kTanlSize,     // Width of eta                                                                                        
    kEtaMagSize = 3,                                             // Width of eta magnitude (signed)                                                                     
  };
  
  typedef edm::Handle<TkLightMesonWordCollection> TkLightMesonWordCollectionHandle;

  // ----------member functions ----------------------
  /*  void printDebugInfo(const TTTrackCollectionHandle& l1PosKaonTracksHandle,
		      const TTTrackCollectionHandle& l1NegKaonTracksHandle,
                      const TTTrackRefCollectionUPtr& vTTTrackOutput,
                      const TTTrackRefCollectionUPtr& vTTTrackEmulationOutput) const;
		      void printTrackInfo(edm::LogInfo& log, const TTTrackRef& track, bool printEmulation = false) const;*/
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  //  void beginJob();
  //  void produce(edm::Event&, const edm::EventSetup&) override;
  //void endJob();

  // ----------selectors -----------------------------
  // Based on recommendations from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGenericSelectors

  // ----------member data ---------------------------
  const edm::EDGetTokenT<TkLightMesonWordCollection> l1PhiMesonWordToken_;
    
  const std::string outputCollectionName_;
  const edm::ParameterSet cutSet_;
  const double dRmax_, tkpairMmin_, tkpairMmax_;
  //bool processSimulatedTracks_, processEmulatedTracks_;
  int debug_;

};

//
// constructors and destructor
//
L1BsMesonSelectionEmulationProducer::L1BsMesonSelectionEmulationProducer(const edm::ParameterSet& iConfig)
  : l1PhiMesonWordToken_(consumes<TkLightMesonWordCollection>(iConfig.getParameter<edm::InputTag>("l1PhiMesonWordInputTag"))),
      outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")),
    cutSet_(iConfig.getParameter<edm::ParameterSet>("cutSet")),
    dRmax_(cutSet_.getParameter<double>("dRmax")),
    //    dzmax_(cutSet_.getParameter<double>("dzmax")),
    tkpairMmin_(cutSet_.getParameter<double>("tkpairMmin")),
    tkpairMmax_(cutSet_.getParameter<double>("tkpairMmax")),
      debug_(iConfig.getParameter<int>("debug")) {
  // Confirm the the configuration makes sense
  produces<l1t::TkLightMesonWordCollection>(outputCollectionName_);
  //produces<TkPhiCandidateRefVector>(outputCollectionName_);
}

L1BsMesonSelectionEmulationProducer::~L1BsMesonSelectionEmulationProducer() {}

//
// member functions
//
/*
void L1BsMesonSelectionEmulationProducer::printDebugInfo(const TTTrackCollectionHandle& l1PosKaonTracksHandle,
						 const TTTrackCollectionHandle& l1NegKaonTracksHandle,
                                              const TTTrackRefCollectionUPtr& vTTTrackOutput,
                                              const TTTrackRefCollectionUPtr& vTTTrackEmulationOutput) const {
  edm::LogInfo log("L1BsMesonSelectionEmulationProducer");
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

void L1BsMesonSelectionEmulationProducer::printTrackInfo(edm::LogInfo& log, const TTTrackRef& track, bool printEmulation) const {
  log << "\t(" << track->momentum().perp() << ", " << track->momentum().eta() << ", " << track->momentum().phi() << ", "
      << track->getStubRefs().size() << ", " << track->stubPtConsistency() << ", " << track->chi2ZRed() << ", "
      << track->chi2XYRed() << ", " << track->z0() << ")\n";

  if (printEmulation) {
    ap_uint<TTTrack_TrackWord::TrackBitWidths::kPtSize> ptEmulationBits = track->getTrackWord()(
        TTTrack_TrackWord::TrackBitLocations::kRinvMSB - 1, TTTrack_TrackWord::TrackBitLocations::kRinvLSB);
    ap_ufixed<TTTrack_TrackWord::TrackBitWidths::kPtSize, TTTrack_TrackWord::TrackBitWidths::kPtMagSize> ptEmulation;
    ptEmulation.V = ptEmulationBits.range();
    TTTrack_TrackWord::tanl_t etaEmulationBits = track->getTanlWord();
    ap_fixed<TTTrack_TrackWord::TrackBitWidths::kEtaSize, TTTrack_TrackWord::TrackBitWidths::kEtaMagSize> etaEmulation;
    etaEmulation.V = etaEmulationBits.range();
    log << "\t\t(" << ptEmulation.to_double() << ", " << etaEmulation.to_double() << ", " << track->getPhi() << ", "
        << track->getNStubs() << ", " << track->getBendChi2() << ", " << track->getChi2RZ() << ", " << track->getChi2RPhi()
        << ", " << track->getZ0() << ")\n";
  }
}
*/
// ------------ method called to produce the data  ------------
void L1BsMesonSelectionEmulationProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
//void L1BsMesonSelectionEmulationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  unique_ptr<l1t::TkLightMesonWordCollection> L1BsMesonEmulationOutput(new l1t::TkLightMesonWordCollection);
  //unique_ptr<TkPhiCandidateRefVector> L1PhiMesonOutputRefVec(new TkPhiCandidateRefVector);

  TkLightMesonWordCollectionHandle l1PhiMesonWordHandle;

  iEvent.getByToken(l1PhiMesonWordToken_, l1PhiMesonWordHandle);

  size_t nPhiMesonOutputApproximate = l1PhiMesonWordHandle->size();
  size_t nBsMesonOutputApproximate = 2*nPhiMesonOutputApproximate;

  L1BsMesonEmulationOutput->reserve(nBsMesonOutputApproximate);
  //L1PhiMesonOutputRefVec->reserve(nPhiMesonOutputApproximate);
  
  for (size_t i = 0; i < nPhiMesonOutputApproximate; i++) {
    const auto& tkPhiMesonWord1 = l1PhiMesonWordHandle->at(i);
        
    double trkptPhi1 = tkPhiMesonWord1.pt();
    double trketaPhi1 = tkPhiMesonWord1.glbeta();
    double trkphiPhi1 = tkPhiMesonWord1.glbphi();
    double trkz0Phi1 = tkPhiMesonWord1.z0();

    double trkpxPhi1 = trkptPhi1*cos(trkphiPhi1);
    double trkpyPhi1 = trkptPhi1*sin(trkphiPhi1);
    double trkpzPhi1 = trkptPhi1*sinh(trketaPhi1);
    
    for (size_t j = i+1; j < nPhiMesonOutputApproximate; j++) {
    const auto& tkPhiMesonWord2 = l1PhiMesonWordHandle->at(j);
        
    double trkptPhi2 = tkPhiMesonWord2.pt();
    double trketaPhi2 = tkPhiMesonWord2.glbeta();
    double trkphiPhi2 = tkPhiMesonWord2.glbphi();
    double trkz0Phi2 = tkPhiMesonWord2.z0();

    double trkpxPhi2 = trkptPhi2*cos(trkphiPhi2);
    double trkpyPhi2 = trkptPhi2*sin(trkphiPhi2);
    double trkpzPhi2 = trkptPhi2*sinh(trketaPhi2);
    
    double trkdrpairBs = sqrt(pow((trkphiPhi1 - trkphiPhi2),2) + pow((trketaPhi1 - trketaPhi2),2));
      // write mass calculation here , for hardware specially

    double trkmasspairBs = sqrt(2*trkptPhi1*trkptPhi2*(cosh(trketaPhi1 - trketaPhi2)-cos(trkphiPhi1 - trkphiPhi2)));

      if (trkdrpairBs > dRmax_) continue; 
      if (trkmasspairBs < tkpairMmin_ || trkmasspairBs > tkpairMmax_) continue; // do it before

      double trkpxBs = trkpxPhi2 + trkpxPhi1;
      double trkpyBs = trkpyPhi2 + trkpyPhi1;
      double trkpzBs = trkpzPhi2 + trkpzPhi1;
      
      l1t::TkLightMesonWord::valid_t trkvalidBs =   tkPhiMesonWord1.valid() && tkPhiMesonWord2.valid();
      l1t::TkLightMesonWord::pt_t trkptBs = sqrt(pow(trkpxBs,2) + pow(trkpyBs,2)); // use Pow()
      l1t::TkLightMesonWord::glbphi_t trkphiBs = atan(trkpyBs/trkpxBs);
      l1t::TkLightMesonWord::glbeta_t trketaBs = asinh(trkpzBs/sqrt(pow(trkpxBs,2) + pow(trkpyBs,2)));
      l1t::TkLightMesonWord::z0_t trkz0Bs = trkz0Phi1 + trkz0Phi2;
      l1t::TkLightMesonWord::mass_t trkmassBs = sqrt(2*trkptPhi1*trkptPhi2*(cosh(trketaPhi1 - trketaPhi2)-cos(trkphiPhi1 - trkphiPhi2)));
      l1t::TkLightMesonWord::type_t trktypeBs = l1t::TkLightMesonWord::TkLightMesonTypes::kBsType;
      l1t::TkLightMesonWord::ntracks_t trkntracksBs = 2;
      l1t::TkLightMesonWord::unassigned_t trkunassignedBs = 0;
      
      l1t::TkLightMesonWord trkBsWord(trkvalidBs, trkptBs, trkphiBs, trketaBs, trkz0Bs, trkmassBs, trktypeBs, trkntracksBs, trkunassignedBs);
      
      L1BsMesonEmulationOutput->push_back(trkBsWord);

      //      std::cout << __PRETTY_FUNCTION__ << __LINE__ << std::endl;
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
  iEvent.put(std::move(L1BsMesonEmulationOutput), outputCollectionName_);

}

//void L1BsMesonSelectionEmulationProducer::beginJob() {}

//void L1BsMesonSelectionEmulationProducer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1BsMesonSelectionEmulationProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //L1BsMesonSelectionEmulationProducer
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1PhiMesonWordInputTag", edm::InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"));
  desc.add<std::string>("outputCollectionName", "Level1TTKaonTracksSelected");
  {
    edm::ParameterSetDescription descCutSet;
    descCutSet.add<double>("dRmax", 0.12)->setComment("dr must be less than this value, []");
    //    descCutSet.add<double>("dxymax", 1.0)->setComment("dxy must be less than this value, [cm]");
    //descCutSet.add<double>("dzmax", 1.0)->setComment("dz must be less than this value, [cm]");
    descCutSet.add<double>("tkpairMmin", 1.0)->setComment("tkpair mass must be greater than this value, [GeV]");
    descCutSet.add<double>("tkpairMmax", 1.03)->setComment("tkpair mass must be less than this value, [GeV]");
    desc.add<edm::ParameterSetDescription>("cutSet", descCutSet);

  }
  desc.add<int>("debug", 0)->setComment("Verbosity levels: 0, 1, 2, 3");
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1BsMesonSelectionEmulationProducer);
