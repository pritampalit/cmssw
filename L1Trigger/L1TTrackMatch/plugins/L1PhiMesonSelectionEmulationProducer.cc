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

class L1PhiMesonSelectionEmulationProducer : public edm::global::EDProducer<> {
public:
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  explicit L1PhiMesonSelectionEmulationProducer(const edm::ParameterSet&);
  ~L1PhiMesonSelectionEmulationProducer() override;

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
  //  void beginJob();
  //  void produce(edm::Event&, const edm::EventSetup&) override;
  //void endJob();

  // ----------selectors -----------------------------
  // Based on recommendations from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGenericSelectors

  // ----------member data ---------------------------
  const edm::EDGetTokenT<TTTrackRefCollection> l1PosKaonTracksToken_;
  const edm::EDGetTokenT<TTTrackRefCollection> l1NegKaonTracksToken_;
  const std::string outputCollectionName_;
  const edm::ParameterSet cutSet_;
  const double dRmax_, tkpairMmin_, tkpairMmax_;
  //bool processSimulatedTracks_, processEmulatedTracks_;
  int debug_;

};

//
// constructors and destructor
//
L1PhiMesonSelectionEmulationProducer::L1PhiMesonSelectionEmulationProducer(const edm::ParameterSet& iConfig)
  : l1PosKaonTracksToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("l1PosKaonTracksInputTag"))),
    l1NegKaonTracksToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("l1NegKaonTracksInputTag"))),
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

L1PhiMesonSelectionEmulationProducer::~L1PhiMesonSelectionEmulationProducer() {}

//
// member functions
//
/*
void L1PhiMesonSelectionEmulationProducer::printDebugInfo(const TTTrackCollectionHandle& l1PosKaonTracksHandle,
						 const TTTrackCollectionHandle& l1NegKaonTracksHandle,
                                              const TTTrackRefCollectionUPtr& vTTTrackOutput,
                                              const TTTrackRefCollectionUPtr& vTTTrackEmulationOutput) const {
  edm::LogInfo log("L1PhiMesonSelectionEmulationProducer");
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

void L1PhiMesonSelectionEmulationProducer::printTrackInfo(edm::LogInfo& log, const TTTrackRef& track, bool printEmulation) const {
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
void L1PhiMesonSelectionEmulationProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
//void L1PhiMesonSelectionEmulationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  unique_ptr<l1t::TkLightMesonWordCollection> L1PhiMesonEmulationOutput(new l1t::TkLightMesonWordCollection);
  //unique_ptr<TkPhiCandidateRefVector> L1PhiMesonOutputRefVec(new TkPhiCandidateRefVector);

  TTTrackCollectionHandle l1PosKaonTracksHandle;
  TTTrackCollectionHandle l1NegKaonTracksHandle;

  iEvent.getByToken(l1PosKaonTracksToken_, l1PosKaonTracksHandle);
  iEvent.getByToken(l1NegKaonTracksToken_, l1NegKaonTracksHandle);
  size_t nPosKaonOutputApproximate = l1PosKaonTracksHandle->size();
  size_t nNegKaonOutputApproximate = l1NegKaonTracksHandle->size();
  size_t nPhiMesonOutputApproximate = nPosKaonOutputApproximate + nNegKaonOutputApproximate;

  L1PhiMesonEmulationOutput->reserve(nPhiMesonOutputApproximate);
  //L1PhiMesonOutputRefVec->reserve(nPhiMesonOutputApproximate);
  
  for (size_t i = 0; i < nPosKaonOutputApproximate; i++) {
    const auto& trackPosKaonRef = l1PosKaonTracksHandle->at(i);
    const auto& trackPosKaon = *trackPosKaonRef;
    
    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> inputTrkPtPos;
    inputTrkPtPos.V = trackPosKaon.getTrackWord()(TTTrack_TrackWord::TrackBitLocations::kRinvMSB - 1,TTTrack_TrackWord::TrackBitLocations::kRinvLSB);
    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> trkptPos = inputTrkPtPos;

    ap_int<TrackBitWidths::kEtaSize> trketainputPos;
    trketainputPos.V = trackPosKaon.getTrackWord()(TTTrack_TrackWord::TrackBitLocations::kTanlMSB, TTTrack_TrackWord::TrackBitLocations::kTanlLSB);
    ap_ufixed<64, 32> etaphi_conv = 1.0 / ETAPHI_LSB;
    ap_int<TrackBitWidths::kEtaSize> trketaPos = etaphi_conv * trketainputPos;
    //ap_int<TrackBitWidths::kEtaSize> trketaPos = trketainputPos / 2.;

    ap_int<TTTrack_TrackWord::TrackBitWidths::kPhiSize> trkphiinputPos;
    trkphiinputPos.V = trackPosKaon.getTrackWord()(TTTrack_TrackWord::TrackBitLocations::kPhiMSB, TTTrack_TrackWord::TrackBitLocations::kPhiLSB);
    ap_int<TTTrack_TrackWord::TrackBitWidths::kPhiSize> trkphiPos = etaphi_conv * trkphiinputPos;

    ap_int<TTTrack_TrackWord::TrackBitWidths::kZ0Size> trkz0inputPos;
    trkz0inputPos.V = trackPosKaon.getTrackWord()(TTTrack_TrackWord::TrackBitLocations::kZ0MSB, TTTrack_TrackWord::TrackBitLocations::kZ0LSB);
    ap_ufixed<64, 32> z0_conv = 1.0 / Z0_LSB;
    ap_int<TTTrack_TrackWord::TrackBitWidths::kZ0Size> trkz0Pos = z0_conv * trkz0inputPos;
    //ap_int<TTTrack_TrackWord::TrackBitWidths::kZ0Size> trkz0Pos = trkz0inputPos / 0.05;

    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> trkpxPos = trkptPos.to_double()*cos(trkphiPos.to_double());
    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> trkpyPos = trkptPos.to_double()*sin(trkphiPos.to_double());
    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> trkpzPos = trkptPos.to_double()*sinh(trketaPos.to_double());
    
    for (size_t j = 0; j < nNegKaonOutputApproximate; j++) {
      const auto& trackNegKaonRef = l1NegKaonTracksHandle->at(j);
      const auto& trackNegKaon = *trackNegKaonRef;
      
      ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> inputTrkPtNeg;
      inputTrkPtNeg.V = trackNegKaon.getTrackWord()(TTTrack_TrackWord::TrackBitLocations::kRinvMSB - 1,TTTrack_TrackWord::TrackBitLocations::kRinvLSB);
      ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> trkptNeg = inputTrkPtNeg;
      
      ap_int<TrackBitWidths::kEtaSize> trketainputNeg;
      trketainputNeg.V = trackNegKaon.getTrackWord()(TTTrack_TrackWord::TrackBitLocations::kTanlMSB, TTTrack_TrackWord::TrackBitLocations::kTanlLSB);
      //    ap_ufixed<64, 32> etaphi_conv = 1.0 / ETAPHI_LSB;
      ap_int<TrackBitWidths::kEtaSize> trketaNeg = etaphi_conv * trketainputNeg;
      
      ap_int<TTTrack_TrackWord::TrackBitWidths::kPhiSize> trkphiinputNeg;
      trkphiinputNeg.V = trackNegKaon.getTrackWord()(TTTrack_TrackWord::TrackBitLocations::kPhiMSB, TTTrack_TrackWord::TrackBitLocations::kPhiLSB);
      ap_int<TTTrack_TrackWord::TrackBitWidths::kPhiSize> trkphiNeg = etaphi_conv * trkphiinputNeg;
      
      ap_int<TTTrack_TrackWord::TrackBitWidths::kZ0Size> trkz0inputNeg;
      trkz0inputNeg.V = trackNegKaon.getTrackWord()(TTTrack_TrackWord::TrackBitLocations::kZ0MSB, TTTrack_TrackWord::TrackBitLocations::kZ0LSB);
      //ap_ufixed<64, 32> z0_conv = 1.0 / Z0_LSB;
      ap_int<TTTrack_TrackWord::TrackBitWidths::kZ0Size> trkz0Neg = z0_conv * trkz0inputNeg;
      
      ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> trkpxNeg = trkptNeg.to_double()*cos(trkphiNeg.to_double());
      ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> trkpyNeg = trkptNeg.to_double()*sin(trkphiNeg.to_double());
      ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> trkpzNeg = trkptNeg.to_double()*sinh(trketaNeg.to_double());
      
      double trkdrpairPhi = sqrt((trkphiPos.to_double() - trkphiNeg.to_double())*(trkphiPos.to_double() - trkphiNeg.to_double()) + (trketaPos.to_double() - trketaNeg.to_double())*(trketaPos.to_double() - trketaNeg.to_double()));
      // write mass calculation here , for hardware specially

      double trkmasspairPhi = sqrt(2*trkptPos.to_double()*trkptNeg.to_double()*(cosh(trketaPos.to_double() - trketaNeg.to_double())-cos(trkphiPos.to_double() - trkphiNeg.to_double())));
      if (trkdrpairPhi > dRmax_) continue; 
      if (trkmasspairPhi < tkpairMmin_ || trkmasspairPhi > tkpairMmax_) continue; // do it before

      ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> trkpxPhi = trkpxNeg + trkpxPos;
      ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> trkpyPhi = trkpyNeg + trkpyPos;
      ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> trkpzPhi = trkpzNeg + trkpzPos;
      
      l1t::TkLightMesonWord::valid_t trkvalidPhi =   trackPosKaon.getValid() && trackNegKaon.getValid();
      l1t::TkLightMesonWord::pt_t trkptPhi = sqrt(pow(trkpxPhi.to_double(),2) + pow(trkpyPhi.to_double(),2)); // use Pow()
      l1t::TkLightMesonWord::glbphi_t trkphiPhi = atan(trkpyPhi.to_double()/trkpxPhi.to_double());
      l1t::TkLightMesonWord::glbeta_t trketaPhi = asinh(trkpzPhi.to_double()/trkptPhi.to_double());
      l1t::TkLightMesonWord::z0_t trkz0Phi = trkz0Pos + trkz0Neg;
      l1t::TkLightMesonWord::mass_t trkmassPhi = sqrt(2*trkptPos.to_double()*trkptNeg.to_double()*(cosh(trketaPos.to_double() - trketaNeg.to_double())-cos(trkphiPos.to_double() - trkphiNeg.to_double())));
      l1t::TkLightMesonWord::type_t trktypePhi = l1t::TkLightMesonWord::TkLightMesonTypes::kPhiType;
      l1t::TkLightMesonWord::ntracks_t trkntracksPhi = 2;
      l1t::TkLightMesonWord::unassigned_t trkunassignedPhi = 0;
      
      l1t::TkLightMesonWord trkPhiWord(trkvalidPhi, trkptPhi, trkphiPhi, trketaPhi, trkz0Phi, trkmassPhi, trktypePhi, trkntracksPhi, trkunassignedPhi);
      
      L1PhiMesonEmulationOutput->push_back(trkPhiWord);
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
  iEvent.put(std::move(L1PhiMesonEmulationOutput), outputCollectionName_);

}

//void L1PhiMesonSelectionEmulationProducer::beginJob() {}

//void L1PhiMesonSelectionEmulationProducer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1PhiMesonSelectionEmulationProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //L1PhiMesonSelectionEmulationProducer
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1PosKaonTracksInputTag", edm::InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"));
  desc.add<edm::InputTag>("l1NegKaonTracksInputTag", edm::InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"));
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
DEFINE_FWK_MODULE(L1PhiMesonSelectionEmulationProducer);
