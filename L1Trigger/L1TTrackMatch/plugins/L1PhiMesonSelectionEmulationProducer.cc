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

//#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuAlgo.h"
//#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuTrackTransform.h"

//#include "hls_math.h"
//
// class declaration

//
using namespace std;
using namespace edm;
using namespace l1t;
//using namespace l1tmetemu;

namespace l1tphimesonemu {
  const unsigned int kInternalPhiWidth{8};
  const unsigned int kGlobalPhiExtra{4};
  static constexpr double minPhi0{-0.7853981696};
  typedef ap_uint<TTTrack_TrackWord::TrackBitWidths::kPhiSize> global_phi_t;
  //typedef ap_uint<15> global_phislice_t;
  const unsigned int kGlobalPhiBins = 1 << kInternalPhiWidth;
  //  const unsigned int kGlobalPhiTotalBins = 1 << TTTrack_TrackWord::TrackBitWidths::kPhiSize;
  //const unsigned int kGlobalPhiTotalBins = 1 << 12;
  const double kStepPhi = (2 * -minPhi0) / kGlobalPhiBins;
  const unsigned int kNSector{9};
  const unsigned int kNQuadrants{4};
  /*  double unpackSignedValue(unsigned int bits, unsigned int nBits){
    int isign = 1;
    unsigned int digitized_maximum = (1 << nBits) - 1;
    if (bits & (1 << (nBits - 1))) {  // check the sign                                                                                                                  
      isign = -1;
      bits = (1 << (nBits + 1)) - bits;  // if negative, flip everything for two's complement encoding                                                                   
    }
    return (double(bits & digitized_maximum)) * isign;
    }*/

  double undigitizeSignedValue(unsigned int twosValue, unsigned int nBits) {
    // Check that none of the bits above the nBits-1 bit, in a range of [0, nBits-1], are set.
    // This makes sure that it isn't possible for the value represented by `twosValue` to be
    //  any bigger than ((1 << nBits) - 1).
    assert((twosValue >> nBits) == 0);

    // Convert from twos compliment to C++ signed integer (normal digitized value)
    int digitizedValue = twosValue;
    if (twosValue & (1 << (nBits - 1))) {  // check if the twosValue is negative
      digitizedValue -= (1 << nBits);
    }

    // Convert to floating point value
    return (double(digitizedValue) + 0.5);
  }

}

class L1PhiMesonSelectionEmulationProducer : public edm::global::EDProducer<> {
public:
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  explicit L1PhiMesonSelectionEmulationProducer(const edm::ParameterSet&);
  ~L1PhiMesonSelectionEmulationProducer() override;

  static constexpr double kmass = 0.493;
  double ETAPHI_LSB = M_PI / (1 << 12);
  double Z0_LSB = 0.05;
  
  l1tphimesonemu::global_phi_t localToGlobalPhi(const TTTrack_TrackWord::phi_t& local_phi, const l1tphimesonemu::global_phi_t& sector_shift) const;

  //  std::vector<l1tphimesonemu::global_phi_t> const getPhiQuad() { return phiQuadrants; }
  //std::vector<l1tphimesonemu::global_phi_t> const getPhiShift() { return phiShift; }
  std::vector<l1tphimesonemu::global_phi_t> generatePhiSliceLUT (unsigned int N){
    float sliceCentre = 0.0;
    std::vector<l1tphimesonemu::global_phi_t> phiLUT;
    for (unsigned int q = 0; q <= N; q++) {
      phiLUT.push_back((l1tphimesonemu::global_phi_t)(sliceCentre / l1tphimesonemu::kStepPhi));
      //    std::cout << "Number of iterations : " << q << "\tslice centre : " << sliceCentre << "\tkStepPhi : " << l1tphimesonemu::kStepPhi << "\tDivision : " << sliceCentre / l1tphimesonemu::kStepPhi << "\tglobal phi_t of division : " << (l1tphimesonemu::global_phi_t)(sliceCentre / l1tphimesonemu::kStepPhi) << "\tphi Lut : " << phiLUT.at(q) << std::endl;

      sliceCentre += 2 * M_PI / N;


    }

    return phiLUT;
  } ;


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

  std::vector<l1tphimesonemu::global_phi_t> phiQuadrants;
  std::vector<l1tphimesonemu::global_phi_t>  phiShift;

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
  //  void generateLUTs(){
  phiQuadrants = L1PhiMesonSelectionEmulationProducer::generatePhiSliceLUT(l1tphimesonemu::kNQuadrants);
  phiShift = L1PhiMesonSelectionEmulationProducer::generatePhiSliceLUT(l1tphimesonemu::kNSector);
  //} ;

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

/*
void L1PhiMesonSelectionEmulationProducer::generateLUTs() { 
  phiQuadrants = L1PhiMesonSelectionEmulationProducer::generatePhiSliceLUT(l1tphimesonemu::kNQuadrants);
  phiShift = L1PhiMesonSelectionEmulationProducer::generatePhiSliceLUT(l1tphimesonemu::kNSector);
  }*/

l1tphimesonemu::global_phi_t L1PhiMesonSelectionEmulationProducer::localToGlobalPhi(const TTTrack_TrackWord::phi_t& local_phi, const l1tphimesonemu::global_phi_t& sector_shift) const {
  int PhiMin = 0;
  int PhiMax = phiQuadrants.back();

  //  std::cout << "PhiMax : " << PhiMax << std::endl;

  int phiMultiplier = TTTrack_TrackWord::TrackBitWidths::kPhiSize - l1tphimesonemu::kInternalPhiWidth;

  int tempPhi = floor(l1tphimesonemu::undigitizeSignedValue(local_phi, TTTrack_TrackWord::TrackBitWidths::kPhiSize) / pow(2, phiMultiplier)) + sector_shift;
  //int tempPhi = floor(l1tphimesonemu::undigitizeSignedValue(local_phi, TTTrack_TrackWord::TrackBitWidths::kPhiSize)) + sector_shift;
  
  // std::cout << "local phi word before conversion : " << local_phi.to_string(2) << "\t local phi value uint :  << " << local_phi.to_uint() << std::endl;

  //  int tempPhi = floor(l1tphimesonemu::undigitizeSignedValue(local_phi.to_uint(), TTTrack_TrackWord::TrackBitWidths::kPhiSize)) + sector_shift;

  //std::cout << "undigitizesigned value : " << l1tphimesonemu::undigitizeSignedValue(local_phi.to_uint(), TTTrack_TrackWord::TrackBitWidths::kPhiSize) << "\tfloor : " << l1tphimesonemu::undigitizeSignedValue(local_phi.to_uint(), TTTrack_TrackWord::TrackBitWidths::kPhiSize) / pow(2, phiMultiplier) << std::endl;

  //  std::cout << "tempPhi floor value : " << tempPhi << "\tsectorshift : " << sector_shift << "\ttemp phi word without phi min : " << l1tphimesonemu::global_phi_t(tempPhi).to_string(2) << std::endl;

  if (tempPhi < PhiMin) {
    tempPhi = tempPhi + PhiMax;
  } else if (tempPhi > PhiMax) {
    tempPhi = tempPhi - PhiMax;
    }  // else                                                                                                                                                             
  //  tempPhi = tempPhi;                                                                                                                                                 
  //  std::cout << "tempPhi final value : " << tempPhi << std::endl;

  l1tphimesonemu::global_phi_t globalPhi = l1tphimesonemu::global_phi_t(tempPhi);
  
  //  std::cout << "global phi : " << globalPhi << std::endl;

  return globalPhi;
}
/*
std::vector<l1tphimesonemu::global_phi_t> L1PhiMesonSelectionEmulationProducer::generatePhiSliceLUT(const unsigned int N) {
  float sliceCentre = 0.0;
  std::vector<l1tphimesonemu::global_phi_t> phiLUT;
  for (unsigned int q = 0; q <= N; q++) {
    phiLUT.push_back((l1tphimesonemu::global_phi_t)(sliceCentre / l1tphimesonemu::kStepPhi));
    sliceCentre += 2 * M_PI / N;
  }
  return phiLUT;
  }*/

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
  
  ap_ufixed<64, 32> etaphi_conv = 1.0 / ETAPHI_LSB;
  ap_ufixed<64, 32> z0_conv = 1.0 * Z0_LSB;

  for (size_t i = 0; i < nPosKaonOutputApproximate; i++) {
    const auto& trackPosKaonRef = l1PosKaonTracksHandle->at(i);
    const auto& trackPosKaon = *trackPosKaonRef;
    
    ap_uint<TrackBitWidths::kPtSize> ptEmulationBitsPos = trackPosKaon.getTrackWord()(TTTrack_TrackWord::TrackBitLocations::kRinvMSB - 1, TTTrack_TrackWord::TrackBitLocations::kRinvLSB);
    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize> ptEmulationPos;
    ptEmulationPos.V = ptEmulationBitsPos.range();
    double trkptPos = ptEmulationPos.to_double();

    TTTrack_TrackWord::tanl_t etaEmulationBitsPos = trackPosKaon.getTanlWord();
    ap_fixed<TrackBitWidths::kEtaSize, TrackBitWidths::kEtaMagSize> etaEmulationPos;
    etaEmulationPos.V = etaEmulationBitsPos.range();
    double trketaPos = etaEmulationPos.to_double();

    //    double trkphiLocalPos = trackPosKaon.getPhi();
    l1tphimesonemu::global_phi_t trkphiEmuPos = L1PhiMesonSelectionEmulationProducer::localToGlobalPhi(trackPosKaon.getPhiWord(), phiShift[trackPosKaon.phiSector()]);
    double trkphiPos = trkphiEmuPos*l1tphimesonemu::kStepPhi ;

    double trkz0Pos = trackPosKaon.getZ0();

    double trkpxPos = trkptPos*cos(trkphiPos);
    double trkpyPos = trkptPos*sin(trkphiPos);
    double trkpzPos = trkptPos*sinh(trketaPos);

    //    std::cout << "nNegKaonOutputApproximate : " << nNegKaonOutputApproximate << std::endl;

    //    std::cout << "poskaon pt : " << trkptPos << "\tlocalphi : " << trackPosKaon.getPhi() << "\tlocalphi word : " << trackPosKaon.getPhiWord().to_string(2) << "\t globalphi poskaon : " << trkphiEmuPos.to_uint() << "\t globalphi poskaon word : " << trkphiEmuPos.to_string(2) << "\tAlexx global phi : " << trkphiPos << std::endl;
    
    for (size_t j = 0; j < nNegKaonOutputApproximate; j++) {
      const auto& trackNegKaonRef = l1NegKaonTracksHandle->at(j);
      const auto& trackNegKaon = *trackNegKaonRef;
      
      ap_uint<TrackBitWidths::kPtSize> ptEmulationBitsNeg = trackNegKaon.getTrackWord()(TTTrack_TrackWord::TrackBitLocations::kRinvMSB - 1, TTTrack_TrackWord::TrackBitLocations::kRinvLSB);
    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize> ptEmulationNeg;
    ptEmulationNeg.V = ptEmulationBitsNeg.range();
    double trkptNeg = ptEmulationNeg.to_double();

    TTTrack_TrackWord::tanl_t etaEmulationBitsNeg = trackNegKaon.getTanlWord();
    ap_fixed<TrackBitWidths::kEtaSize, TrackBitWidths::kEtaMagSize> etaEmulationNeg;
    etaEmulationNeg.V = etaEmulationBitsNeg.range();
    double trketaNeg = etaEmulationNeg.to_double();

    //double trkphiLocalNeg = trackNegKaon.getPhi();
    l1tphimesonemu::global_phi_t trkphiEmuNeg = L1PhiMesonSelectionEmulationProducer::localToGlobalPhi(trackNegKaon.getPhiWord(), phiShift[trackNegKaon.phiSector()]);
    double trkphiNeg = trkphiEmuNeg*l1tphimesonemu::kStepPhi;

    //std::cout << "global pt neg : " << trkptNeg << "global phi neg : " << trkphiNeg << std::endl;

    double trkz0Neg = trackNegKaon.getZ0();

    double trkpxNeg = trkptNeg*cos(trkphiNeg);
    double trkpyNeg = trkptNeg*sin(trkphiNeg);
    double trkpzNeg = trkptNeg*sinh(trketaNeg);

    double convdPhi = trkphiPos - trkphiNeg;

    if (convdPhi < 0 ) convdPhi = convdPhi + 2*M_PI;
    else if (convdPhi > 2*M_PI) convdPhi = convdPhi - 2*M_PI;
      
    double trkdrpairPhi = sqrt(pow(convdPhi,2) + pow((trketaPos - trketaNeg),2));
      // write mass calculation here , for hardware specially

    /*    std::cout << "Phi Emu trkptPos : " << trkptPos << "\t trkptNeg : " << trkptNeg << "\t trketaPos : " << trketaPos << "\t trketaNeg : " << trketaNeg << "\t trkPhiPos : " << trkphiPos << "\t trkPhiNeg : " << trkphiNeg << std::endl;
    std::cout << "Phi Emu Pos Kaon track cos(trkphiPos) : " << cos(trkphiPos) << "\tsin(trkphiPos) : " << sin(trkphiPos) << "\tsinh(trketaPos) : " << sinh(trketaPos) << "\tcosh(trketaPos)" << std::endl;
    std::cout << "Phi Emu Neg Kaon track cos(trkphiNeg) : " << cos(trkphiNeg) << "\tsin(trkphiNeg) : " << sin(trkphiNeg) << "\tsinh(trketaNeg) : " << sinh(trketaNeg) << "\tcosh(trketaPos)" << std::endl;

    std::cout << "Phi Emu Pos Kaon track px : " << trkpxPos << "\tpy : " << trkpyPos << "\tpz : " << trkpzPos << std::endl;
    std::cout << "Phi Emu Neg Kaon track px : " << trkpxNeg << "\tpy : " << trkpyNeg << "\tpz : " << trkpzNeg << std::endl;

    std::cout << "Phi Emu cosh in mass : " << cosh(trketaPos - trketaNeg) << "\t cos in mass : " << cos(trkphiPos - trkphiNeg) << std::endl;
    std::cout << "Phi Emu diff in dR for eta : " << (trketaPos - trketaNeg) << "\t for phi : " << (trkphiPos - trkphiNeg) << std::endl;
    std::cout << "Phi Emu diff square in dR for eta : " << pow((trketaPos - trketaNeg),2) << "\t for phi : " << pow((trkphiPos - trkphiNeg),2) << std::endl;*/

    double trkmasspairPhi = sqrt(2*trkptPos*trkptNeg*(cosh(trketaPos - trketaNeg)-cos(trkphiPos - trkphiNeg)));
    
    if (i == 0 && j == 2) {
      //std::cout << "Emu pos tk pt : " << trkptPos << "\t postk eta : " << trketaPos << "\t postk phi : " << trkphiPos << std::endl;
      //std::cout << "Emu neg tk pt : " << trkptNeg << "\t negtk eta : " << trketaNeg << "\t negtk phi : " << trkphiNeg << std::endl;
      std::cout << "Phi Emu trkptPos : " << trkptPos << "\t trkptNeg : " << trkptNeg << "\t trketaPos : " << trketaPos << "\t trketaNeg : " << trketaNeg << "\t trkPhiPos : " << trkphiPos << "\t trkPhiNeg : " << trkphiNeg << std::endl;
      std::cout << "Phi Emu Pos Kaon track cos(trkphiPos) : " << cos(trkphiPos) << "\tsin(trkphiPos) : " << sin(trkphiPos) << "\tsinh(trketaPos) : " << sinh(trketaPos) << "\tcosh(trketaPos)" << std::endl;
      std::cout << "Phi Emu Neg Kaon track cos(trkphiNeg) : " << cos(trkphiNeg) << "\tsin(trkphiNeg) : " << sin(trkphiNeg) << "\tsinh(trketaNeg) : " << sinh(trketaNeg) << "\tcosh(trketaPos)" << std::endl;
      
      std::cout << "Phi Emu Pos Kaon track px : " << trkpxPos << "\tpy : " << trkpyPos << "\tpz : " << trkpzPos << std::endl;
      std::cout << "Phi Emu Neg Kaon track px : " << trkpxNeg << "\tpy : " << trkpyNeg << "\tpz : " << trkpzNeg << std::endl;
      
      std::cout << "Phi Emu cosh in mass : " << cosh(trketaPos - trketaNeg) << "\t cos in mass : " << cos(trkphiPos - trkphiNeg) << std::endl;
      std::cout << "Phi Emu diff in dR for eta : " << (trketaPos - trketaNeg) << "\t for phi : " << (trkphiPos - trkphiNeg) << std::endl;
      std::cout << "Phi Emu diff square in dR for eta : " << pow((trketaPos - trketaNeg),2) << "\t for phi : " << pow((trkphiPos - trkphiNeg),2) << std::endl;
      std::cout << "Emu Masspair phi : " << trkmasspairPhi << "\ttrk dr pair : " << trkdrpairPhi << std::endl;
    }


    //std::cout << "trkmass pair phi beforfe mass cut : " << trkmasspairPhi << std::endl;
    
    //std::cout << "trkdr pair phi beforfe mass cut : " << trkdrpairPhi << std::endl;
    if (trkdrpairPhi > dRmax_) continue; 
    if (trkmasspairPhi < tkpairMmin_ || trkmasspairPhi > tkpairMmax_) continue; // do it before
    
    double trkpxPhi = trkpxNeg + trkpxPos;
    double trkpyPhi = trkpyNeg + trkpyPos;
    double trkpzPhi = trkpzNeg + trkpzPos;


    /* std::cout << "phi emul pt in double format : " << sqrt(pow(trkpxPhi,2) + pow(trkpyPhi,2)) << std::endl;
    std::cout << "phi emul eta in double format : " << asinh(trkpzPhi/sqrt(pow(trkpxPhi,2) + pow(trkpyPhi,2))) << std::endl;

    std::cout << "trk phi emulation mass inside analyzer (original and emulation ) double format: " << sqrt(2*trkptPos*trkptNeg*(cosh(trketaPos - trketaNeg)-cos(trkphiPos - trkphiNeg)))  <<  std::endl;*/

      l1t::TkLightMesonWord::valid_t trkvalidPhi =   trackPosKaon.getValid() && trackNegKaon.getValid();
      l1t::TkLightMesonWord::pt_t trkptPhi = sqrt(pow(trkpxPhi,2) + pow(trkpyPhi,2)); // use Pow()
      l1t::TkLightMesonWord::glbphi_t trkphiPhi = atan2(trkpyPhi,trkpxPhi) / ETAPHI_LSB;
      //l1t::TkLightMesonWord::glbphi_t trkphiPhi = asin(trkpyPhi/sqrt(pow(trkpxPhi,2) + pow(trkpyPhi,2))) / ETAPHI_LSB;
      l1t::TkLightMesonWord::glbeta_t trketaPhi = asinh(trkpzPhi/sqrt(pow(trkpxPhi,2) + pow(trkpyPhi,2))) / ETAPHI_LSB;
      l1t::TkLightMesonWord::z0_t trkz0Phi = ((trkz0Pos + trkz0Neg) / Z0_LSB)* 0.5;
      l1t::TkLightMesonWord::mass_t trkmassPhi = sqrt(2*trkptPos*trkptNeg*(cosh(trketaPos - trketaNeg)-cos(trkphiPos - trkphiNeg)));
      l1t::TkLightMesonWord::type_t trktypePhi = l1t::TkLightMesonWord::TkLightMesonTypes::kPhiType;
      l1t::TkLightMesonWord::ntracks_t trkntracksPhi = 2;
      l1t::TkLightMesonWord::unassigned_t trkunassignedPhi = 0;

      /*      std::cout << "phi emul pt just before booking trkphiWord (original and emulation ): : " << trkptPhi.to_double() << std::endl;
      std::cout << "phi emul eta just before booking trkphiWord (original and emulation ): : " << trketaPhi.to_double() << std::endl;
      std::cout << "phi emul phi just before booking trkphiWord (original and emulation ): : " << trkphiPhi.to_double() << std::endl;
      
      std::cout << "trkPhiword before" << std::endl;

      //std::cout << __PRETTY_FUNCTION__ << __LINE__ << std::endl;

      std::cout << "trk phi emulation mass inside analyzer (original and emulation ) before phi booking: " << trkmassPhi  <<  std::endl;*/

      l1t::TkLightMesonWord trkPhiWord(trkvalidPhi, trkptPhi, trkphiPhi, trketaPhi, trkz0Phi, trkmassPhi, trktypePhi, trkntracksPhi, trkunassignedPhi);
      
      /*      std::cout << "trkPhiword after" << std::endl;
      std::cout << "trk phi emulation pt inside analyzer (original and emulation ): " << trkPhiWord.pt()  <<  std::endl;
      std::cout << "trk phi emulation eta inside analyzer (original and emulation ): " << trkPhiWord.glbeta()  <<  std::endl;
      std::cout << "trk phi emulation phi inside analyzer (original and emulation ): " << trkPhiWord.glbphi()  <<  std::endl;
      std::cout << "trk phi emulation mass inside analyzer (original and emulation ): " << trkPhiWord.mass()  <<  std::endl;*/

      L1PhiMesonEmulationOutput->push_back(trkPhiWord);

      
    }
  }

  std::cout << "# Phi meson emu : " << L1PhiMesonEmulationOutput->size() << std::endl;

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
