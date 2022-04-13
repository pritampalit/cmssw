#ifndef DQM_TrackingMonitorSource_RecoTrackMonitor_h
#define DQM_TrackingMonitorSource_RecoTrackMonitor_h

#include <string>
#include <vector>
#include <map>
#include <set>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

class BeamSpot;
class TrackCollection;
class VertexCollection;
class TrackingRecHit;

class RecoTrackMonitor : public DQMEDAnalyzer {
 public:
  RecoTrackMonitor( const edm::ParameterSet& );

 protected:

  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& eSetup);
  void processHit(const TrackingRecHit& recHit, edm::EventSetup const& iSetup, const TrackerGeometry& tkGeom, double wfac=1);
  void processClusters(const edm::Handle<edmNew::DetSetVector<SiStripCluster>>& clusterHandle, edm::EventSetup const& iSetup, const TrackerGeometry& tkGeom, double wfac=1);
  void addClusterToMap(uint32_t detid, const SiStripCluster* cluster);
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &);

 private:

  edm::ParameterSet parameters_;

  std::string moduleName_;
  std::string folderName_;
  const edm::InputTag trackTag_;
  const edm::InputTag bsTag_;
  const edm::InputTag vertexTag_;
  const edm::InputTag puSummaryTag_;
  const edm::InputTag clusterTag_;
  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puSummaryToken_; 
  const edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster> > clusterToken_;

  const std::string trackQuality_;
  const bool doPUCorrection_;
  const bool isMC_;
  const bool haveAllHistograms_;
  const std::string puScaleFactorFile_;
  const bool verbose_;

  MonitorElement* residualXPBH_;
  MonitorElement* residualXPEH_;
  MonitorElement* residualXTIBH_;
  MonitorElement* residualXTOBH_;
  MonitorElement* residualXTIDH_;
  MonitorElement* residualXTECH_;
  MonitorElement* residualYPBH_;
  MonitorElement* residualYPEH_;
  MonitorElement* residualYTIBH_;
  MonitorElement* residualYTOBH_;
  MonitorElement* residualYTIDH_;
  MonitorElement* residualYTECH_;
 
  MonitorElement* hitmap_pixelH_;
  MonitorElement* hitEta_pixelH_;
  MonitorElement* hitPhi_pixelH_;
  MonitorElement* hitmap2_pixelH_;
  MonitorElement* hitmap_stripH_;
  MonitorElement* hitEta_stripH_;
  MonitorElement* hitPhi_stripH_;
  MonitorElement* hitmap2_stripH_;

  MonitorElement* hitmap_PixB_LayerH_[4];
  MonitorElement* hitEta_PixB_LayerH_[4];
  MonitorElement* hitPhi_PixB_LayerH_[4];
  MonitorElement* hitmap2_PixB_LayerH_[4];

  MonitorElement* hitmap_TIB_LayerH_[4];
  MonitorElement* hitEta_TIB_LayerH_[4];
  MonitorElement* hitPhi_TIB_LayerH_[4];
  MonitorElement* hitmap2_TIB_LayerH_[4];

  MonitorElement* hitmap_TOB_LayerH_[6]; 
  MonitorElement* hitEta_TOB_LayerH_[6];
  MonitorElement* hitPhi_TOB_LayerH_[6];
  MonitorElement* hitmap2_TOB_LayerH_[6];

  // Exclusive Quantities

  MonitorElement* nHitsTIBSVsEtaH_;
  MonitorElement* nHitsTOBSVsEtaH_;
  MonitorElement* nHitsTECSVsEtaH_;
  MonitorElement* nHitsTIDSVsEtaH_;
  MonitorElement* nHitsStripSVsEtaH_;

  MonitorElement* nHitsTIBDVsEtaH_;
  MonitorElement* nHitsTOBDVsEtaH_;
  MonitorElement* nHitsTECDVsEtaH_;
  MonitorElement* nHitsTIDDVsEtaH_;
  MonitorElement* nHitsStripDVsEtaH_;

  MonitorElement* hOnTrkClusChargeThinH_;
  MonitorElement* hOnTrkClusWidthThinH_;
  MonitorElement* hOnTrkClusChargeThickH_;
  MonitorElement* hOnTrkClusWidthThickH_;

  MonitorElement* hOffTrkClusChargeThinH_;
  MonitorElement* hOffTrkClusWidthThinH_;
  MonitorElement* hOffTrkClusChargeThickH_;
  MonitorElement* hOffTrkClusWidthThickH_;

  unsigned long long m_cacheID_;

  std::vector<float> vpu_;
  std::map<uint32_t, std::set<const SiStripCluster*> > clusterMap_;
};
#endif
