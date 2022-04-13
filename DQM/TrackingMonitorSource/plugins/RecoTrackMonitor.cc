#include "TFile.h"
#include "TH1.h"
#include "TMath.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TPRegexp.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "Geometry/CommonDetUnit/interface/GluedGeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "RecoLocalTracker/SiStripClusterizer/interface/SiStripClusterInfo.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DQM/TrackingMonitorSource/interface/RecoTrackMonitor.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/TrackerGeomDet.h"

// -----------------------------
// constructors and destructor
// -----------------------------
RecoTrackMonitor::RecoTrackMonitor(const edm::ParameterSet& ps): 
  parameters_(ps),
  moduleName_(parameters_.getUntrackedParameter<std::string>("moduleName", "StandaloneTrackMonitor")),
  folderName_(parameters_.getUntrackedParameter<std::string>("folderName", "highPurityTracks")),
  trackTag_(parameters_.getUntrackedParameter<edm::InputTag>("trackInputTag", edm::InputTag("generalTracks"))),
  bsTag_(parameters_.getUntrackedParameter<edm::InputTag>("offlineBeamSpot", edm::InputTag("offlineBeamSpot"))),
  vertexTag_(parameters_.getUntrackedParameter<edm::InputTag>("vertexTag", edm::InputTag("offlinePrimaryVertices"))),
  //primaryVertexInputTag_(parameters_.getUntrackedParameter<edm::InputTag>("primaryVertexInputTag", edm::InputTag("primaryVertex"))),
  puSummaryTag_(parameters_.getUntrackedParameter<edm::InputTag>("puTag", edm::InputTag("addPileupInfo"))),
  clusterTag_(parameters_.getUntrackedParameter<edm::InputTag>("clusterTag", edm::InputTag("siStripClusters"))),
  trackToken_(consumes<reco::TrackCollection>(trackTag_)),
  bsToken_(consumes<reco::BeamSpot>(bsTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  //pvToken_(consumes<reco::VertexCollection>(primaryVertexInputTag_)),
  puSummaryToken_(consumes<std::vector<PileupSummaryInfo> >(puSummaryTag_)),
  clusterToken_(consumes<edmNew::DetSetVector<SiStripCluster> >(clusterTag_)),
  trackQuality_(parameters_.getUntrackedParameter<std::string>("trackQuality", "highPurity")),
  doPUCorrection_(parameters_.getUntrackedParameter<bool>("doPUCorrection", false)),
  isMC_(parameters_.getUntrackedParameter<bool>("isMC", false)),
  haveAllHistograms_(parameters_.getUntrackedParameter<bool>("haveAllHistograms", false)),
  puScaleFactorFile_(parameters_.getUntrackedParameter<std::string>("puScaleFactorFile", "PileupScaleFactor.root")),
  verbose_(parameters_.getUntrackedParameter<bool>("verbose", false))
{
  nHitsTIBSVsEtaH_ = nullptr;
  nHitsTOBSVsEtaH_ = nullptr;
  nHitsTECSVsEtaH_ = nullptr;
  nHitsTIDSVsEtaH_ = nullptr;
  nHitsStripSVsEtaH_ = nullptr;

  nHitsTIBDVsEtaH_ = nullptr;
  nHitsTOBDVsEtaH_ = nullptr;
  nHitsTECDVsEtaH_ = nullptr;
  nHitsTIDDVsEtaH_ = nullptr;
  nHitsStripDVsEtaH_ = nullptr;

  hOnTrkClusChargeThinH_ = nullptr;
  hOnTrkClusWidthThinH_ = nullptr;
  hOnTrkClusChargeThickH_ = nullptr;
  hOnTrkClusWidthThickH_ = nullptr;
  
  hOffTrkClusChargeThinH_ = nullptr;
  hOffTrkClusWidthThinH_ = nullptr;
  hOffTrkClusChargeThickH_ = nullptr;
  hOffTrkClusWidthThickH_ = nullptr;
  hitmap_pixelH_ = nullptr;
  hitmap_stripH_ = nullptr;

  // Read pileup weight factors
  if (isMC_ && doPUCorrection_) {
    vpu_.clear();
    TFile* f1 = TFile::Open(puScaleFactorFile_.c_str());
    TH1F* h1 = dynamic_cast<TH1F*>(f1->Get("pileupweight"));
    for (int i = 1; i <= h1->GetNbinsX(); ++i) vpu_.push_back(h1->GetBinContent(i));
    f1->Close();
  }

}

void RecoTrackMonitor::bookHistograms(DQMStore::IBooker &iBook, edm::Run const& iRun, edm::EventSetup const& iSetup) {
  edm::ParameterSet TrackEtaHistoPar = parameters_.getParameter<edm::ParameterSet>("trackEtaH");
  edm::ParameterSet TrackPtHistoPar = parameters_.getParameter<edm::ParameterSet>("trackPtH");

  std::string currentFolder = moduleName_ + "/" + folderName_ ;
  iBook.setCurrentFolder(currentFolder.c_str());

  // The following are common with the official tool
  if (haveAllHistograms_) {
    residualXPBH_ = iBook.book1D("residualXPixelBarrel", "Residual in X in Pixel Barrel", 20, -5.0, 5.0);
    residualXPEH_ = iBook.book1D("residualXPixelEndcap", "Residual in X in Pixel Endcap", 20, -5.0, 5.0);
    residualXTIBH_ = iBook.book1D("residualXStripTIB", "Residual in X in Strip TIB", 20, -5.0, 5.0);
    residualXTOBH_ = iBook.book1D("residualXStripTOB", "Residual in X in Strip TOB", 20, -5.0, 5.0);
    residualXTECH_ = iBook.book1D("residualXStripTEC", "Residual in X in Strip TEC", 20, -5.0, 5.0);
    residualXTIDH_ = iBook.book1D("residualXStripTID", "Residual in X in Strip TID", 20, -5.0, 5.0);
    residualYPBH_ = iBook.book1D("residualYPixelBarrel", "Residual in Y in Pixel Barrel", 20, -5.0, 5.0);
    residualYPEH_ = iBook.book1D("residualYPixelEndcap", "Residual in Y in Pixel Endcap", 20, -5.0, 5.0);
    residualYTIBH_ = iBook.book1D("residualYStripTIB", "Residual in Y in Strip TIB", 20, -5.0, 5.0);
    residualYTOBH_ = iBook.book1D("residualYStripTOB", "Residual in Y in Strip TOB", 20, -5.0, 5.0);
    residualYTECH_ = iBook.book1D("residualYStripTEC", "Residual in Y in Strip TEC", 20, -5.0, 5.0);
    residualYTIDH_ = iBook.book1D("residualYStripTID", "Residual in Y in Strip TID", 20, -5.0, 5.0);
  }
  // Exclusive histograms
  nHitsTIBSVsEtaH_ = iBook.bookProfile("nHitsTIBSVsEta", "Number of Hits in TIB Vs Eta (Single-sided)",
                                       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
                                       TrackEtaHistoPar.getParameter<double>("Xmin"),
                                       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsTOBSVsEtaH_ = iBook.bookProfile("nHitsTOBSVsEta", "Number of Hits in TOB Vs Eta (Single-sided)",
                                       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
                                       TrackEtaHistoPar.getParameter<double>("Xmin"),
                                       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsTECSVsEtaH_ = iBook.bookProfile("nHitsTECSVsEta", "Number of Hits in TEC Vs Eta (Single-sided)",
                                       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
                                       TrackEtaHistoPar.getParameter<double>("Xmin"),
                                       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsTIDSVsEtaH_ = iBook.bookProfile("nHitsTIDSVsEta", "Number of Hits in TID Vs Eta (Single-sided)",
                                       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
                                       TrackEtaHistoPar.getParameter<double>("Xmin"),
                                       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");

  nHitsStripSVsEtaH_ = iBook.bookProfile("nHitsStripSVsEta", "Number of Strip Hits Vs Eta (Single-sided)",
                                         TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
                                         TrackEtaHistoPar.getParameter<double>("Xmin"),
                                         TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");

  nHitsTIBDVsEtaH_ = iBook.bookProfile("nHitsTIBDVsEta", "Number of Hits in TIB Vs Eta (Double-sided)",
                                       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
                                       TrackEtaHistoPar.getParameter<double>("Xmin"),
                                       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsTOBDVsEtaH_ = iBook.bookProfile("nHitsTOBDVsEta", "Number of Hits in TOB Vs Eta (Double-sided)",
                                       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
                                       TrackEtaHistoPar.getParameter<double>("Xmin"),
                                       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsTECDVsEtaH_ = iBook.bookProfile("nHitsTECDVsEta", "Number of Hits in TEC Vs Eta (Double-sided)",
                                       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
                                       TrackEtaHistoPar.getParameter<double>("Xmin"),
                                       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsTIDDVsEtaH_ = iBook.bookProfile("nHitsTIDDVsEta", "Number of Hits in TID Vs Eta (Double-sided)",
                                       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
                                       TrackEtaHistoPar.getParameter<double>("Xmin"),
                                       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsStripDVsEtaH_ = iBook.bookProfile("nHitsStripDVsEta", "Number of Strip Hits Vs Eta (Double-sided)",
                                         TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
                                         TrackEtaHistoPar.getParameter<double>("Xmin"),
                                         TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");

  // On and off-track cluster properties                                                                                                                
  hOnTrkClusChargeThinH_ = iBook.book1D("hOnTrkClusChargeThin", "On-track Cluster Charge (Thin Sensor)", 100, 0, 1000);
  hOnTrkClusWidthThinH_ = iBook.book1D("hOnTrkClusWidthThin", "On-track Cluster Width (Thin Sensor)", 20, -0.5, 19.5);
  hOnTrkClusChargeThickH_ = iBook.book1D("hOnTrkClusChargeThick", "On-track Cluster Charge (Thick Sensor)", 100, 0, 1000);
  hOnTrkClusWidthThickH_ = iBook.book1D("hOnTrkClusWidthThick", "On-track Cluster Width (Thick Sensor)", 20, -0.5, 19.5);

  hOffTrkClusChargeThinH_ = iBook.book1D("hOffTrkClusChargeThin", "Off-track Cluster Charge (Thin Sensor)", 100, 0, 1000);
  hOffTrkClusWidthThinH_ = iBook.book1D("hOffTrkClusWidthThin", "Off-track Cluster Width (Thin Sensor)", 20, -0.5, 19.5);
  hOffTrkClusChargeThickH_ = iBook.book1D("hOffTrkClusChargeThick", "Off-track Cluster Charge (Thick Sensor)", 100, 0, 1000);
  hOffTrkClusWidthThickH_ = iBook.book1D("hOffTrkClusWidthThick", "Off-track Cluster Width (Thick Sensor)", 20, -0.5, 19.5);

  // RecHits
  hitmap_pixelH_ = iBook.book2D("hitmap_pixel","2D hitefficiency map in Global eta-phi plane of Pixel", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_pixelH_ = iBook.book1D("hitEta_pixel","Eta distribution in rechit Pixel", 60, -3.0, 3.0);
  hitPhi_pixelH_ = iBook.book1D("hitPhi_pixel","Phi distribution in rechit Pixel", 70, -3.5, 3.5);
  hitmap2_pixelH_ = iBook.book2D("hitmap2_pixel","2D hitefficiency map in Global z-phi plane of Pixel", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  hitmap_stripH_ = iBook.book2D("hitmap_strip","2D hitefficiency map in Global eta-phi plane of Strip", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_stripH_ = iBook.book1D("hitEta_strip","Eta distribution in rechit Strip", 60, -3.0, 3.0);
  hitPhi_stripH_ = iBook.book1D("hitPhi_strip","Phi distribution in rechit Strip", 70, -3.5, 3.5);  
  hitmap2_stripH_ = iBook.book2D("hitmap2_strip","2D hitefficiency map in Global z-phi plane of STrip", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  // Pixel Barrel
  hitmap_PixB_LayerH_[0] = iBook.book2D("hitmap_PixB_L1","2D hitefficiency map in Global eta-phi PixBL1", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_PixB_LayerH_[0] = iBook.book1D("hitEta_PixB_L1","Eta distribution in rechit Pixel barrel L1", 60, -3.0, 3.0);
  hitPhi_PixB_LayerH_[0] = iBook.book1D("hitPhi_PixB_L1","Phi distribution in rechit Pixel barrel L1", 70, -3.5, 3.5);
  hitmap2_PixB_LayerH_[0] = iBook.book2D("hitmap2_PixB_L1","2D hitefficiency map in Global z-phi plane of PixBL1", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  hitmap_PixB_LayerH_[1] = iBook.book2D("hitmap_PixB_L2","2D hitefficiency map in Global eta-phi PixBL2", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_PixB_LayerH_[1] = iBook.book1D("hitEta_PixB_L2","Eta distribution in rechit Pixel barrel L2", 60, -3.0, 3.0);
  hitPhi_PixB_LayerH_[1] = iBook.book1D("hitPhi_PixB_L2","Phi distribution in rechit Pixel barrel L2", 70, -3.5, 3.5);
  hitmap2_PixB_LayerH_[1] = iBook.book2D("hitmap2_PixB_L2","2D hitefficiency map in Global z-phi plane of PixBL2", 120, -20.0, 20.0 , 200, -3.2, 3.2); 

  hitmap_PixB_LayerH_[2] = iBook.book2D("hitmap_PixB_L3","2D hitefficiency map in Global eta-phi PixBL3", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_PixB_LayerH_[2] = iBook.book1D("hitEta_PixB_L3","Eta distribution in rechit Pixel barrel L3", 60, -3.0, 3.0);
  hitPhi_PixB_LayerH_[2] = iBook.book1D("hitPhi_PixB_L3","Phi distribution in rechit Pixel barrel L3", 70, -3.5, 3.5);
  hitmap2_PixB_LayerH_[2] = iBook.book2D("hitmap2_PixB_L3","2D hitefficiency map in Global z-phi plane of PixBL3", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  hitmap_PixB_LayerH_[3] = iBook.book2D("hitmap_PixB_L4","2D hitefficiency map in Global eta-phi PixBL4", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_PixB_LayerH_[3] = iBook.book1D("hitEta_PixB_L4","Eta distribution in rechit Pixel barrel L4", 60, -3.0, 3.0);
  hitPhi_PixB_LayerH_[3] = iBook.book1D("hitPhi_PixB_L4","Phi distribution in rechit Pixel barrel L4", 70, -3.5, 3.5);
  hitmap2_PixB_LayerH_[3] = iBook.book2D("hitmap2_PixB_L4","2D hitefficiency map in Global z-phi plane of PixBL4", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  // TIB
  hitmap_TIB_LayerH_[0] = iBook.book2D("hitmap_TIB_L1","2D hitefficiency map in Global eta-phi TIBL1", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_TIB_LayerH_[0] = iBook.book1D("hitEta_TIB_L1","Eta distribution in rechit TIBL1", 60, -3.0, 3.0);
  hitPhi_TIB_LayerH_[0] = iBook.book1D("hitPhi_TIB_L1","Phi distribution in rechit TIBL1", 70, -3.5, 3.5);
  hitmap2_TIB_LayerH_[0] = iBook.book2D("hitmap2_TIB_L1","2D hitefficiency map in Global z-phi plane of TIBL1", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  hitmap_TIB_LayerH_[1] = iBook.book2D("hitmap_TIB_L2","2D hitefficiency map in Global eta-phi TIBL2", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_TIB_LayerH_[1] = iBook.book1D("hitEta_TIB_L2","Eta distribution in rechit TIBL2", 60, -3.0, 3.0);
  hitPhi_TIB_LayerH_[1] = iBook.book1D("hitPhi_TIB_L2","Phi distribution in rechit TIBL2", 70, -3.5, 3.5);
  hitmap2_TIB_LayerH_[1] = iBook.book2D("hitmap2_TIB_L2","2D hitefficiency map in Global z-phi plane of TIBL2", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  hitmap_TIB_LayerH_[2] = iBook.book2D("hitmap_TIB_L3","2D hitefficiency map in Global eta-phi TIBL3", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_TIB_LayerH_[2] = iBook.book1D("hitEta_TIB_L3","Eta distribution in rechit TIBL3", 60, -3.0, 3.0);
  hitPhi_TIB_LayerH_[2] = iBook.book1D("hitPhi_TIB_L3","Phi distribution in rechit TIBL3", 70, -3.5, 3.5);
  hitmap2_TIB_LayerH_[2] = iBook.book2D("hitmap2_TIB_L3","2D hitefficiency map in Global z-phi plane of TIBL3", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  hitmap_TIB_LayerH_[3] = iBook.book2D("hitmap_TIB_L4","2D hitefficiency map in Global eta-phi TIBL4", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_TIB_LayerH_[3] = iBook.book1D("hitEta_TIB_L4","Eta distribution in rechit TIBL4", 60, -3.0, 3.0);
  hitPhi_TIB_LayerH_[3] = iBook.book1D("hitPhi_TIB_L4","Phi distribution in rechit TIBL4", 70, -3.5, 3.5);
  hitmap2_TIB_LayerH_[3] = iBook.book2D("hitmap2_TIB_L4","2D hitefficiency map in Global z-phi plane of TIBL4", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  // TOB
  hitmap_TOB_LayerH_[0] = iBook.book2D("hitmap_TOB_L1","2D hitefficiency map in Global eta-phi TOBL1", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_TOB_LayerH_[0] = iBook.book1D("hitEta_TOB_L1","Eta distribution in rechit TOBL1", 60, -3.0, 3.0);
  hitPhi_TOB_LayerH_[0] = iBook.book1D("hitPhi_TOB_L1","Phi distribution in rechit TOBL1", 70, -3.5, 3.5);
  hitmap2_TOB_LayerH_[0] = iBook.book2D("hitmap2_TOB_L1","2D hitefficiency map in Global z-phi plane of TOBL1", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  hitmap_TOB_LayerH_[1] = iBook.book2D("hitmap_TOB_L2","2D hitefficiency map in Global eta-phi TOBL2", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_TOB_LayerH_[1] = iBook.book1D("hitEta_TOB_L2","Eta distribution in rechit TOBL2", 60, -3.0, 3.0);
  hitPhi_TOB_LayerH_[1] = iBook.book1D("hitPhi_TOB_L2","Phi distribution in rechit TOBL2", 70, -3.5, 3.5);
  hitmap2_TOB_LayerH_[1] = iBook.book2D("hitmap2_TOB_L2","2D hitefficiency map in Global z-phi plane of TOBL2", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  hitmap_TOB_LayerH_[2] = iBook.book2D("hitmap_TOB_L3","2D hitefficiency map in Global eta-phi TOBL3", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_TOB_LayerH_[2] = iBook.book1D("hitEta_TOB_L3","Eta distribution in rechit TOBL3", 60, -3.0, 3.0);
  hitPhi_TOB_LayerH_[2] = iBook.book1D("hitPhi_TOB_L3","Phi distribution in rechit TOBL3", 70, -3.5, 3.5);
  hitmap2_TOB_LayerH_[2] = iBook.book2D("hitmap2_TOB_L3","2D hitefficiency map in Global z-phi plane of TOBL3", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  hitmap_TOB_LayerH_[3] = iBook.book2D("hitmap_TOB_L4","2D hitefficiency map in Global eta-phi TOBL4", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_TOB_LayerH_[3] = iBook.book1D("hitEta_TOB_L4","Eta distribution in rechit TOBL4", 60, -3.0, 3.0);
  hitPhi_TOB_LayerH_[3] = iBook.book1D("hitPhi_TOB_L4","Phi distribution in rechit TOBL4", 70, -3.5, 3.5);
  hitmap2_TOB_LayerH_[3] = iBook.book2D("hitmap2_TOB_L4","2D hitefficiency map in Global z-phi plane of TOBL4", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  hitmap_TOB_LayerH_[4] = iBook.book2D("hitmap_TOB_L5","2D hitefficiency map in Global eta-phi TOBL5", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_TOB_LayerH_[4] = iBook.book1D("hitEta_TOB_L5","Eta distribution in rechit TOBL5", 60, -3.0, 3.0);
  hitPhi_TOB_LayerH_[4] = iBook.book1D("hitPhi_TOB_L5","Phi distribution in rechit TOBL5", 70, -3.5, 3.5);
  hitmap2_TOB_LayerH_[4] = iBook.book2D("hitmap2_TOB_L5","2D hitefficiency map in Global z-phi plane of TOBL5", 120, -20.0, 20.0 , 200, -3.2, 3.2);

  hitmap_TOB_LayerH_[5] = iBook.book2D("hitmap_TOB_L6","2D hitefficiency map in Global eta-phi TOBL6", 60, -3.0, 3.0 , 70, -3.5, 3.5);
  hitEta_TOB_LayerH_[5] = iBook.book1D("hitEta_TOB_L6","Eta distribution in rechit TOBL6", 60, -3.0, 3.0);
  hitPhi_TOB_LayerH_[5] = iBook.book1D("hitPhi_TOB_L6","Phi distribution in rechit TOBL6", 70, -3.5, 3.5);
  hitmap2_TOB_LayerH_[5] = iBook.book2D("hitmap2_TOB_L6","2D hitefficiency map in Global z-phi plane of TOBL6", 120, -20.0, 20.0 , 200, -3.2, 3.2);  
}

void RecoTrackMonitor::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  if (verbose_) std::cout << "Begin RecoTrackMonitor" << std::endl;

  // Get event setup (to get global transformation)                                                                                                     
  edm::ESHandle<TrackerGeometry> geomHandle;                                                                                                          
  iSetup.get<TrackerDigiGeometryRecord>().get(geomHandle);                                                                                            
  const TrackerGeometry& tkGeom = (*geomHandle);                                                                                                      

  edm::ESHandle<TrackerTopology> tTopoHandle_;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle_);

  std::string geomType_;
  edm::ESHandle<TrackerGeometry> geomHandle_;
  edm::ESWatcher<TrackerDigiGeometryRecord> theTkDigiGeomWatcher;
  if (theTkDigiGeomWatcher.check(iSetup)) {
    iSetup.get<TrackerDigiGeometryRecord>().get(geomType_, geomHandle_);
  }  
  if (!geomHandle_.isValid()) return;

  const TrackerTopology* const tTopo = tTopoHandle_.product();  
  const TrackerGeometry* tGeom = geomHandle_.product();

  // SiStripClusters
  edm::Handle<edmNew::DetSetVector<SiStripCluster>> clusterHandle;
  iEvent.getByToken(clusterToken_, clusterHandle);

  // Primary vertex collection
  edm::Handle<reco::VertexCollection> vertexColl;
  iEvent.getByToken(vertexToken_, vertexColl);
  if (vertexColl->size() > 0) {
    const reco::Vertex& pv = (*vertexColl)[0];

    // Beam spot
    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(bsToken_, beamSpot);

    // Track collection
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(trackToken_, tracks);

    // Access PU information
    double wfac = 1.0;  // for data
    if (!iEvent.isRealData()) {
      edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
      iEvent.getByToken(puSummaryToken_, PupInfo);
     
      if (verbose_) edm::LogInfo("RecoTrackMonitor") << "nPUColl = " << PupInfo->size();
      for (auto const& v: *PupInfo) {
	int bx = v.getBunchCrossing();
	if (bx == 0) {
	  int ntrueInt = v.getTrueNumInteractions();
	  int nVertex = (vertexColl.isValid() ? vertexColl->size() : 0);
	  if (doPUCorrection_) {
	    if (nVertex > -1 && nVertex < int(vpu_.size())) wfac = vpu_.at(nVertex);
	    else wfac = 0.0;
	  }
	}
      }
    }
    if (verbose_) edm::LogInfo("RecoTrackMonitor") << "PU reweight factor = " << wfac;
    if (verbose_) std::cout << "PU scale factor" << wfac << std::endl;

    if (!vertexColl.isValid())
      edm::LogError("RecoTrackMonitor") << "Error! Failed to get reco::Vertex Collection, " << vertexTag_;
    
    int ntracks = 0;
    if (tracks.isValid()) {
      edm::LogInfo("RecoTrackMonitor") << "Total # of Tracks: " << tracks->size();
      if (verbose_) edm::LogInfo("RecoTrackMonitor") <<"Total # of Tracks: " << tracks->size();
      reco::Track::TrackQuality quality = reco::Track::qualityByName(trackQuality_);
      if (verbose_) std::cout <<"Total # of Tracks: " << tracks->size() << std::endl;
    
      for (auto const& track: *tracks) {
	if (!track.quality(quality)) continue;
	++ntracks;     
 
	double eta = track.eta();

	// From here, the purpose is to access RECO instead of AOD
	double residualXPB = 0, residualXPE = 0, residualXTIB = 0, residualXTOB = 0, residualXTEC = 0, residualXTID = 0;
	double residualYPB = 0, residualYPE = 0, residualYTIB = 0, residualYTOB = 0, residualYTEC = 0, residualYTID = 0;

	int nStripTIBS = 0, nStripTOBS = 0, nStripTECS = 0, nStripTIDS = 0;
	int nStripTIBD = 0, nStripTOBD = 0, nStripTECD = 0, nStripTIDD = 0;

	reco::TrackResiduals residuals = track.residuals();
	int i = 0;
	for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it,++i) {
	  const TrackingRecHit& hitn = (**it);   
     	  if (hitn.isValid()) {
	    DetId detId = hitn.geographicalId();
	    int subdetId = hitn.geographicalId().subdetId();

            const GeomDet *geomDet = tGeom->idToDet(detId);
	    if (!geomDet) continue;
	    Global3DPoint gPos = geomDet->surface().toGlobal(hitn.localPosition());
	    //   std::cout << "z position of hit: " << gPos.z() << std::endl;

	    int layer = tTopo->layer(detId);
            int lbin = layer - 1;             
            if (subdetId == PixelSubdetector::PixelBarrel || subdetId == PixelSubdetector::PixelEndcap) {
	      hitmap_pixelH_->Fill(gPos.eta(), gPos.phi());
              hitEta_pixelH_->Fill(gPos.eta());
	      hitPhi_pixelH_->Fill(gPos.phi());
	      hitmap2_pixelH_-> Fill(gPos.z(), gPos.phi());
              if (subdetId == PixelSubdetector::PixelBarrel) {
		//std::cout << "PixB: " << layer << std::endl; 
		hitmap_PixB_LayerH_[lbin]->Fill(gPos.eta(), gPos.phi());
		hitEta_PixB_LayerH_[lbin]->Fill(gPos.eta());
		hitPhi_PixB_LayerH_[lbin]->Fill(gPos.phi());
		hitmap2_PixB_LayerH_[lbin]->Fill(gPos.z(), gPos.phi());		
              }
	    }
	    else if (subdetId == StripSubdetector::TIB || subdetId == StripSubdetector::TOB || 
                     subdetId == StripSubdetector::TEC || subdetId == StripSubdetector::TID) {
	      hitmap_stripH_->Fill(gPos.eta(), gPos.phi());
	      hitEta_stripH_->Fill(gPos.eta());
	      hitPhi_stripH_->Fill(gPos.phi());
	      hitmap2_stripH_->Fill(gPos.z(), gPos.phi());
	      if (subdetId == StripSubdetector::TIB) {
		//std::cout << "TIB: " << layer << std::endl; 
		hitmap_TIB_LayerH_[lbin]->Fill(gPos.eta(), gPos.phi());
		hitEta_TIB_LayerH_[lbin]->Fill(gPos.eta());
		hitPhi_TIB_LayerH_[lbin]->Fill(gPos.phi());
		hitmap2_TIB_LayerH_[lbin]->Fill(gPos.z(), gPos.phi());
              }
	      if (subdetId == StripSubdetector::TOB) {
		//std::cout << "TOB: " << layer << std::endl; 
		hitmap_TOB_LayerH_[lbin]->Fill(gPos.eta(), gPos.phi());
		hitEta_TOB_LayerH_[lbin]->Fill(gPos.eta());
		hitPhi_TOB_LayerH_[lbin]->Fill(gPos.phi());
		hitmap2_TOB_LayerH_[lbin]->Fill(gPos.z(), gPos.phi());
              }
	    }
	    if      (subdetId == PixelSubdetector::PixelBarrel) {
	      residualXPB  = residuals.residualX(i);
	      residualYPB  = residuals.residualY(i);
	    }
	    else if (subdetId == PixelSubdetector::PixelEndcap) {
	      residualXPE  = residuals.residualX(i);
	      residualYPE  = residuals.residualY(i);
	    }
	    else if (subdetId == StripSubdetector::TIB) {
	      residualXTIB = residuals.residualX(i);
	      residualYTIB = residuals.residualY(i);
	    }
	    else if (subdetId == StripSubdetector::TOB) {
	      residualXTOB = residuals.residualX(i);
	      residualYTOB = residuals.residualY(i);
	    }
	    else if (subdetId == StripSubdetector::TEC) {
	      residualXTEC = residuals.residualX(i);
	      residualYTEC = residuals.residualY(i);
	    }
	    else if (subdetId == StripSubdetector::TID) {
	      residualXTID = residuals.residualX(i);
	      residualYTID = residuals.residualY(i);
	    }
	    
	    if (hitn.geographicalId().det() == DetId::Tracker) {
	      // Find on-track clusters
	      processHit(hitn, iSetup, tkGeom, wfac);
            
	      const DetId detId(hitn.geographicalId());
	      const SiStripDetId stripId(detId);
	      if (0) std::cout << "Hit Dimension: " << hitn.dimension()
			       << ", isGlued: " << stripId.glued()
			       << ", isStereo: " << stripId.stereo()
			       << std::endl;
            
	      if (stripId.glued()) {
		if      (subdetId == StripSubdetector::TIB) ++nStripTIBD;
		else if (subdetId == StripSubdetector::TOB) ++nStripTOBD;
		else if (subdetId == StripSubdetector::TEC) ++nStripTECD;
		else if (subdetId == StripSubdetector::TID) ++nStripTIDD;
	      }
	      else {
		if      (subdetId == StripSubdetector::TIB) ++nStripTIBS;
		else if (subdetId == StripSubdetector::TOB) ++nStripTOBS;
		else if (subdetId == StripSubdetector::TEC) ++nStripTECS;
		else if (subdetId == StripSubdetector::TID) ++nStripTIDS;
	      }
	    }
	  }
	  residualXPBH_->Fill(residualXPB);
	  residualXPEH_->Fill(residualXPE);
	  residualXTIBH_->Fill(residualXTIB);
	  residualXTOBH_->Fill(residualXTOB);
	  residualXTECH_->Fill(residualXTEC);
	  residualXTIDH_->Fill(residualXTID);
	  residualYPBH_->Fill(residualYPB);
	  residualYPEH_->Fill(residualYPE);
	  residualYTIBH_->Fill(residualYTIB);
	  residualYTOBH_->Fill(residualYTOB);
	  residualYTECH_->Fill(residualYTEC);
	  residualYTIDH_->Fill(residualYTID);
	}
	nHitsTIBSVsEtaH_->Fill(eta, nStripTIBS);
	nHitsTOBSVsEtaH_->Fill(eta, nStripTOBS);
	nHitsTECSVsEtaH_->Fill(eta, nStripTECS);
	nHitsTIDSVsEtaH_->Fill(eta, nStripTIDS);
	nHitsStripSVsEtaH_->Fill(eta, nStripTIBS+nStripTOBS+nStripTECS+nStripTIDS);
                                                                     
	nHitsTIBDVsEtaH_->Fill(eta, nStripTIBD);
	nHitsTOBDVsEtaH_->Fill(eta, nStripTOBD);
	nHitsTECDVsEtaH_->Fill(eta, nStripTECD);
	nHitsTIDDVsEtaH_->Fill(eta, nStripTIDD);
	nHitsStripDVsEtaH_->Fill(eta, nStripTIBD+nStripTOBD+nStripTECD+nStripTIDD);
      }
    }
    else {
      edm::LogError("RecoTrackMonitor") << "Error! Failed to get reco::Track collection, " << trackTag_;
    }
    // off-track cluster properties                                                                                                                       
    if (clusterHandle.isValid()) processClusters(clusterHandle, iSetup, tkGeom, wfac);
    else edm::LogError("RecoTrackMonitor") << "ClusterCollection " << clusterTag_ << " not valid!!" << std::endl;

    if (verbose_) std::cout << "Ends RecoTrackMonitor successfully" << std::endl;
  }
}

void RecoTrackMonitor::processClusters(const edm::Handle<edmNew::DetSetVector<SiStripCluster>>& clusterHandle, 
                                       edm::EventSetup const& iSetup, const TrackerGeometry& tkGeom, double wfac) {
  // Loop on Dets
  for (edmNew::DetSetVector<SiStripCluster>::const_iterator dsvit  = clusterHandle->begin(); 
                                                            dsvit != clusterHandle->end();
                                                          ++dsvit) {
    uint32_t detId = dsvit->id();
    std::map<uint32_t, std::set<const SiStripCluster*> >::iterator jt = clusterMap_.find(detId);
    bool detid_found = (jt != clusterMap_.end()) ? true : false;
    
    // Loop on Clusters
    for (edmNew::DetSet<SiStripCluster>::const_iterator clusit  = dsvit->begin(); 
     	                                                clusit != dsvit->end();         
	                                              ++clusit) {
      if (detid_found) {
	std::set<const SiStripCluster*>& s = jt->second;
	if (s.find(&*clusit) != s.end()) continue;
      }
      
      SiStripClusterInfo info(*clusit, iSetup, detId);
      float charge = info.charge();
      float width = info.width();
      
      const GeomDetUnit* detUnit = tkGeom.idToDetUnit(detId);
      float thickness =  detUnit->surface().bounds().thickness(); // unit cm
      if (thickness > 0.035) {
	hOffTrkClusChargeThickH_->Fill(charge, wfac);
	hOffTrkClusWidthThickH_->Fill(width, wfac);
      }
      else {
	hOffTrkClusChargeThinH_->Fill(charge, wfac);
	hOffTrkClusWidthThinH_->Fill(width, wfac);
      }
    }
  }
}
void RecoTrackMonitor::processHit(const TrackingRecHit& recHit, edm::EventSetup const& iSetup, const TrackerGeometry& tkGeom, double wfac)
{
  uint32_t detid = recHit.geographicalId();
  const GeomDetUnit* detUnit = tkGeom.idToDetUnit(detid);
  float thickness =  detUnit->surface().bounds().thickness(); // unit cm      

  auto const& thit = static_cast<BaseTrackerRecHit const&>(recHit);
  if (!thit.isValid()) return;

  auto const& clus = thit.firstClusterRef();
  if (!clus.isValid()) return;
  if (!clus.isStrip()) return;

  if (thit.isMatched()) {
    const SiStripMatchedRecHit2D& matchedHit = dynamic_cast<const SiStripMatchedRecHit2D&>(recHit);

    auto& clusterM = matchedHit.monoCluster();
    SiStripClusterInfo infoM(clusterM, iSetup, detid);
    if (thickness > 0.035) {
      hOnTrkClusChargeThickH_->Fill(infoM.charge(), wfac);
      hOnTrkClusWidthThickH_->Fill(infoM.width(), wfac);
    }
    else {
      hOnTrkClusChargeThinH_->Fill(infoM.charge(), wfac);
      hOnTrkClusWidthThinH_->Fill(infoM.width(), wfac);
    }
    addClusterToMap(detid, &clusterM);

    auto& clusterS = matchedHit.stereoCluster();
    SiStripClusterInfo infoS(clusterS, iSetup, detid);
    if (thickness > 0.035) {
      hOnTrkClusChargeThickH_->Fill(infoS.charge(), wfac);
      hOnTrkClusWidthThickH_->Fill(infoS.width(), wfac );
    }
    else {
      hOnTrkClusChargeThinH_->Fill(infoS.charge(), wfac);
      hOnTrkClusWidthThinH_->Fill(infoS.width(), wfac);
    }
    addClusterToMap(detid, &clusterS);
  }
  else {
    auto& cluster = clus.stripCluster();
    SiStripClusterInfo info(cluster, iSetup, detid);
    if (thickness > 0.035) {
      hOnTrkClusChargeThickH_->Fill(info.charge(), wfac);
      hOnTrkClusWidthThickH_->Fill(info.width(), wfac);
    }
    else {
      hOnTrkClusChargeThinH_->Fill(info.charge(), wfac);
      hOnTrkClusWidthThinH_->Fill(info.width(), wfac);
    }
    addClusterToMap(detid, &cluster);
  }
}

void RecoTrackMonitor::addClusterToMap(uint32_t detid, const SiStripCluster* cluster) {
  std::map<uint32_t, std::set<const SiStripCluster*> >::iterator it = clusterMap_.find(detid);
  if (it == clusterMap_.end()) {
    std::set<const SiStripCluster*> s;
    s.insert(cluster);
    clusterMap_.insert(std::pair<uint32_t, std::set<const SiStripCluster*> >(detid, s));
  }
  else {
    std::set<const SiStripCluster*>& s = it->second;
    s.insert(cluster);
  }
}

void RecoTrackMonitor::endLuminosityBlock(edm::LuminosityBlock const& lumiBlock, edm::EventSetup const& eSetup){
}
// Define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoTrackMonitor);
