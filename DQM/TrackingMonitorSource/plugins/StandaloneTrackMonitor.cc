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
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "RecoLocalTracker/SiStripClusterizer/interface/SiStripClusterInfo.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DQM/TrackingMonitorSource/interface/StandaloneTrackMonitor.h"
#include "DQM/TrackingMonitor/interface/TrackBuildingAnalyzer.h"
#include <string>

namespace {
  bool trackSelected(unsigned char mask, unsigned char qual) { return mask & 1 << qual; } 
}

// -----------------------------
// constructors and destructor
// -----------------------------
StandaloneTrackMonitor::StandaloneTrackMonitor(const edm::ParameterSet& ps): 
  parameters_(ps),
  moduleName_(parameters_.getUntrackedParameter<std::string>("moduleName", "StandaloneTrackMonitor")),
  folderName_(parameters_.getUntrackedParameter<std::string>("folderName", "highPurityTracks")),
  trackTag_(parameters_.getUntrackedParameter<edm::InputTag>("trackInputTag", edm::InputTag("generalTracks"))),
  bsTag_(parameters_.getUntrackedParameter<edm::InputTag>("offlineBeamSpot", edm::InputTag("offlineBeamSpot"))),
  vertexTag_(parameters_.getUntrackedParameter<edm::InputTag>("vertexTag", edm::InputTag("offlinePrimaryVertices"))),
  puSummaryTag_(parameters_.getUntrackedParameter<edm::InputTag>("puTag", edm::InputTag("addPileupInfo"))),
  clusterTag_(parameters_.getUntrackedParameter<edm::InputTag>("clusterTag", edm::InputTag("siStripClusters"))),
  trackToken_(consumes<reco::TrackCollection>(trackTag_)),
  bsToken_(consumes<reco::BeamSpot>(bsTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  puSummaryToken_(consumes<std::vector<PileupSummaryInfo> >(puSummaryTag_)),
  clusterToken_(consumes<edmNew::DetSetVector<SiStripCluster> >(clusterTag_)),
  trackQuality_(parameters_.getUntrackedParameter<std::string>("trackQuality", "highPurity")),
  doPUCorrection_(parameters_.getUntrackedParameter<bool>("doPUCorrection", false)),
  isMC_(parameters_.getUntrackedParameter<bool>("isMC", false)),
  haveAllHistograms_(parameters_.getUntrackedParameter<bool>("haveAllHistograms", false)),
  puScaleFactorFile_(parameters_.getUntrackedParameter<std::string>("puScaleFactorFile", "PileupScaleFactor.root")),
  mvaProducers_(parameters_.getUntrackedParameter<std::vector<std::string> >("MVAProducers")),
  mvaTrackTag_(parameters_.getUntrackedParameter<edm::InputTag>("TrackProducerForMVA")),
  mvaTrackToken_(consumes<edm::View<reco::Track> >(mvaTrackTag_)),
  tcProducer_(parameters_.getUntrackedParameter<edm::InputTag>("TCProducer")),
  algoName_(parameters_.getUntrackedParameter<std::string>("AlgoName")),
  verbose_(parameters_.getUntrackedParameter<bool>("verbose", false))
{
  for (auto v: mvaProducers_) {
    mvaQualityTokens_.push_back(std::make_tuple(consumes<MVACollection>(edm::InputTag(v, "MVAValues")),
						consumes<QualityMaskCollection>(edm::InputTag(v, "QualityMasks")))
			       );
  }

  // for MC only
  nVtxH_ = nullptr;
  nVertexH_ = nullptr;
  bunchCrossingH_ = nullptr;
  nPUH_ = nullptr;
  trueNIntH_ = nullptr;

  nLostHitsVspTH_ = nullptr;
  nLostHitsVsEtaH_ = nullptr;
  nLostHitsVsCosThetaH_ = nullptr;
  nLostHitsVsPhiH_ = nullptr;

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

  nValidHitsVspTH_ = nullptr;
  nValidHitsVsEtaH_ = nullptr;
  nValidHitsVsCosThetaH_ = nullptr;
  nValidHitsVsPhiH_ = nullptr;
  nValidHitsVsnVtxH_ = nullptr;

  nValidHitsPixVsEtaH_ = nullptr;
  nValidHitsPixBVsEtaH_ = nullptr;
  nValidHitsPixEVsEtaH_ = nullptr;
  nValidHitsStripVsEtaH_ = nullptr;
  nValidHitsTIBVsEtaH_ = nullptr;
  nValidHitsTOBVsEtaH_ = nullptr;
  nValidHitsTECVsEtaH_ = nullptr;
  nValidHitsTIDVsEtaH_ = nullptr;

  nValidHitsPixVsPhiH_ = nullptr;
  nValidHitsPixBVsPhiH_ = nullptr;
  nValidHitsPixEVsPhiH_ = nullptr;
  nValidHitsStripVsPhiH_ = nullptr;
  nValidHitsTIBVsPhiH_ = nullptr;
  nValidHitsTOBVsPhiH_ = nullptr;
  nValidHitsTECVsPhiH_ = nullptr;
  nValidHitsTIDVsPhiH_ = nullptr;

  nLostHitsPixVsEtaH_ = nullptr;
  nLostHitsPixBVsEtaH_ = nullptr;
  nLostHitsPixEVsEtaH_ = nullptr;
  nLostHitsStripVsEtaH_ = nullptr;
  nLostHitsTIBVsEtaH_ = nullptr;
  nLostHitsTOBVsEtaH_ = nullptr;
  nLostHitsTECVsEtaH_ = nullptr;
  nLostHitsTIDVsEtaH_ = nullptr;

  nLostHitsPixVsPhiH_ = nullptr;
  nLostHitsPixBVsPhiH_ = nullptr;
  nLostHitsPixEVsPhiH_ = nullptr;
  nLostHitsStripVsPhiH_ = nullptr;
  nLostHitsTIBVsPhiH_ = nullptr;
  nLostHitsTOBVsPhiH_ = nullptr;
  nLostHitsTECVsPhiH_ = nullptr;
  nLostHitsTIDVsPhiH_ = nullptr;

  hOnTrkClusChargeThinH_ = nullptr;
  hOnTrkClusWidthThinH_ = nullptr;
  hOnTrkClusChargeThickH_ = nullptr;
  hOnTrkClusWidthThickH_ = nullptr;
  
  hOffTrkClusChargeThinH_ = nullptr;
  hOffTrkClusWidthThinH_ = nullptr;
  hOffTrkClusChargeThickH_ = nullptr;
  hOffTrkClusWidthThickH_ = nullptr;

  // Read pileup weight factors
  if (isMC_ && doPUCorrection_) {
    vpu_.clear();
    TFile* f1 = TFile::Open(puScaleFactorFile_.c_str());
    TH1F* h1 = dynamic_cast<TH1F*>(f1->Get("pileupweight"));
    for (int i = 1; i <= h1->GetNbinsX(); ++i) vpu_.push_back(h1->GetBinContent(i));
    f1->Close();
  }
}

void StandaloneTrackMonitor::bookHistograms(DQMStore::IBooker &ibook, edm::Run const& iRun, edm::EventSetup const& iSetup) {
  
  edm::ParameterSet TrackEtaHistoPar = parameters_.getParameter<edm::ParameterSet>("trackEtaH");
  edm::ParameterSet TrackPtHistoPar = parameters_.getParameter<edm::ParameterSet>("trackPtH");
  edm::ParameterSet TrackMVAHistoPar = parameters_.getParameter<edm::ParameterSet>("trackMVAH");

  std::string currentFolder = moduleName_ + "/" + folderName_ ;
  ibook.setCurrentFolder(currentFolder.c_str());

  std::vector<std::string> mvaProducers = parameters_.getUntrackedParameter<std::vector<std::string> >("MVAProducers");
  edm::InputTag tcProducer = parameters_.getUntrackedParameter<edm::InputTag>("TCProducer");
  //MVA
  const std::string& CatagoryName = algoName_;

  for (size_t i = 1, end = mvaProducers.size(); i <= end; ++i) {
    auto num = std::to_string(i);
    std::string pfix;

    if (i == 1) {
      trackMVAsHP.push_back(nullptr);
      trackMVAsHPVsPtProfile.push_back(nullptr);
      trackMVAsHPVsEtaProfile.push_back(nullptr);
    } else {
      pfix = " (not loose-selected)";
      std::string pfix2 = " (not HP-selected)";
      histname = "TrackMVA" + num + "HP_" + tcProducer.label() + "_";
      trackMVAsHP.push_back(ibook.book1D(histname + CatagoryName, histname + CatagoryName + pfix2, TrackMVAHistoPar.getParameter<int32_t>("Xbins"), TrackMVAHistoPar.getParameter<double>("Xmin"), TrackMVAHistoPar.getParameter<double>("Xmax")));
      trackMVAsHP.back()->setAxisTitle("Track selection MVA" + num, 1);
      trackMVAsHP.back()->setAxisTitle("Number of tracks", 2);

      histname = "TrackMVA" + num + "HPVsPtProfile_" + tcProducer.label() + "_";
			    trackMVAsHPVsPtProfile.push_back(ibook.bookProfile(histname + CatagoryName,histname + CatagoryName + pfix2,TrackPtHistoPar.getParameter<int32_t>("Xbins"),TrackPtHistoPar.getParameter<double>("Xmin"),TrackPtHistoPar.getParameter<double>("Xmax"),TrackMVAHistoPar.getParameter<int32_t>("Xbins"),TrackMVAHistoPar.getParameter<double>("Xmin"),TrackMVAHistoPar.getParameter<double>("Xmax")));
      trackMVAsHPVsPtProfile.back()->setAxisTitle("Track p_{T} (GeV/c)", 1);
      trackMVAsHPVsPtProfile.back()->setAxisTitle("Track selection MVA" + num, 2);

      histname = "TrackMVA" + num + "HPVsEtaProfile_" + tcProducer.label() + "_";
			    trackMVAsHPVsEtaProfile.push_back(ibook.bookProfile(histname + CatagoryName, histname + CatagoryName + pfix2, TrackEtaHistoPar.getParameter<int32_t>("Xbins"), TrackEtaHistoPar.getParameter<double>("Xmin"), TrackEtaHistoPar.getParameter<double>("Xmax"), TrackMVAHistoPar.getParameter<int32_t>("Xbins"), TrackMVAHistoPar.getParameter<double>("Xmin"), TrackMVAHistoPar.getParameter<double>("Xmax")));
      trackMVAsHPVsEtaProfile.back()->setAxisTitle("Track #eta", 1);
      trackMVAsHPVsEtaProfile.back()->setAxisTitle("Track selection MVA" + num, 2);
    }

    histname = "TrackMVA" + num + "_" + tcProducer.label() + "_";
      trackMVAs.push_back(ibook.book1D(histname + CatagoryName, histname + CatagoryName + pfix, TrackMVAHistoPar.getParameter<int32_t>("Xbins"), TrackMVAHistoPar.getParameter<double>("Xmin"), TrackMVAHistoPar.getParameter<double>("Xmax")));
    trackMVAs.back()->setAxisTitle("Track selection MVA" + num, 1);
    trackMVAs.back()->setAxisTitle("Number of tracks", 2);

    histname = "TrackMVA" + num + "VsPtProfile_" + tcProducer.label() + "_";
    trackMVAsVsPtProfile.push_back(ibook.bookProfile(histname + CatagoryName,histname + CatagoryName + pfix,TrackPtHistoPar.getParameter<int32_t>("Xbins"),TrackPtHistoPar.getParameter<double>("Xmin"),TrackPtHistoPar.getParameter<double>("Xmax"),TrackMVAHistoPar.getParameter<int32_t>("Xbins"),TrackMVAHistoPar.getParameter<double>("Xmin"),TrackMVAHistoPar.getParameter<double>("Xmax")));
    trackMVAsVsPtProfile.back()->setAxisTitle("Track p_{T} (GeV/c)", 1);
    trackMVAsVsPtProfile.back()->setAxisTitle("Track selection MVA" + num, 2);

    histname = "TrackMVA" + num + "VsEtaProfile_" + tcProducer.label() + "_";
    trackMVAsVsEtaProfile.push_back(ibook.bookProfile(histname + CatagoryName, histname + CatagoryName + pfix, TrackEtaHistoPar.getParameter<int32_t>("Xbins"), TrackEtaHistoPar.getParameter<double>("Xmin"), TrackEtaHistoPar.getParameter<double>("Xmax"), TrackMVAHistoPar.getParameter<int32_t>("Xbins"), TrackMVAHistoPar.getParameter<double>("Xmin"), TrackMVAHistoPar.getParameter<double>("Xmax")));
    trackMVAsVsEtaProfile.back()->setAxisTitle("Track #eta", 1);
    trackMVAsVsEtaProfile.back()->setAxisTitle("Track selection MVA" + num, 2);
  }


  // MVA


  // The following are common with the official tool
  if (haveAllHistograms_) {
    trackEtaH_ = ibook.book1D("trackEta", "Track Eta", 
			      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
			      TrackEtaHistoPar.getParameter<double>("Xmin"), 
			      TrackEtaHistoPar.getParameter<double>("Xmax"));
    
    trackEtaerrH_ = ibook.book1D("trackEtaerr", "Track Eta Error", 50,0.0,1.0);
    trackCosThetaH_ = ibook.book1D("trackCosTheta", "Track Cos(Theta)", 50,-1.0,1.0);
    trackThetaerrH_ = ibook.book1D("trackThetaerr", "Track Theta Error", 50,0.0,1.0);
    trackPhiH_ = ibook.book1D("trackPhi", "Track Phi", 70,-3.5,3.5);
    trackPhierrH_ = ibook.book1D("trackPhierr", "Track Phi Error", 50,0.0,1.0);
    
    trackPH_ = ibook.book1D("trackP", "Track 4-momentum", 50,0.0,10.0);
    trackPtH_ = ibook.book1D("trackPt", "Track Pt", 
			     TrackPtHistoPar.getParameter<int32_t>("Xbins"),
			     TrackPtHistoPar.getParameter<double>("Xmin"),
			     TrackPtHistoPar.getParameter<double>("Xmax"));
    trackPtUpto2GeVH_ = ibook.book1D("trackPtUpto2GeV", "Track Pt upto 2GeV",100,0,2.0);
    trackPtOver10GeVH_ = ibook.book1D("trackPtOver10GeV","Track Pt greater than 10 GeV",100,0,100.0);
    trackPterrH_ = ibook.book1D("trackPterr", "Track Pt Error",100,0.0,100.0);
    trackqOverpH_ = ibook.book1D("trackqOverp", "q Over p",40,-10.0,10.0);
    trackqOverperrH_ = ibook.book1D("trackqOverperr","q Over p Error",50,0.0,25.0);
    trackChargeH_ = ibook.book1D("trackCharge", "Track Charge", 50, -5, 5);
    trackChi2H_ = ibook.book1D("trackChi2", "Chi2",100,0.0,100.0);
    tracknDOFH_ = ibook.book1D("tracknDOF", "nDOF",100,0.0,100.0);
    trackChi2ProbH_ = ibook.book1D("trackChi2Prob", "Chi2prob",50,0.0,1.0);
    trackChi2oNDFH_ = ibook.book1D("trackChi2oNDF", "Chi2oNDF",100,0.0,100.0);
    trackd0H_ = ibook.book1D("trackd0", "Track d0",100,-1,1);
    trackChi2bynDOFH_ = ibook.book1D("trackChi2bynDOF", "Chi2 Over nDOF",100,0.0,10.0);
    trackalgoH_ = ibook.book1D("trackalgo", "Track Algo",50,0.0,50.0);
    trackorigalgoH_ = ibook.book1D("trackorigalgo", "Track Original Algo",50,0.0,50.0);
    DistanceOfClosestApproachToPVH_ = ibook.book1D("DistanceOfClosestApproachToPV", "DistanceOfClosestApproachToPV",100,-0.1,0.1);
    DistanceOfClosestApproachToPVVsPhiH_ = ibook.bookProfile("DistanceOfClosestApproachToPVVsPhi", "DistanceOfClosestApproachToPVVsPhi",100,-3.5,3.5,0.0,0.0,"g");
    xPointOfClosestApproachVsZ0wrtPVH_ = ibook.bookProfile("xPointOfClosestApproachVsZ0wrtPV", "xPointOfClosestApproachVsZ0wrtPV",120,-60,60,0.0,0.0,"g");
    yPointOfClosestApproachVsZ0wrtPVH_ = ibook.bookProfile("yPointOfClosestApproachVsZ0wrtPV", "yPointOfClosestApproachVsZ0wrtPV",120,-60,60,0.0,0.0,"g");
    
    sip3dToPVH_ = ibook.book1D("sip3dToPV", "signed Impact Point 3d To PV",200,-10,10);
    sip2dToPVH_ = ibook.book1D("sip2dToPV", "signed Impact Point 2d To PV",200,-10,10);
    sipDxyToPVH_ = ibook.book1D("sipDxyToPV", "signed Impact Point dxy To PV",100,-10,10);
    sipDzToPVH_ = ibook.book1D("sipDzToPV", "signed Impact Point dz To PV",100,-10,10);
    
    nvalidTrackerHitsH_ = ibook.book1D("nvalidTrackerhits", "No. of Valid Tracker Hits",45,0.5,45.5);
    nvalidPixelHitsH_ = ibook.book1D("nvalidPixelHits", "No. of Valid Hits in Pixel",7,-0.5,6.5);
    nvalidPixelBHitsH_ = ibook.book1D("nvalidPixelBarrelHits", "No. of Valid Hits in Pixel Barrel",5,-0.5,4.5);
    nvalidPixelEHitsH_ = ibook.book1D("nvalidPixelEndcapHits", "No. of Valid Hits in Pixel Endcap",5,-0.5,5.5);
    nvalidStripHitsH_ = ibook.book1D("nvalidStripHits", "No. of Valid Hits in Strip",35,-0.5,34.5);
    nvalidTIBHitsH_ = ibook.book1D("nvalidTIBHits", "No. of Valid Hits in Strip TIB",5,-0.5,4.5);
    nvalidTOBHitsH_ = ibook.book1D("nvalidTOBHits", "No. of Valid Hits in Strip TOB",10,-0.5,9.5);
    nvalidTIDHitsH_ = ibook.book1D("nvalidTIDHits", "No. of Valid Hits in Strip TID",5,-0.5,4.5);
    nvalidTECHitsH_ = ibook.book1D("nvalidTECHits", "No. of Valid Hits in Strip TEC",10,-0.5,9.5);

    nlostTrackerHitsH_ = ibook.book1D("nlostTrackerhits", "No. of Lost Tracker Hits",15,-0.5,14.5);
    nlostPixelHitsH_ = ibook.book1D("nlostPixelHits", "No. of Lost Hits in Pixel",8,-0.5,7.5);
    nlostPixelBHitsH_ = ibook.book1D("nlostPixelBarrelHits", "No. of Lost Hits in Pixel Barrel",5,-0.5,4.5);
    nlostPixelEHitsH_ = ibook.book1D("nlostPixelEndcapHits", "No. of Lost Hits in Pixel Endcap",4,-0.5,3.5);
    nlostStripHitsH_ = ibook.book1D("nlostStripHits", "No. of Lost Hits in Strip",10,-0.5,9.5);
    nlostTIBHitsH_ = ibook.book1D("nlostTIBHits", "No. of Lost Hits in Strip TIB",5,-0.5,4.5);
    nlostTOBHitsH_ = ibook.book1D("nlostTOBHits", "No. of Lost Hits in Strip TOB",10,-0.5,9.5);
    nlostTIDHitsH_ = ibook.book1D("nlostTIDHits", "No. of Lost Hits in Strip TID",5,-0.5,4.5);
    nlostTECHitsH_ = ibook.book1D("nlostTECHits", "No. of Lost Hits in Strip TEC",10,-0.5,9.5);    

    trkLayerwithMeasurementH_ = ibook.book1D("trkLayerwithMeasurement", "No. of Layers per Track",20,0.0,20.0);
    pixelLayerwithMeasurementH_ = ibook.book1D("pixelLayerwithMeasurement", "No. of Pixel Layers per Track",10,0.0,10.0);
    pixelBLayerwithMeasurementH_ = ibook.book1D("pixelBLayerwithMeasurement", "No. of Pixel Barrel Layers per Track",5,0.0,5.0);
    pixelELayerwithMeasurementH_ = ibook.book1D("pixelELayerwithMeasurement", "No. of Pixel Endcap Layers per Track",5,0.0,5.0);
    stripLayerwithMeasurementH_ = ibook.book1D("stripLayerwithMeasurement", "No. of Strip Layers per Track",20,0.0,20.0);
    stripTIBLayerwithMeasurementH_ = ibook.book1D("stripTIBLayerwithMeasurement", "No. of Strip TIB Layers per Track",10,0.0,10.0);
    stripTOBLayerwithMeasurementH_ = ibook.book1D("stripTOBLayerwithMeasurement", "No. of Strip TOB Layers per Track",10,0.0,10.0);
    stripTIDLayerwithMeasurementH_ = ibook.book1D("stripTIDLayerwithMeasurement", "No. of Strip TID Layers per Track",5,0.0,5.0);
    stripTECLayerwithMeasurementH_ = ibook.book1D("stripTECLayerwithMeasurement", "No. of Strip TEC Layers per Track",15,0.0,15.0);

    nlostHitsH_ = ibook.book1D("nlostHits", "No. of Lost Hits",10,-0.5,9.5);

    beamSpotXYposH_ = ibook.book1D("beamSpotXYpos", "XY position of beam spot",40,-4.0,4.0);
    beamSpotXYposerrH_ = ibook.book1D("beamSpotXYposerr", "Error in XY position of beam spot",20,0.0,4.0);
    beamSpotZposH_ = ibook.book1D("beamSpotZpos", "Z position of beam spot",100,-20.0,20.0);
    beamSpotZposerrH_ = ibook.book1D("beamSpotZposerr", "Error in Z position of beam spot", 50, 0.0, 5.0);

    vertexXposH_ = ibook.book1D("vertexXpos", "Vertex X position", 50, -1.0, 1.0);
    vertexYposH_ = ibook.book1D("vertexYpos", "Vertex Y position", 50, -1.0, 1.0);
    vertexZposH_ = ibook.book1D("vertexZpos", "Vertex Z position", 100,-20.0,20.0);
    nVertexH_ = ibook.book1D("nVertex", "# of vertices", 60, -0.5, 59.5);
    nVtxH_ = ibook.book1D("nVtx", "# of vtxs", 60, -0.5, 59.5);
    

    nMissingInnerHitBH_ = ibook.book1D("nMissingInnerHitB", "No. missing inner hit per Track in Barrel", 6, -0.5, 5.5);
    nMissingInnerHitEH_ = ibook.book1D("nMissingInnerHitE", "No. missing inner hit per Track in Endcap", 6, -0.5, 5.5);
    nMissingOuterHitBH_ = ibook.book1D("nMissingOuterHitB", "No. missing outer hit per Track in Barrel", 11, -0.5, 10.5);
    nMissingOuterHitEH_ = ibook.book1D("nMissingOuterHitE", "No. missing outer hit per Track in Endcap", 11, -0.5, 10.5);

    residualXPBH_ = ibook.book1D("residualXPixelBarrel", "Residual in X in Pixel Barrel", 100, -0.1, 0.1);
    residualXPEH_ = ibook.book1D("residualXPixelEndcap", "Residual in X in Pixel Endcap", 100, -0.04, 0.04);
    residualXTIBH_ = ibook.book1D("residualXStripTIB", "Residual in X in Strip TIB", 100, -0.05, 0.05);
    residualXTOBH_ = ibook.book1D("residualXStripTOB", "Residual in X in Strip TOB", 100, -0.04, 0.04);
    residualXTECH_ = ibook.book1D("residualXStripTEC", "Residual in X in Strip TEC", 100, -0.04, 0.04);
    residualXTIDH_ = ibook.book1D("residualXStripTID", "Residual in X in Strip TID", 100, -0.04, 0.04);
    residualYPBH_ = ibook.book1D("residualYPixelBarrel", "Residual in Y in Pixel Barrel", 100, -0.1, 0.1);
    residualYPEH_ = ibook.book1D("residualYPixelEndcap", "Residual in Y in Pixel Endcap", 100, -0.04, 0.04);
    residualYTIBH_ = ibook.book1D("residualYStripTIB", "Residual in Y in Strip TIB", 100, -4., 4.);
    residualYTOBH_ = ibook.book1D("residualYStripTOB", "Residual in Y in Strip TOB", 100, -4., 4.);
    residualYTECH_ = ibook.book1D("residualYStripTEC", "Residual in Y in Strip TEC", 100, -4., 4.);
    residualYTIDH_ = ibook.book1D("residualYStripTID", "Residual in Y in Strip TID", 100, -4., 4.);

    nTracksH_ = ibook.book1D("nTracks", "No. of Tracks", 100, -0.5, 999.5);
  }
  if (isMC_) {
    bunchCrossingH_ = ibook.book1D("bunchCrossing", "Bunch Crosssing", 60, 0, 60.0);
    nPUH_ = ibook.book1D("nPU", "No of Pileup", 100, 0, 100.0);
    trueNIntH_ = ibook.book1D("trueNInt", "True no of Interactions", 100, 0, 100.0);
  }
  // Exclusive histograms
  
  nLostHitByLayerH_ = ibook.book1D("nLostHitByLayer", "No. of Lost Hit per Layer", 29, 0.5, 29.5);

  nLostHitByLayerPixH_ = ibook.book1D("nLostHitByLayerPix", "No. of Lost Hit per Layer for Pixel detector", 7, 0.5, 7.5);

  nLostHitByLayerStripH_ = ibook.book1D("nLostHitByLayerStrip", "No. of Lost Hit per Layer for SiStrip detector", 22, 0.5, 22.5);

  //  nLostHitsVsnvtxH_ = ibook.bookProfile("nLostHitsVsnvtx", "Number of Lost Hits Vs nvtx",
  //					60,-0.5,59.5,0.0,0.0,"g");
  
  nLostHitsVspTH_ = ibook.bookProfile("nLostHitsVspT", "Number of Lost Hits Vs pT",
				      TrackPtHistoPar.getParameter<int32_t>("Xbins"),
				      TrackPtHistoPar.getParameter<double>("Xmin"),
				      TrackPtHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsVsEtaH_ = ibook.bookProfile("nLostHitsVsEta", "Number of Lost Hits Vs Eta", 
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"), 
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsVsCosThetaH_ = ibook.bookProfile("nLostHitsVsCosTheta", "Number of Lost Hits Vs Cos(Theta)",50,-1.0,1.0,0.0,0.0,"g");
  nLostHitsVsPhiH_ = ibook.bookProfile("nLostHitsVsPhi", "Number of Lost Hits Vs Phi",100,-3.5,3.5,0.0,0.0,"g");

  nHitsTIBSVsEtaH_ = ibook.bookProfile("nHitsTIBSVsEta", "Number of Hits in TIB Vs Eta (Single-sided)",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsTOBSVsEtaH_ = ibook.bookProfile("nHitsTOBSVsEta", "Number of Hits in TOB Vs Eta (Single-sided)",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");   
  nHitsTECSVsEtaH_ = ibook.bookProfile("nHitsTECSVsEta", "Number of Hits in TEC Vs Eta (Single-sided)",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsTIDSVsEtaH_ = ibook.bookProfile("nHitsTIDSVsEta", "Number of Hits in TID Vs Eta (Single-sided)",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  
  nHitsStripSVsEtaH_ = ibook.bookProfile("nHitsStripSVsEta", "Number of Strip Hits Vs Eta (Single-sided)",
					 TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					 TrackEtaHistoPar.getParameter<double>("Xmin"),
					 TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  
  nHitsTIBDVsEtaH_ = ibook.bookProfile("nHitsTIBDVsEta", "Number of Hits in TIB Vs Eta (Double-sided)",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsTOBDVsEtaH_ = ibook.bookProfile("nHitsTOBDVsEta", "Number of Hits in TOB Vs Eta (Double-sided)",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");   
  nHitsTECDVsEtaH_ = ibook.bookProfile("nHitsTECDVsEta", "Number of Hits in TEC Vs Eta (Double-sided)",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsTIDDVsEtaH_ = ibook.bookProfile("nHitsTIDDVsEta", "Number of Hits in TID Vs Eta (Double-sided)",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nHitsStripDVsEtaH_ = ibook.bookProfile("nHitsStripDVsEta", "Number of Strip Hits Vs Eta (Double-sided)",
					 TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					 TrackEtaHistoPar.getParameter<double>("Xmin"),
					 TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  
  nValidHitsVspTH_ = ibook.bookProfile("nValidHitsVspT", "Number of Valid Hits Vs pT",
				  TrackPtHistoPar.getParameter<int32_t>("Xbins"),
				  TrackPtHistoPar.getParameter<double>("Xmin"),
				  TrackPtHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsVsnVtxH_ = ibook.bookProfile("nValidHitsVsnVtx", "Number of Valid Hits Vs Number of Vertex", 100,0,50,0.0,0.0,"g");
  nValidHitsVsEtaH_ = ibook.bookProfile("nValidHitsVsEta", "Number of Hits Vs Eta", 
				   TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				   TrackEtaHistoPar.getParameter<double>("Xmin"),
				   TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  
  nValidHitsVsCosThetaH_ = ibook.bookProfile("nValidHitsVsCosTheta", "Number of Valid Hits Vs Cos(Theta)", 50,-1.0,1.0,0.0,0.0,"g");
  nValidHitsVsPhiH_ = ibook.bookProfile("nValidHitsVsPhi", "Number of Valid Hits Vs Phi", 100,-3.5,3.5,0.0,0.0,"g");
  
  nValidHitsPixVsEtaH_ = ibook.bookProfile("nValidHitsPixVsEta", "Number of Valid Hits in Pixel Vs Eta",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsPixBVsEtaH_ = ibook.bookProfile("nValidHitsPixBVsEta", "Number of Valid Hits in Pixel Barrel Vs Eta",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsPixEVsEtaH_ = ibook.bookProfile("nValidHitsPixEVsEta", "Number of Valid Hits in Pixel Endcap Vs Eta",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsStripVsEtaH_ = ibook.bookProfile("nValidHitsStripVsEta", "Number of Valid Hits in SiStrip Vs Eta",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsTIBVsEtaH_ = ibook.bookProfile("nValidHitsTIBVsEta", "Number of Valid Hits in TIB Vs Eta",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsTOBVsEtaH_ = ibook.bookProfile("nValidHitsTOBVsEta", "Number of Valid Hits in TOB Vs Eta",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");   
  nValidHitsTECVsEtaH_ = ibook.bookProfile("nValidHitsTECVsEta", "Number of Valid Hits in TEC Vs Eta",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsTIDVsEtaH_ = ibook.bookProfile("nValidHitsTIDVsEta", "Number of Valid Hits in TID Vs Eta",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");

  nValidHitsPixVsPhiH_ = ibook.bookProfile("nValidHitsPixVsPhi", "Number of Valid Hits in Pixel Vs Phi",
					   TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					   TrackEtaHistoPar.getParameter<double>("Xmin"),
					   TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsPixBVsPhiH_ = ibook.bookProfile("nValidHitsPixBVsPhi", "Number of Valid Hits in Pixel Barrel Vs Phi",
					    TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					    TrackEtaHistoPar.getParameter<double>("Xmin"),
					    TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsPixEVsPhiH_ = ibook.bookProfile("nValidHitsPixEVsPhi", "Number of Valid Hits in Pixel Endcap Vs Phi",
					    TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					    TrackEtaHistoPar.getParameter<double>("Xmin"),
					    TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsStripVsPhiH_ = ibook.bookProfile("nValidHitsStripVsPhi", "Number of Valid Hits in SiStrip Vs Phi",
					     TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					     TrackEtaHistoPar.getParameter<double>("Xmin"),
					     TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsTIBVsPhiH_ = ibook.bookProfile("nValidHitsTIBVsPhi", "Number of Valid Hits in TIB Vs Phi",
					   TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					   TrackEtaHistoPar.getParameter<double>("Xmin"),
					   TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsTOBVsPhiH_ = ibook.bookProfile("nValidHitsTOBVsPhi", "Number of Valid Hits in TOB Vs Phi",
					   TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					   TrackEtaHistoPar.getParameter<double>("Xmin"),
					   TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsTECVsPhiH_ = ibook.bookProfile("nValidHitsTECVsPhi", "Number of Valid Hits in TEC Vs Phi",
					   TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					   TrackEtaHistoPar.getParameter<double>("Xmin"),
					   TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nValidHitsTIDVsPhiH_ = ibook.bookProfile("nValidHitsTIDVsPhi", "Number of Valid Hits in TID Vs Phi",
					   TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					   TrackEtaHistoPar.getParameter<double>("Xmin"),
					   TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");

  nLostHitsPixVsEtaH_ = ibook.bookProfile("nLostHitsPixVsEta", "Number of Lost Hits in Pixel Vs Eta",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsPixBVsEtaH_ = ibook.bookProfile("nLostHitsPixBVsEta", "Number of Lost Hits in Pixel Barrel Vs Eta",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsPixEVsEtaH_ = ibook.bookProfile("nLostHitsPixEVsEta", "Number of Lost Hits in Pixel Endcap Vs Eta",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsStripVsEtaH_ = ibook.bookProfile("nLostHitsStripVsEta", "Number of Lost Hits in SiStrip Vs Eta",
				       TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				       TrackEtaHistoPar.getParameter<double>("Xmin"),
				       TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsTIBVsEtaH_ = ibook.bookProfile("nLostHitsTIBVsEta", "Number of Lost Hits in TIB Vs Eta",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsTOBVsEtaH_ = ibook.bookProfile("nLostHitsTOBVsEta", "Number of Lost Hits in TOB Vs Eta",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");   
  nLostHitsTECVsEtaH_ = ibook.bookProfile("nLostHitsTECVsEta", "Number of Lost Hits in TEC Vs Eta",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsTIDVsEtaH_ = ibook.bookProfile("nLostHitsTIDVsEta", "Number of Lost Hits in TID Vs Eta",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");

  nLostHitsPixVsPhiH_ = ibook.bookProfile("nLostHitsPixVsPhi", "Number of Lost Hits in Pixel Vs Phi",
					   TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					   TrackEtaHistoPar.getParameter<double>("Xmin"),
					   TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsPixBVsPhiH_ = ibook.bookProfile("nLostHitsPixBVsPhi", "Number of Lost Hits in Pixel Barrel Vs Phi",
					    TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					    TrackEtaHistoPar.getParameter<double>("Xmin"),
					    TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsPixEVsPhiH_ = ibook.bookProfile("nLostHitsPixEVsPhi", "Number of Lost Hits in Pixel Endcap Vs Phi",
					    TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					    TrackEtaHistoPar.getParameter<double>("Xmin"),
					    TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsStripVsPhiH_ = ibook.bookProfile("nLostHitsStripVsPhi", "Number of Lost Hits in SiStrip Vs Phi",
					     TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					     TrackEtaHistoPar.getParameter<double>("Xmin"),
					     TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsTIBVsPhiH_ = ibook.bookProfile("nLostHitsTIBVsPhi", "Number of Lost Hits in TIB Vs Phi",
					   TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					   TrackEtaHistoPar.getParameter<double>("Xmin"),
					   TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsTOBVsPhiH_ = ibook.bookProfile("nLostHitsTOBVsPhi", "Number of Lost Hits in TOB Vs Phi",
					   TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					   TrackEtaHistoPar.getParameter<double>("Xmin"),
					   TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsTECVsPhiH_ = ibook.bookProfile("nLostHitsTECVsPhi", "Number of Lost Hits in TEC Vs Phi",
					   TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					   TrackEtaHistoPar.getParameter<double>("Xmin"),
					   TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  nLostHitsTIDVsPhiH_ = ibook.bookProfile("nLostHitsTIDVsPhi", "Number of Lost Hits in TID Vs Phi",
					   TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
					   TrackEtaHistoPar.getParameter<double>("Xmin"),
					   TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");

  trackChi2oNDFVsEtaH_ = ibook.bookProfile("trackChi2oNDFVsEta", "chi2/ndof of Tracks Vs Eta",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");

  trackChi2oNDFVsPhiH_ = ibook.bookProfile("trackChi2oNDFVsPhi", "chi2/ndof of Tracks Vs Phi",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  
  trackChi2probVsEtaH_ = ibook.bookProfile("trackChi2probVsEta", "chi2 probability of Tracks Vs Eta",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");

  trackChi2probVsPhiH_ = ibook.bookProfile("trackChi2probVsPhi", "chi2 probability of Tracks Vs Phi",
				      TrackEtaHistoPar.getParameter<int32_t>("Xbins"),
				      TrackEtaHistoPar.getParameter<double>("Xmin"),
				      TrackEtaHistoPar.getParameter<double>("Xmax"),0.0,0.0,"g");
  
  


  // On and off-track cluster properties
  hOnTrkClusChargeThinH_ = ibook.book1D("hOnTrkClusChargeThin", "On-track Cluster Charge (Thin Sensor)", 100, 0, 1000);
  hOnTrkClusWidthThinH_ = ibook.book1D("hOnTrkClusWidthThin", "On-track Cluster Width (Thin Sensor)", 20, -0.5, 19.5);
  hOnTrkClusChargeThickH_ = ibook.book1D("hOnTrkClusChargeThick", "On-track Cluster Charge (Thick Sensor)", 100, 0, 1000);
  hOnTrkClusWidthThickH_ = ibook.book1D("hOnTrkClusWidthThick", "On-track Cluster Width (Thick Sensor)", 20, -0.5, 19.5);
  
  hOffTrkClusChargeThinH_ = ibook.book1D("hOffTrkClusChargeThin", "Off-track Cluster Charge (Thin Sensor)", 100, 0, 1000);
  hOffTrkClusWidthThinH_ = ibook.book1D("hOffTrkClusWidthThin", "Off-track Cluster Width (Thin Sensor)", 20, -0.5, 19.5);
  hOffTrkClusChargeThickH_ = ibook.book1D("hOffTrkClusChargeThick", "Off-track Cluster Charge (Thick Sensor)", 100, 0, 1000);
  hOffTrkClusWidthThickH_ = ibook.book1D("hOffTrkClusWidthThick", "Off-track Cluster Width (Thick Sensor)", 20, -0.5, 19.5);
}
void StandaloneTrackMonitor::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  if (verbose_) std::cout << "Begin StandaloneTrackMonitor" << std::endl;

  //std::cout << "DEBU 1" << std::endl;  
  // Get event setup (to get global transformation)                                  
  edm::ESHandle<TrackerGeometry> geomHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(geomHandle);
  const TrackerGeometry& tkGeom = (*geomHandle);
  
  //std::cout << "DEBU 2" << std::endl;  
  // Primary vertex collection
  edm::Handle<reco::VertexCollection> vertexColl;
  iEvent.getByToken(vertexToken_, vertexColl);
  if (!vertexColl.isValid()) {
    std::cerr << "Error! Failed to get reco::Vertex Collection, for " << vertexTag_ << std::endl;
    edm::LogError("DqmTrackStudy") << "Error! Failed to get reco::Vertex Collection, " << vertexTag_;
  }
  if (vertexColl->size() < 1) {
    std::cerr << "No good vertex in the event!!" << std::endl;
    return;
  }
  const reco::Vertex& pv = (*vertexColl)[0];

  //std::cout << "DEBU 3" << std::endl;  
  // Beam spot
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(bsToken_, beamSpot);
  if (!beamSpot.isValid()) std::cerr << "Beamspot for input tag: " << bsTag_ << " not found!!" << std::endl;

  //std::cout << "DEBU 4" << std::endl;  
  // Track collection
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackToken_, tracks);
  if (!tracks.isValid()) std::cerr << "TrackCollection for input tag: " << trackTag_ << " not found!!" << std::endl;

  //std::cout << "DEBU 5" << std::endl;  
  // Access PU information
  double wfac = 1.0;  // for data
  if (!iEvent.isRealData()) {
    edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(puSummaryToken_, PupInfo);
    
    if (verbose_) edm::LogInfo("StandaloneTrackMonitor") << "nPUColl = " << PupInfo->size();
    if (PupInfo.isValid()) {
      for (auto const& v : *PupInfo) {
	int bx = v.getBunchCrossing();
	if (bunchCrossingH_) bunchCrossingH_->Fill(bx);
	if (bx == 0) {
	  if (nPUH_) nPUH_->Fill(v.getPU_NumInteractions());
	  int ntrueInt = v.getTrueNumInteractions();
	  int nVertex = (vertexColl.isValid() ? vertexColl->size() : 0);
	  if (trueNIntH_) trueNIntH_->Fill(ntrueInt);
	  if (doPUCorrection_) {
	    if (nVertex > -1 && nVertex < int(vpu_.size())) wfac = vpu_.at(nVertex);
	    else wfac = 0.0;
	    //if (ntrueInt > -1 && ntrueInt < int(vpu_.size())) wfac = vpu_.at(ntrueInt);
	  }
	}
      }
    }
    else 
      std::cerr << "PUSummary for input tag: " << puSummaryTag_ << " not found!!" << std::endl;
  }
  if (verbose_) edm::LogInfo("StandaloneTrackMonitor") << "PU reweight factor = " << wfac;
  if (verbose_) std::cout << "PU scale factor" << wfac << std::endl;
  
  if (haveAllHistograms_) {
    int nvtx = (vertexColl.isValid() ? vertexColl->size() : 0);
    nVertexH_->Fill(nvtx, wfac);
    nVtxH_->Fill(nvtx);
  }
  
  //std::cout << "DEBU 7" << std::endl;  
  // Get MVA and quality mask collections
  std::vector<const MVACollection*> mvaCollections;
  std::vector<const QualityMaskCollection*> qualityMaskCollections;
  for (const auto& tokenTpl : mvaQualityTokens_) {
    edm::Handle<MVACollection> hmva;
    iEvent.getByToken(std::get<0>(tokenTpl), hmva);
    
    edm::Handle<QualityMaskCollection> hqual;
    iEvent.getByToken(std::get<1>(tokenTpl), hqual);
    
    if (hmva.isValid() && hqual.isValid()) {
      mvaCollections.push_back(hmva.product());
      qualityMaskCollections.push_back(hqual.product());
    }
    else {
      std::cerr << "Track MVA related handles not accessible!" << std::endl;
    }
  }
  
  //std::cout << "DEBU 8" << std::endl;  
  edm::Handle<edm::View<reco::Track> > htracks;
  iEvent.getByToken(mvaTrackToken_, htracks);
  if (htracks.isValid()) {
    edm::LogInfo("StandaloneTrackMonitor") << "Total # of Tracks: " << htracks->size();
    if (verbose_) edm::LogInfo("StandaloneTrackMonitor") <<"Total # of Tracks: " << htracks->size();
    
    int ntracks = htracks->size();
    auto nmva = mvaCollections.size();
    
    std::cout << "#mva Tracks: " << ntracks << std::endl;
    for (auto iTrack = 0; iTrack < ntracks; ++iTrack) {
      const auto& track = (*tracks)[iTrack];
      
      // Fill MVA1 histos with all tracks, MVA2 histos only with tracks
      // not selected by MVA1 etc
      bool selectedLoose = false;
      bool selectedHP = false;
      
      const auto pt = track.pt();
      const auto eta = track.eta();
      
      for (size_t iMVA = 0; iMVA < nmva; ++iMVA) {
	const auto mva = (*(mvaCollections[iMVA]))[iTrack];
	if (!selectedLoose) {
	  trackMVAs[iMVA]->Fill(mva);
	  trackMVAsVsPtProfile[iMVA]->Fill(pt, mva);
	  trackMVAsVsEtaProfile[iMVA]->Fill(eta, mva);
	}
	if (iMVA >= 1 && !selectedHP) {
	  trackMVAsHP[iMVA]->Fill(mva);
	  trackMVAsHPVsPtProfile[iMVA]->Fill(pt, mva);
	  trackMVAsHPVsEtaProfile[iMVA]->Fill(eta, mva);
	}
	
	const auto qual = (*(qualityMaskCollections)[iMVA])[iTrack];
	selectedLoose |= trackSelected(qual, reco::TrackBase::loose);
	selectedHP |= trackSelected(qual, reco::TrackBase::highPurity);
	
	if (selectedLoose && selectedHP)
	  break;
      }
    }
  }
  else {
    std::cerr << "Track collection for " << mvaTrackTag_ << " not found!" << std::endl;
  }
  //std::cout << "DEBU 9" << std::endl;  
  int ntracks = 0;
  if (tracks.isValid()) {
    edm::LogInfo("StandaloneTrackMonitor") << "Total # of Tracks: " << tracks->size();
    if (verbose_) edm::LogInfo("StandaloneTrackMonitor") <<"Total # of Tracks: " << tracks->size();
    reco::Track::TrackQuality quality = reco::Track::qualityByName(trackQuality_);
    if (verbose_) std::cout <<"Total # of Tracks: " << tracks->size() << std::endl;

    for (auto const& track : *tracks) {
      if (!track.quality(quality)) continue;
      ++ntracks;     
 
      double eta = track.eta();
      double theta = track.theta();
      double phi = track.phi();
      double pt = track.pt();

      const reco::HitPattern& hitp = track.hitPattern();
      double nValidTrackerHits = hitp.numberOfValidTrackerHits();
      double nValidPixelHits = hitp.numberOfValidPixelHits();
      double nValidPixelBHits = hitp.numberOfValidPixelBarrelHits();
      double nValidPixelEHits = hitp.numberOfValidPixelEndcapHits();
      double nValidStripHits = hitp.numberOfValidStripHits();
      double nValidTIBHits = hitp.numberOfValidStripTIBHits();
      double nValidTOBHits = hitp.numberOfValidStripTOBHits();
      double nValidTIDHits = hitp.numberOfValidStripTIDHits();
      double nValidTECHits = hitp.numberOfValidStripTECHits();

      int missingInnerHit = hitp.numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
      int missingOuterHit = hitp.numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS);

      nValidHitsVspTH_->Fill(pt, nValidTrackerHits);
      nValidHitsVsEtaH_->Fill(eta, nValidTrackerHits);
      nValidHitsVsCosThetaH_->Fill(std::cos(theta), nValidTrackerHits);
      nValidHitsVsPhiH_->Fill(phi, nValidTrackerHits);
      nValidHitsVsnVtxH_->Fill(vertexColl->size(), nValidTrackerHits);

      nValidHitsPixVsEtaH_->Fill(eta, nValidPixelHits);
      nValidHitsPixBVsEtaH_->Fill(eta, nValidPixelBHits);
      nValidHitsPixEVsEtaH_->Fill(eta, nValidPixelEHits);
      nValidHitsStripVsEtaH_->Fill(eta, nValidStripHits);
      nValidHitsTIBVsEtaH_->Fill(eta, nValidTIBHits);
      nValidHitsTOBVsEtaH_->Fill(eta, nValidTOBHits);
      nValidHitsTECVsEtaH_->Fill(eta, nValidTECHits);
      nValidHitsTIDVsEtaH_->Fill(eta, nValidTIDHits);

      nValidHitsPixVsPhiH_->Fill(phi, nValidPixelHits);
      nValidHitsPixBVsPhiH_->Fill(phi, nValidPixelBHits);
      nValidHitsPixEVsPhiH_->Fill(phi, nValidPixelEHits);
      nValidHitsStripVsPhiH_->Fill(phi, nValidStripHits);
      nValidHitsTIBVsPhiH_->Fill(phi, nValidTIBHits);
      nValidHitsTOBVsPhiH_->Fill(phi, nValidTOBHits);
      nValidHitsTECVsPhiH_->Fill(phi, nValidTECHits);
      nValidHitsTIDVsPhiH_->Fill(phi, nValidTIDHits);

      int nLostHits = track.numberOfLostHits();
      int nLostTrackerHits = hitp.numberOfLostTrackerHits(reco::HitPattern::TRACK_HITS);
      int nLostPixHits = hitp.numberOfLostPixelHits(reco::HitPattern::TRACK_HITS);
      int nLostPixBHits = hitp.numberOfLostPixelBarrelHits(reco::HitPattern::TRACK_HITS);
      int nLostPixEHits = hitp.numberOfLostPixelEndcapHits(reco::HitPattern::TRACK_HITS);
      int nLostStripHits = hitp.numberOfLostStripHits(reco::HitPattern::TRACK_HITS);
      int nLostStripTIBHits = hitp.numberOfLostStripTIBHits(reco::HitPattern::TRACK_HITS);
      int nLostStripTIDHits = hitp.numberOfLostStripTIDHits(reco::HitPattern::TRACK_HITS);
      int nLostStripTOBHits = hitp.numberOfLostStripTOBHits(reco::HitPattern::TRACK_HITS);
      int nLostStripTECHits = hitp.numberOfLostStripTECHits(reco::HitPattern::TRACK_HITS);
      
      nLostHitsVspTH_->Fill(pt, nLostTrackerHits);
      nLostHitsVsEtaH_->Fill(eta, nLostTrackerHits);
      nLostHitsVsCosThetaH_->Fill(std::cos(theta), nLostTrackerHits);
      nLostHitsVsPhiH_->Fill(phi, nLostTrackerHits);

      nLostHitsPixVsEtaH_->Fill(eta, nLostPixHits);
      nLostHitsPixBVsEtaH_->Fill(eta, nLostPixBHits);
      nLostHitsPixEVsEtaH_->Fill(eta, nLostPixEHits);
      nLostHitsStripVsEtaH_->Fill(eta, nLostStripHits);
      nLostHitsTIBVsEtaH_->Fill(eta, nLostStripTIBHits);
      nLostHitsTOBVsEtaH_->Fill(eta, nLostStripTOBHits);
      nLostHitsTECVsEtaH_->Fill(eta, nLostStripTECHits);
      nLostHitsTIDVsEtaH_->Fill(eta, nLostStripTIDHits);

      nLostHitsPixVsPhiH_->Fill(phi, nLostPixHits);
      nLostHitsPixBVsPhiH_->Fill(phi, nLostPixBHits);
      nLostHitsPixEVsPhiH_->Fill(phi, nLostPixEHits);
      nLostHitsStripVsPhiH_->Fill(phi, nLostStripHits);
      nLostHitsTIBVsPhiH_->Fill(phi, nLostStripTIBHits);
      nLostHitsTOBVsPhiH_->Fill(phi, nLostStripTOBHits);
      nLostHitsTECVsPhiH_->Fill(phi, nLostStripTECHits);
      nLostHitsTIDVsPhiH_->Fill(phi, nLostStripTIDHits);


      if (abs(eta) <= 1.4) {
	nMissingInnerHitBH_->Fill(missingInnerHit);
	nMissingOuterHitBH_->Fill(missingOuterHit);
      }
      else {
	nMissingInnerHitEH_->Fill(missingInnerHit);
        nMissingOuterHitEH_->Fill(missingOuterHit);
      }

      // Attention !! Residual is only for RECO dataset
      //std::cout << "DEBU 10" << std::endl;  
      
      double residualXPB = 0, residualXPE = 0, residualXTIB = 0, residualXTOB = 0, residualXTEC = 0, residualXTID = 0;
      double residualYPB = 0, residualYPE = 0, residualYTIB = 0, residualYTOB = 0, residualYTEC = 0, residualYTID = 0; 
          reco::TrackResiduals residuals = track.residuals();
      int i = 0;
      for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it,++i) {
        const TrackingRecHit& hitn = (**it);
	if (hitn.isValid()) {
	  int subdet = hitn.geographicalId().subdetId();

	  if      (subdet == PixelSubdetector::PixelBarrel) {
	    residualXPB  = residuals.residualX(i);
	    residualYPB  = residuals.residualY(i);
	  }
	  else if (subdet == PixelSubdetector::PixelEndcap) {
	    residualXPE  = residuals.residualX(i);
	    residualYPE  = residuals.residualY(i);
	  }
	  else if (subdet == StripSubdetector::TIB) {
	    residualXTIB = residuals.residualX(i);
	    residualYTIB = residuals.residualY(i);
	  }
	  else if (subdet == StripSubdetector::TOB) {
	    residualXTOB = residuals.residualX(i);
	    residualYTOB = residuals.residualY(i);
	  }
	  else if (subdet == StripSubdetector::TEC) {
	    residualXTEC = residuals.residualX(i);
	    residualYTEC = residuals.residualY(i);
	  }
	  else if (subdet == StripSubdetector::TID) {
	    residualXTID = residuals.residualX(i);
	    residualYTID = residuals.residualY(i);
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

      //std::cout << "DEBU 11" << std::endl;  
      for (int i = 0; i < hitp.numberOfAllHits(reco::HitPattern::TRACK_HITS); i++) {
	uint32_t hit = hitp.getHitPattern(reco::HitPattern::TRACK_HITS, i);
	if (hitp.missingHitFilter(hit)) {
	  double losthitBylayer = -1.0;
	  double losthitBylayerPix = -1.0;
	  double losthitBylayerStrip = -1.0;
	  int layer = hitp.getLayer(hit);
	  //int side = hitp.getSide(hit);
	  if (hitp.pixelBarrelHitFilter(hit)) {
	    losthitBylayer = layer;
	    losthitBylayerPix = layer;	    
	  }
	  else if (hitp.pixelEndcapHitFilter(hit)) {
	    //losthitBylayer = (side == 0) ? layer+3 : layer+5;
	    losthitBylayer = layer+4;
	    losthitBylayerPix = layer+4;
	  }
	  else if (hitp.stripTIBHitFilter(hit)) {
	    losthitBylayer = layer + 7;
	    losthitBylayerStrip = layer;
	  }
	  else if (hitp.stripTIDHitFilter(hit)) {
	    //losthitBylayer = (side == 0) ? layer+11 : layer+14;
	    losthitBylayer = layer+11;
	    losthitBylayerStrip = layer+4;
	  }
	  else if (hitp.stripTOBHitFilter(hit)) {
	    losthitBylayer = layer + 14;
	    losthitBylayerStrip = layer + 7;
	  }
	  else if (hitp.stripTECHitFilter(hit)) {
	    //losthitBylayer = (side == 0) ? layer+23 : layer+32;
	    losthitBylayer = layer+20;
	    losthitBylayerStrip = layer+13;
	  }
	  if (losthitBylayer > -1) nLostHitByLayerH_->Fill(losthitBylayer, wfac);
	  if (losthitBylayerPix > -1) nLostHitByLayerPixH_->Fill(losthitBylayerPix, wfac);
	  if (losthitBylayerStrip > -1) nLostHitByLayerStripH_->Fill(losthitBylayerStrip, wfac);
	}
      }
      
      //std::cout << "DEBU 12" << std::endl;  
      if (haveAllHistograms_) {
	double etaError = track.etaError();
	double thetaError = track.thetaError();
	double phiError = track.phiError();
	double p = track.p();
	double ptError = track.ptError();
	double qoverp = track.qoverp();
	double qoverpError = track.qoverpError();
	double charge = track.charge();
	
	double dxy = track.dxy(beamSpot->position());
	double dxyError = track.dxyError();
	double dz = track.dz(beamSpot->position());
	double dzError = track.dzError();
	
	double trkd0 = track.d0();     
	double chi2 = track.chi2();
	double ndof = track.ndof();
	double chi2prob = TMath::Prob(track.chi2(),(int)track.ndof());
	double chi2oNDF = track.normalizedChi2();
	double vx = track.vx();
	double vy = track.vy();
	double vz = track.vz();
	unsigned int track_algo = track.algo();
	unsigned int track_origalgo = track.originalAlgo();
 	double distanceOfClosestApproachToPV = track.dxy(pv.position());
	double xPointOfClosestApproachwrtPV = track.vx()-pv.position().x();
	double yPointOfClosestApproachwrtPV = track.vy()-pv.position().y();
	double positionZ0 = track.dz(pv.position());
	
	edm::ESHandle<TransientTrackBuilder> theB;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
	reco::TransientTrack transTrack = theB->build(track);
	
	double sip3dToPV = 0, sip2dToPV = 0;
	GlobalVector dir(track.px(), track.py(), track.pz());
	std::pair<bool, Measurement1D> ip3d = IPTools::signedImpactParameter3D(transTrack, dir, pv);
	std::pair<bool, Measurement1D> ip2d = IPTools::signedTransverseImpactParameter(transTrack, dir, pv);
	if(ip3d.first) sip3dToPV = ip3d.second.value()/ip3d.second.error();
	if(ip2d.first) sip2dToPV = ip2d.second.value()/ip2d.second.error();
        double sipDxyToPV = track.dxy(pv.position())/track.dxyError();
	double sipDzToPV = track.dz(pv.position())/track.dzError();				             
	
	// Fill the histograms
	trackEtaH_->Fill(eta, wfac);
	trackEtaerrH_->Fill(etaError, wfac);
	trackCosThetaH_->Fill(std::cos(theta), wfac);
	trackThetaerrH_->Fill(thetaError, wfac);
	trackPhiH_->Fill(phi, wfac);
	trackPhierrH_->Fill(phiError, wfac);
	trackPH_->Fill(p, wfac);
	trackPtH_->Fill(pt, wfac);
	if (pt <= 2) trackPtUpto2GeVH_->Fill(pt, wfac);
	if (pt >= 10) trackPtOver10GeVH_->Fill(pt, wfac);
	trackPterrH_->Fill(ptError, wfac);
	trackqOverpH_->Fill(qoverp, wfac);
	trackqOverperrH_->Fill(qoverpError, wfac);
	trackChargeH_->Fill(charge, wfac);
	trackChi2H_->Fill(chi2, wfac);
	trackChi2ProbH_->Fill(chi2prob, wfac);
	trackChi2oNDFH_->Fill(chi2oNDF, wfac);
	trackd0H_->Fill(trkd0, wfac);
        tracknDOFH_->Fill(ndof, wfac);
	trackChi2bynDOFH_->Fill(chi2/ndof, wfac);
	trackalgoH_->Fill(track_algo, wfac);
	trackorigalgoH_->Fill(track_origalgo, wfac);
	trackChi2oNDFVsEtaH_->Fill(eta, chi2oNDF);
	trackChi2oNDFVsPhiH_->Fill(phi, chi2oNDF);
	trackChi2probVsEtaH_->Fill(eta, chi2prob);
	trackChi2probVsPhiH_->Fill(phi, chi2prob);
	
	nlostHitsH_->Fill(nLostHits, wfac);
	nlostTrackerHitsH_->Fill(nLostTrackerHits, wfac);
	
	beamSpotXYposH_->Fill(dxy, wfac);
	beamSpotXYposerrH_->Fill(dxyError, wfac);
	beamSpotZposH_->Fill(dz, wfac);
	beamSpotZposerrH_->Fill(dzError, wfac);
	
	vertexXposH_->Fill(vx, wfac);
	vertexYposH_->Fill(vy, wfac);
	vertexZposH_->Fill(vz, wfac);
	
        DistanceOfClosestApproachToPVH_->Fill(distanceOfClosestApproachToPV);
	DistanceOfClosestApproachToPVVsPhiH_->Fill(phi, distanceOfClosestApproachToPV);
	xPointOfClosestApproachVsZ0wrtPVH_->Fill(positionZ0, xPointOfClosestApproachwrtPV);
	yPointOfClosestApproachVsZ0wrtPVH_->Fill(positionZ0, yPointOfClosestApproachwrtPV);
	
	sip3dToPVH_->Fill(sip3dToPV);
	sip2dToPVH_->Fill(sip2dToPV);
	sipDxyToPVH_->Fill(sipDxyToPV);
	sipDzToPVH_->Fill(sipDzToPV);
       	
        double trackerLayersWithMeasurement = hitp.trackerLayersWithMeasurement();
	double pixelLayersWithMeasurement = hitp.pixelLayersWithMeasurement();
	double pixelBLayersWithMeasurement = hitp.pixelBarrelLayersWithMeasurement();
	double pixelELayersWithMeasurement = hitp.pixelEndcapLayersWithMeasurement();
	double stripLayersWithMeasurement = hitp.stripLayersWithMeasurement();
	double stripTIBLayersWithMeasurement = hitp.stripTIBLayersWithMeasurement();
	double stripTOBLayersWithMeasurement = hitp.stripTOBLayersWithMeasurement();
	double stripTIDLayersWithMeasurement = hitp.stripTIDLayersWithMeasurement();
	double stripTECLayersWithMeasurement = hitp.stripTECLayersWithMeasurement();
	
	trkLayerwithMeasurementH_->Fill(trackerLayersWithMeasurement, wfac);
	pixelLayerwithMeasurementH_->Fill(pixelLayersWithMeasurement, wfac);
	pixelBLayerwithMeasurementH_->Fill(pixelBLayersWithMeasurement, wfac);
	pixelELayerwithMeasurementH_->Fill(pixelELayersWithMeasurement, wfac);
	stripLayerwithMeasurementH_->Fill(stripLayersWithMeasurement, wfac);
	stripTIBLayerwithMeasurementH_->Fill(stripTIBLayersWithMeasurement, wfac);
	stripTOBLayerwithMeasurementH_->Fill(stripTOBLayersWithMeasurement, wfac);
	stripTIDLayerwithMeasurementH_->Fill(stripTIDLayersWithMeasurement, wfac);
	stripTECLayerwithMeasurementH_->Fill(stripTECLayersWithMeasurement, wfac);
	
	nvalidTrackerHitsH_->Fill(nValidTrackerHits, wfac);
	nvalidPixelHitsH_->Fill(nValidPixelHits, wfac);
	nvalidPixelBHitsH_->Fill(nValidPixelBHits, wfac);
	nvalidPixelEHitsH_->Fill(nValidPixelEHits, wfac);
	nvalidStripHitsH_->Fill(nValidStripHits, wfac);
	nvalidTIBHitsH_->Fill(nValidTIBHits, wfac);
	nvalidTOBHitsH_->Fill(nValidTOBHits, wfac);
	nvalidTIDHitsH_->Fill(nValidTIDHits, wfac);
	nvalidTECHitsH_->Fill(nValidTECHits, wfac);
	
	nlostTrackerHitsH_->Fill(nLostTrackerHits, wfac);
	nlostPixelHitsH_->Fill(nLostPixHits, wfac);
	nlostPixelBHitsH_->Fill(nLostPixBHits, wfac);
	nlostPixelEHitsH_->Fill(nLostPixEHits, wfac);
	nlostStripHitsH_->Fill(nLostStripHits, wfac);
	nlostTIBHitsH_->Fill(nLostStripTIBHits, wfac);
	nlostTOBHitsH_->Fill(nLostStripTOBHits, wfac);
	nlostTIDHitsH_->Fill(nLostStripTIDHits, wfac);
	nlostTECHitsH_->Fill(nLostStripTECHits, wfac);
      }

      //std::cout << "DEBU 13" << std::endl;  

      // Attention !! The following loop is only for RECO dataset
      //*/
      int nStripTIBS = 0, nStripTOBS = 0, nStripTECS = 0, nStripTIDS = 0;
      int nStripTIBD = 0, nStripTOBD = 0, nStripTECD = 0, nStripTIDD = 0;
      
      for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it) {
	//std::cout << "DEBU 13_0_0" << std::endl;
	const TrackingRecHit& hit = (**it);
	if (hit.isValid()) {
	  if (hit.geographicalId().det() == DetId::Tracker) {
	    //std::cout << "DEBU 13_0" << std::endl;
	    int subdetId = hit.geographicalId().subdetId();
            //std::cout << "DEBU 13_1" << std::endl;
            // Find on-track clusters
            processHit(hit, iSetup, tkGeom, wfac);
	    //std::cout << "DEBU 13_2" << std::endl;
	    
	    const DetId detId(hit.geographicalId());
	    const SiStripDetId stripId(detId);
	    if (0) std::cout << "Hit Dimension: " << hit.dimension()
			     << ", isGlued: " << stripId.glued()
			     << ", isStereo: " << stripId.stereo()
			     << std::endl;
	    //std::cout << "DEBU 13_3" << std::endl;
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
	    //  std::cout << "DEBU 13_4" << std::endl;
	  }
	}
      }   
      
      nHitsTIBSVsEtaH_->Fill(eta, nStripTIBS);
      nHitsTOBSVsEtaH_->Fill(eta, nStripTOBS);
      nHitsTECSVsEtaH_->Fill(eta, nStripTECS);
      nHitsTIDSVsEtaH_->Fill(eta, nStripTIDS);
      nHitsStripSVsEtaH_->Fill(eta, nStripTIBS+nStripTOBS+nStripTECS+nStripTIDS);
      //std::cout << "DEBU 13_5" << std::endl;
      nHitsTIBDVsEtaH_->Fill(eta, nStripTIBD);
      nHitsTOBDVsEtaH_->Fill(eta, nStripTOBD);
      nHitsTECDVsEtaH_->Fill(eta, nStripTECD);
      nHitsTIDDVsEtaH_->Fill(eta, nStripTIDD);
      nHitsStripDVsEtaH_->Fill(eta, nStripTIBD+nStripTOBD+nStripTECD+nStripTIDD);
      //*/
      //      std::cout << "DEBU 14" << std::endl;  
    }
  }
  else {
    edm::LogError("StandaloneTrackMonitor") << "Error! Failed to get reco::Track collection, " << trackTag_;
  }
  if (haveAllHistograms_) nTracksH_->Fill(ntracks, wfac);
  
  // off track cluster properties
  processClusters(iEvent, iSetup, tkGeom, wfac);
  
  if (verbose_) std::cout << "Ends StandaloneTrackMonitor successfully" << std::endl;
}
void StandaloneTrackMonitor::processClusters(edm::Event const& iEvent, edm::EventSetup const& iSetup, const TrackerGeometry& tkGeom, double wfac)
{
  // SiStripClusters
  edm::Handle<edmNew::DetSetVector<SiStripCluster> > clusterHandle;
  iEvent.getByToken(clusterToken_, clusterHandle);

  if (clusterHandle.isValid()) {
    // Loop on Dets
    for (edmNew::DetSetVector<SiStripCluster>::const_iterator dsvit  = clusterHandle->begin(); 
                                                              dsvit != clusterHandle->end();
	                                                    ++dsvit)
    {
      uint32_t detId = dsvit->id();
      std::map<uint32_t, std::set<const SiStripCluster*> >::iterator jt = clusterMap_.find(detId);
      bool detid_found = (jt != clusterMap_.end()) ? true : false;

      // Loop on Clusters
      for (edmNew::DetSet<SiStripCluster>::const_iterator clusit  = dsvit->begin(); 
                                                          clusit != dsvit->end();	  
                                                        ++clusit)
      {
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
  else {
    edm::LogError("StandaloneTrackMonitor") << "ClusterCollection " << clusterTag_ << " not valid!!" << std::endl;
  }
}
void StandaloneTrackMonitor::processHit(const TrackingRecHit& recHit, edm::EventSetup const& iSetup, const TrackerGeometry& tkGeom, double wfac)
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
void StandaloneTrackMonitor::addClusterToMap(uint32_t detid, const SiStripCluster* cluster) {
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
void StandaloneTrackMonitor::endLuminosityBlock(edm::LuminosityBlock const& lumiBlock, edm::EventSetup const& eSetup){
}
// Define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(StandaloneTrackMonitor);
