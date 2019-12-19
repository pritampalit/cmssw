import FWCore.ParameterSet.Config as cms
#from DQM.TrackingMonitorSource.StandaloneTrackMonitor_cfi import *
from DQM.TrackingMonitorSource.RecoTrackMonitor_cfi import *
# from DQM.TrackingMonitor.V0Monitor_cfi import *

# Primary Vertex Selector
selectedPrimaryVertices = cms.EDFilter("VertexSelector",
                                       src = cms.InputTag('offlinePrimaryVertices'),
                                       cut = cms.string("!isFake && ndof >= 4 && abs(z) < 24 && abs(position.Rho) < 2.0"),
                                       filter = cms.bool(True)
                                       )
# Track Selector
selectedTracks = cms.EDFilter("TrackSelector",
                              src = cms.InputTag('generalTracks'),
                              cut = cms.string("pt > 1.0"),
                              filter = cms.bool(True)
                              )
# HLT path selector
hltPathFilter = cms.EDFilter("HLTPathSelector",
                             processName = cms.string("HLT"),
                             verbose = cms.untracked.bool(False),
                             hltPathsOfInterest = cms.vstring("HLT_ZeroBias_v"),
                             triggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                             triggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT")
                             )

# Added module for V0Monitoring for Ks only
# KshortMonitor = v0Monitor.clone()
# KshortMonitor.FolderName = cms.string("Tracking/V0Monitoring/Ks")
# KshortMonitor.v0         = cms.InputTag('generalV0Candidates:Kshort')
# KshortMonitor.histoPSet.massPSet = cms.PSet(
#   nbins = cms.int32 ( 100 ),
#   xmin  = cms.double( 0.400),
#   xmax  = cms.double( 0.600),
# )

# For MinBias
recoTrackMonitorMC = recoTrackMonitor.clone(
    puScaleFactorFile = cms.untracked.string("PileupScaleFactor_316060_wrt_nVertex_ZeroBias.root"),
    doPUCorrection    = cms.untracked.bool(True),
    isMC              = cms.untracked.bool(True)
    )
recoValidationMinbias = cms.Sequence(
    hltPathFilter
    * selectedTracks
    * selectedPrimaryVertices
    * recoTrackMonitor	 )
recoValidationMinbiasMC = cms.Sequence(
    hltPathFilter
    * selectedTracks
    * selectedPrimaryVertices
    * recoTrackMonitorMC )
