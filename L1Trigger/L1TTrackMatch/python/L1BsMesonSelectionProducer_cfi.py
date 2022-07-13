import FWCore.ParameterSet.Config as cms

L1BsMesonSelectionProducer = cms.EDProducer('L1BsMesonSelectionProducer',
  l1PhiCandsInputTag = cms.InputTag("L1PhiMesonSelectionProducer","Level1TTPhiMesonSelected"),
  outputCollectionName = cms.string("Level1TTBsMesonSelected"),
  cutSet = cms.PSet(
                   dRmax = cms.double(1.0),
                   dRmin = cms.double(0.2),
                   dxymax = cms.double(1.0),
                   dzmax = cms.double(1.0),
                   phipairMmin = cms.double(5.29),
                   phipairMmax = cms.double(5.48)
                   ),
  debug = cms.int32(0)
)

#L1PhiMesonSelectionProducerExtended = L1PhiMesonSelectionProducer.clone(
  #l1TracksInputTag = cms.InputTag("L1TrackSelectionProducerExtended", "Level1TTTracksExtendedSelected"),
  #l1TracksInputTag = cms.InputTag("L1TrackNullSelectionProducerExtended", "Level1TTTracksExtendedNullSelected"),
  #outputCollectionName = "Level1TTKaonTracksExtendedSelected",
#)

