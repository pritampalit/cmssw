import FWCore.ParameterSet.Config as cms

L1BsMesonSelectionEmulationProducer = cms.EDProducer('L1BsMesonSelectionEmulationProducer',
  l1PhiMesonWordInputTag = cms.InputTag("L1PhiMesonSelectionEmulationProducer","Level1TTPhiMesonSelectedEmulation"),
  outputCollectionName = cms.string("Level1TTBsMesonSelectedEmulation"),
  cutSet = cms.PSet(
                   dRmax = cms.double(0.12),
                   #dxymax = cms.double(1.0),
                   #dzmax = cms.double(1.0),
                   tkpairMmin = cms.double(1.0),
                   tkpairMmax = cms.double(1.03)
                   ),
  debug = cms.int32(0)
)

#L1PhiMesonSelectionProducerExtended = L1PhiMesonSelectionProducer.clone(
  #l1TracksInputTag = cms.InputTag("L1TrackSelectionProducerExtended", "Level1TTTracksExtendedSelected"),
  #l1TracksInputTag = cms.InputTag("L1TrackNullSelectionProducerExtended", "Level1TTTracksExtendedNullSelected"),
  #outputCollectionName = "Level1TTKaonTracksExtendedSelected",
#)

