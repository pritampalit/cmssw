import FWCore.ParameterSet.Config as cms

L1PhiMesonSelectionEmulationProducer = cms.EDProducer('L1PhiMesonSelectionEmulationProducer',
  l1PosKaonTracksInputTag = cms.InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedEmulationPositivecharge"),
  l1NegKaonTracksInputTag = cms.InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedEmulationNegativecharge"),
  outputCollectionName = cms.string("Level1TTPhiMesonSelectedEmulation"),
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

