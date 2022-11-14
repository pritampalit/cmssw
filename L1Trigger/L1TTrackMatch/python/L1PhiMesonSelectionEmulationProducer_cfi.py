import FWCore.ParameterSet.Config as cms
#from L1Trigger.VertexFinder.VertexProducer_cff import VertexProducer

L1PhiMesonSelectionEmulationProducer = cms.EDProducer('L1PhiMesonSelectionEmulationProducer',
  l1PosKaonTracksInputTag = cms.InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedEmulationPositivecharge"),
  l1NegKaonTracksInputTag = cms.InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedEmulationNegativecharge"),
#L1VertexInputTag = cms.InputTag("VertexProducer", VertexProducer.l1VertexCollectionName.value()),                                                      
  outputCollectionName = cms.string("Level1TTPhiMesonSelectedEmulation"),
  cutSet = cms.PSet(
                   dRmax = cms.double(0.12),
                   #dxymax = cms.double(1.0),
                   #dzmax = cms.double(1.0),
                   tkpairMmin = cms.double(1.0),
                   tkpairMmax = cms.double(1.03)
                   ),
  debug = cms.int32(0)
  #useGTTinput  = cms.bool( False )
)

#L1PhiMesonSelectionProducerExtended = L1PhiMesonSelectionProducer.clone(
  #l1TracksInputTag = cms.InputTag("L1TrackSelectionProducerExtended", "Level1TTTracksExtendedSelected"),
  #l1TracksInputTag = cms.InputTag("L1TrackNullSelectionProducerExtended", "Level1TTTracksExtendedNullSelected"),
  #outputCollectionName = "Level1TTKaonTracksExtendedSelected",
#)

