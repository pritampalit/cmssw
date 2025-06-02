import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.jetsAK4_CHS_cff import *
from PhysicsTools.NanoAOD.jetsAK8_cff import *
from PhysicsTools.NanoAOD.jetMC_cff import *
from PhysicsTools.NanoAOD.jetConstituents_cff import *
from PhysicsTools.NanoAOD.muons_cff import *
from PhysicsTools.NanoAOD.taus_cff import *
from PhysicsTools.NanoAOD.boostedTaus_cff import *
from PhysicsTools.NanoAOD.electrons_cff import *
from PhysicsTools.NanoAOD.lowPtElectrons_cff import *
from PhysicsTools.NanoAOD.photons_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.extraflags_cff import *
from PhysicsTools.NanoAOD.ttbarCategorization_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.particlelevel_cff import *
from PhysicsTools.NanoAOD.genWeightsTable_cfi import *
from PhysicsTools.NanoAOD.tauSpinnerTable_cfi import *
from PhysicsTools.NanoAOD.genVertex_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.met_cff import *
from PhysicsTools.NanoAOD.triggerObjects_cff import *
from PhysicsTools.NanoAOD.isotracks_cff import *
from PhysicsTools.NanoAOD.protons_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.NanoAOD.fsrPhotons_cff import *
from PhysicsTools.NanoAOD.softActivity_cff import *


nanoMetadataCHS = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)

linkedObjectsCHS = cms.EDProducer("PATObjectCrossLinker",
   jets=cms.InputTag("finalJets"),
   muons=cms.InputTag("finalMuons"),
   electrons=cms.InputTag("finalElectrons"),
   lowPtElectrons=cms.InputTag("finalLowPtElectrons"),
   taus=cms.InputTag("finalTaus"),
   boostedTaus=cms.InputTag("finalBoostedTaus"),
   photons=cms.InputTag("finalPhotons"),
   vertices=cms.InputTag("slimmedSecondaryVertices")
)
'''
from PhysicsTools.NanoAOD.lhcInfoProducer_cfi import lhcInfoProducer
lhcInfoTable = lhcInfoProducer.clone()
(~run3_common).toModify(
    lhcInfoTable, useNewLHCInfo=False
)
'''


#jetTableCHS.variables.pt.precision=10

if hasattr(jetTable.variables, 'chFPV0EF'):
    del jetTable.variables.chFPV0EF

if hasattr(jetTable.variables, 'hfadjacentEtaStripsSize'):
    del jetTable.variables.hfadjacentEtaStripsSize

if hasattr(jetTable.variables, 'hfcentralEtaStripSize'):
    del jetTable.variables.hfcentralEtaStripSize

if hasattr(jetTable.variables, 'hfsigmaEtaEta'):
    del jetTable.variables.hfsigmaEtaEta

if hasattr(jetTable.variables, 'hfsigmaPhiPhi'):
    del jetTable.variables.hfsigmaPhiPhi

if hasattr(jetTable.variables, 'puId'):
    del jetTable.variables.puId

if hasattr(jetTable.variables, 'puIdDisc'):
    del jetTable.variables.puIdDisc
    
if hasattr(jetTable.variables, 'qgl'):
    del jetTable.variables.qgl
    
jetTable.externalVariables = cms.PSet() 


corrT1METJetTableCHS = corrT1METJetTable.clone(
    name = cms.string("CorrT1METJetCHS"),
    )
jetForMETTaskCHS =  cms.Task(basicJetsForMetForT1METNano,corrT1METJetTableCHS)

jetTableCHS = jetTable.clone(
    src = cms.InputTag("linkedObjectsCHS","jets"),
    name = cms.string("JeCHSt"),
    )


    
    
nanoTableTaskCommonCHS = cms.Task(
    cms.Task(nanoMetadataCHS),
    jetTask, jetForMETTaskCHS,
    #jetAK8Task, jetConstituentsTask,
    #extraFlagsProducersTask, muonTask, tauTask, boostedTauTask,
    #electronTask , lowPtElectronTask, photonTask,
    #vertexTask, isoTrackTask, jetAK8LepTask,  # must be after all the leptons
    #softActivityTask,
    cms.Task(linkedObjectsCHS),
    ##//jetTablesTask
    cms.Task(jetTableCHS)
    ##cms.Task(jetTableCHS)
    #jetAK8TablesTask, jetConstituentsTablesTask,
    #muonTablesTask, fsrTablesTask, tauTablesTask, boostedTauTablesTask,
    #electronTablesTask, lowPtElectronTablesTask, photonTablesTask,
    #globalTablesTask, vertexTablesTask, metTablesTask, extraFlagsTableTask,
    #isoTrackTablesTask,softActivityTablesTask
)

nanoSequenceCommonCHS = cms.Sequence(nanoTableTaskCommonCHS)

'''
nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectTablesTask)
nanoSequenceOnlyData = cms.Sequence(cms.Sequence(protonTablesTask) + lhcInfoTable)

nanoSequenceCHS = cms.Sequence(nanoSequenceCommonCHS + nanoSequenceOnlyData + nanoSequenceOnlyFullSim)

nanoTableTaskFS = cms.Task(
    genParticleTask, particleLevelTask, jetMCTask, muonMCTask, electronMCTask, lowPtElectronMCTask, photonMCTask,
    tauMCTask, boostedTauMCTask,
    metMCTable, ttbarCatMCProducersTask, globalTablesMCTask, ttbarCategoryTableTask,
    genWeightsTableTask, genVertexTablesTask, genParticleTablesTask, genProtonTablesTask, particleLevelTablesTask, tauSpinnerTableTask
)

nanoSequenceFS_CHS = cms.Sequence(nanoSequenceCommonCHS + cms.Sequence(nanoTableTaskFS))

# GenVertex only stored in newer MiniAOD
nanoSequenceMC_CHS = nanoSequenceFS_CHS.copy()
nanoSequenceMC_CHS.insert(nanoSequenceFS_CHS.index(nanoSequenceCommonCHS)+1,nanoSequenceOnlyFullSim)
'''
