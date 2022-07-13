############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

############################################################
# edit options here
############################################################
L1TRK_INST ="MyL1TrackJets" ### if not in input DIGRAW then we make them in the above step
process = cms.Process(L1TRK_INST)

#L1TRKALGO = 'HYBRID'  #baseline, 4par fit
# L1TRKALGO = 'HYBRID_DISPLACED'  #extended, 5par fit
L1TRKALGO = 'HYBRID_PROMPTANDDISP'

DISPLACED = ''

############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '123X_mcRun4_realistic_v3', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.INFO.limit = cms.untracked.int32(0) # default: 0

############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

readFiles = cms.untracked.vstring(
"file:/afs/cern.ch/work/p/ppalit/public/HLSL1Trigger/EvtGenforBstoPhiPhi/test_benjamin_23434/step3.root",
###"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/c6df2819-ed05-4b98-8f92-81b7d1b1092e.root",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/3f476d95-1ef7-4be6-977b-6bcd1a7c5678.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/68d651da-4cb7-4bf4-b002-66aecc57a2bc.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/db0e0ce2-4c5a-4988-9dbd-52066e40b9d2.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/257a9712-0a96-47b7-897e-f5d980605e46.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/bee31399-8559-4243-b539-cae1ea897def.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/24629540-2377-4168-9ae5-518ddd4c43a9.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/e31ba8f0-332a-4a1a-8bc0-91a12a5fe3db.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/17902198-4db6-4fcc-9e8c-787991b4db32.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/4d433f3f-69e3-4912-bcdf-98a09d6d9e3d.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/59978fc3-0c28-455d-ac6a-a5a3bedee7fb.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/e4c534e1-31af-4bae-9914-0d4280e669fd.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/0c3bd223-4ad6-404d-b01e-761999abe00a.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/905c35d2-b35a-42c2-8a7f-a9979e3b8977.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/9a6ee436-2599-4f71-b310-9b9a2f929256.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/21123900-2086-4017-8167-bee36c4a1863.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/c693a68e-2af8-4894-a4c5-7aaee4778081.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/a156163a-e8a1-4584-9b78-8a692360b291.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/e9c823bc-75e4-474a-b49d-5f543df628e5.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/26208c7e-2b97-4d0c-aa02-579c09633b7b.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/2c646b27-4d06-4074-91c2-07991c0c17ed.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/abaa055c-de3a-4bea-8abf-ed16a4e359d1.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/075bb084-3865-4267-acd5-3ef4040648b0.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/91b02f36-794c-45dd-af73-49b53a64beff.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/696d9175-baad-42b3-b1f9-0b607526fe9c.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/3c134007-e362-45ee-b9b2-1a1a038061e2.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/f72f655d-5419-41ee-a1fb-c967d2a81e65.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/9f6ef9ae-2040-4746-a349-160a808f7ca0.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/f1e0a69f-73f6-4990-8671-ac3f00bc3ea8.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/f4bc045e-0eb9-4ba3-86a2-1cbe84021692.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/a5a8a4ca-9710-4711-9ca8-3b59c05e608e.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/86475540-5db7-4711-9012-b47538333232.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/c40a607f-7053-4600-8c6b-42302b45bff1.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/e091ed70-6a36-4893-b74b-845ff46a5d67.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/22aa6854-dbae-463d-846c-59f50143f6b7.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/82f8a879-b693-4615-8431-a046cb90659e.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/d345a34a-1e4f-4caf-9886-02f62b4b30cc.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/b50f9a35-25e1-4b36-9c04-56bbcfa43723.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/f1f02a41-73a3-4275-b02b-c28f2f78c0b3.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/31204f72-0194-4c33-bf5e-29ed0393ca97.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/2a2e4eaa-2a8d-4744-ab18-0c64cf39a287.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/c3f9e0c7-8ee7-4d93-b4ee-34a4f5c9b55f.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/06290813-320e-4dbb-8028-513013b77be3.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/7eef8b7f-6223-4b81-8a7b-693b4782b1bc.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/f0fdc961-deaa-4ba5-857a-706f90654a65.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/693d4240-af4b-420d-bc75-a6c1e15509b5.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/508674b5-880b-4e50-9bff-e5377fbeb4df.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/c7112d2d-0554-4607-9a1c-af5c1c3a058d.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/fe282d88-8ed6-4444-b5b7-374835486b57.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/212546db-1579-4f0e-a475-3549587de0c3.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/42560f29-3410-47df-b06e-7bf830d7e7d3.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/959077d7-283e-4a84-a604-65ca6061b172.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/e61d05de-34ef-4242-81a7-e6a96a1a68fb.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/4ad997c8-e970-424b-bc2e-4a48874d9c02.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/c6088cad-beaa-4afc-9998-4f37ed5a9fad.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/11e77350-eb0c-453f-8431-089b35dcfef3.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/4343bff4-b3db-4f3f-acf5-4ae7a2121793.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/4d6ff32d-4b6e-4b2f-9bfa-8077b1096c80.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/2cc5e74b-d976-4995-8dee-90369b249f42.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/97c93469-4806-4157-9afb-e9c127e144ce.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/08f6ec7d-531a-4a4c-88d7-265b90b11a59.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/b24dd6e1-2eb7-4c8e-a0aa-2606d1e44061.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/3df297ce-3c8d-4219-974e-6a6b891c9d60.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/352914e8-d52d-4f16-a903-1640537df8a7.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/c4f2c922-a9fa-4f24-a362-4e70bb44c194.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/02dfe367-a422-4a86-ac29-8efb56a97410.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/9457463c-3dc2-4d54-b6d1-6d93aa4ced4d.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/9b77f121-25b6-4374-9677-8e6509afb340.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/94250dd8-3979-460a-a179-bb593b6fbf37.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/2bddf8cf-601d-45be-9bf0-4f4467e99fef.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/a976264d-878a-4410-879a-bc9a9531ab19.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/5e7adab3-4400-4cce-aacd-2b24969f46cf.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/428a7e91-603d-41e0-895e-57fd4382fc89.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/87f69f01-336b-448e-a8f8-fd716ecbe56b.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/79f6c8ba-45bc-4f8a-a25c-890c62a5ef8c.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/f173806b-bac6-44a0-905d-00ee05cc2979.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/9b40cc77-fd24-40bb-bf9c-561388bc68a1.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/430d077b-e54d-4039-a181-0b1ae9e58aa9.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/3ec1a85b-8f7a-479d-841a-fe9478d4d016.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/4bbc340c-96d0-4535-bf17-f5d0af61b840.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/b75f1ab3-0483-4188-ba2b-3bead2d0fc44.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/01d90f1e-c5d7-4146-90df-78657ed6cf1b.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/31afd3e3-a148-423d-b884-3ce9c3af805b.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/2892cd44-6c49-450e-ace8-625d8400b933.root",",
#"/store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_123X_mcRun4_realistic_v3_2026D77PU200-v1/2580000/9257b49f-d498-43ed-865d-c77e36fba487.root",",    
#"file:#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_GEN_SIM_DIGI_RAW.root","
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_1.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_10.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_100.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_11.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_12.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_13.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_14.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_15.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_16.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_17.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_18.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_19.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_2.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_20.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_21.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_22.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_23.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_24.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_25.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_26.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_27.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_28.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_29.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_3.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_30.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_31.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_32.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_33.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_34.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_35.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_36.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_37.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_38.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_39.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_4.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_40.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_41.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_42.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_43.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_44.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_45.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_46.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_47.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_48.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_49.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_5.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_50.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_51.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_52.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_53.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_54.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_55.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_56.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_57.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_58.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_59.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_6.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_60.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_61.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_62.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_63.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_64.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_65.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_66.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_67.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_68.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_69.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_7.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_70.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_71.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_72.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_73.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_74.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_75.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_76.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_77.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_78.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_79.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_8.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_80.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_81.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_82.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_83.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_84.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_85.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_86.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_87.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_88.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_89.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_9.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_90.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_91.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_92.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_93.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_94.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_95.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_96.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_97.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_98.root",
#"file:/eos/user/g/gsaha4/Exotic/EvtGenforBstoPhiPhi/step3_output_99.root",
#"file:/afs/cern.ch/work/p/ppalit/public/HLSL1Trigger/EvtGenforBstoPhiPhi/step3/step3_withFEVT.root"
)
secFiles = cms.untracked.vstring()

process.source = cms.Source ("PoolSource",
                            fileNames = readFiles,
                            secondaryFileNames = secFiles,
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            )


process.TFileService = cms.Service("TFileService", fileName = cms.string('bstophiphiFromttbarinput.root'), closeFileFast = cms.untracked.bool(True))


############################################################
# L1 tracking: remake stubs?
############################################################

process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")

from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")

process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)
process.TTClusterStubTruth = cms.Path(process.TrackTriggerAssociatorClustersStubs)


# DTC emulation
process.load('L1Trigger.TrackerDTC.ProducerES_cff')
process.load('L1Trigger.TrackerDTC.ProducerED_cff')
process.dtc = cms.Path(process.TrackerDTCProducer)#*process.TrackerDTCAnalyzer)

process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")
process.load("L1Trigger.L1TTrackMatch.L1TrackSelectionProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1TrackJetProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1PhiMesonSelectionProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1BsMesonSelectionProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1GTTInputProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1KaonTrackSelectionProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1TrackJetEmulationProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1TrackFastJetProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1TrackerEtMissProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1TrackerEtMissEmulatorProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1TkHTMissProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1TkHTMissEmulatorProducer_cfi")
process.load('L1Trigger.VertexFinder.VertexProducer_cff')


############################################################
# Primary vertex
############################################################
process.L1VertexFinder = process.VertexProducer.clone()
process.pPV = cms.Path(process.L1VertexFinder)
process.L1VertexFinderEmulator = process.VertexProducer.clone()
process.L1VertexFinderEmulator.VertexReconstruction.Algorithm = "fastHistoEmulation"
process.L1VertexFinderEmulator.l1TracksInputTag = cms.InputTag("L1GTTInputProducer","Level1TTTracksConverted")
process.pPVemu = cms.Path(process.L1VertexFinderEmulator)

process.L1TrackFastJets.L1PrimaryVertexTag = cms.InputTag("L1VertexFinder", "l1vertices")
process.L1TrackFastJetsExtended.L1PrimaryVertexTag = cms.InputTag("L1VertexFinder", "l1vertices")
process.L1TrackJets.L1PVertexCollection = cms.InputTag("L1VertexFinder", "l1vertices")
process.L1TrackJetsExtended.L1PVertexCollection = cms.InputTag("L1VertexFinder", "l1vertices")
process.L1TrackerEtMiss.L1VertexInputTag = cms.InputTag("L1VertexFinder", "l1vertices")
process.L1TrackerHTMiss.L1VertexInputTag = cms.InputTag("L1VertexFinder", "l1vertices")
process.L1TrackerEtMissExtended.L1VertexInputTag = cms.InputTag("L1VertexFinder", "l1vertices")
process.L1TrackerHTMissExtended.L1VertexInputTag = cms.InputTag("L1VertexFinder", "l1vertices")
process.L1TrackerEmuEtMiss.L1VertexInputTag = cms.InputTag("L1VertexFinderEmulator", "l1verticesEmulation")


# HYBRID: prompt tracking
if (L1TRKALGO == 'HYBRID'):
    process.TTTracksEmu = cms.Path(process.L1HybridTracks)
    process.TTTracksEmuWithTruth = cms.Path(process.L1HybridTracksWithAssociators)
    process.pL1TrackSelection = cms.Path(process.L1TrackSelectionProducer)
    process.pL1KaonTrackSelection = cms.Path(process.L1KaonTrackSelectionProducer)
    #process.pL1KaonTrackSelection = cms.Path(process.L1TrackNullSelectionProducer*process.L1KaonTrackSelectionProducer)
    process.pL1PhiMesonSelection = cms.Path(process.L1PhiMesonSelectionProducer)
    process.pL1BsMesonSelection = cms.Path(process.L1BsMesonSelectionProducer)
    process.pL1TrackJets = cms.Path(process.L1TrackJets)
    process.pL1TrackFastJets=cms.Path(process.L1TrackFastJets)
    process.pL1GTTInput = cms.Path(process.L1GTTInputProducer)
    process.pL1TrackJetsEmu = cms.Path(process.L1TrackJetsEmulation)
    process.pTkMET = cms.Path(process.L1TrackerEtMiss)
    process.pTkMETEmu = cms.Path(process.L1TrackerEmuEtMiss)
    process.pTkMHT = cms.Path(process.L1TrackerHTMiss)
    process.pTkMHTEmulator = cms.Path(process.L1TrackerEmuHTMiss)
    DISPLACED = 'Prompt'

# HYBRID: extended tracking
elif (L1TRKALGO == 'HYBRID_DISPLACED'):
    process.TTTracksEmu = cms.Path(process.L1ExtendedHybridTracks)
    process.TTTracksEmuWithTruth = cms.Path(process.L1ExtendedHybridTracksWithAssociators)
    process.pL1TrackSelection = cms.Path(process.L1TrackSelectionProducerExtended)
    process.pL1KaonTrackSelection = cms.Path(process.L1KaonTrackSelectionProducerExtended)
    #process.pL1KaonTrackSelection = cms.Path(process.L1TrackNullSelectionProducerExtended*process.L1KaonTrackSelectionProducerExtended)
    process.pL1PhiMesonSelection = cms.Path(process.L1PhiMesonSelectionProducer)
    process.pL1BsMesonSelection = cms.Path(process.L1BsMesonSelectionProducer)
    process.pL1TrackJets = cms.Path(process.L1TrackJetsExtended)
    process.pL1TrackFastJets = cms.Path(process.L1TrackFastJetsExtended)
    process.pL1GTTInput = cms.Path(process.L1GTTInputProducerExtended)
    process.pL1TrackJetsEmu = cms.Path(process.L1TrackJetsExtendedEmulation)
    process.pTkMET = cms.Path(process.L1TrackerEtMissExtended)
    #process.pTkMETEmu = cms.Path(process.L1TrackerEmuEtMissExtended) #Doesn't exist
    process.pTkMHT = cms.Path(process.L1TrackerHTMissExtended)
    process.pTkMHTEmulator = cms.Path(process.L1TrackerEmuHTMissExtended)
    DISPLACED = 'Displaced'#

# HYBRID: extended tracking
elif (L1TRKALGO == 'HYBRID_PROMPTANDDISP'):
    process.TTTracksEmu = cms.Path(process.L1PromptExtendedHybridTracks)
    process.TTTracksEmuWithTruth = cms.Path(process.L1PromptExtendedHybridTracksWithAssociators)
    process.pL1TrackSelection = cms.Path(process.L1TrackSelectionProducer*process.L1TrackSelectionProducerExtended)
    process.pL1KaonTrackSelection = cms.Path(process.L1KaonTrackSelectionProducer*process.L1KaonTrackSelectionProducerExtended)
    #process.pL1KaonTrackSelection = cms.Path(process.L1TrackNullSelectionProducer*process.L1TrackNullSelectionProducerExtended*process.L1KaonTrackSelectionProducer*process.L1KaonTrackSelectionProducerExtended)
    process.pL1PhiMesonSelection = cms.Path(process.L1PhiMesonSelectionProducer)
    process.pL1BsMesonSelection = cms.Path(process.L1BsMesonSelectionProducer)
    process.pL1TrackJets = cms.Path(process.L1TrackJets*process.L1TrackJetsExtended)
    process.pL1TrackFastJets = cms.Path(process.L1TrackFastJets*process.L1TrackFastJetsExtended)
    process.pL1GTTInput = cms.Path(process.L1GTTInputProducer*process.L1GTTInputProducerExtended)
    process.pL1TrackJetsEmu = cms.Path(process.L1TrackJetsEmulation*process.L1TrackJetsExtendedEmulation)
    process.pTkMET = cms.Path(process.L1TrackerEtMiss*process.L1TrackerEtMissExtended)
    process.pTkMETEmu = cms.Path(process.L1TrackerEmuEtMiss)
    process.pTkMHT = cms.Path(process.L1TrackerHTMiss*process.L1TrackerHTMissExtended)
    process.pTkMHTEmulator = cms.Path(process.L1TrackerEmuHTMiss*process.L1TrackerEmuHTMissExtended)
    DISPLACED = 'Both'




############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# e.g. single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13
#      pions in jets = 6
#      taus = 15
#      all TPs = 1
############################################################

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackObjectNtupleMaker',
        MyProcess = cms.int32(1),
        DebugMode = cms.bool(False),      # printout lots of debug statements
        SaveAllTracks = cms.bool(True),  # save *all* L1 tracks, not just truth matched to primary particle
        SaveStubs = cms.bool(False),      # save some info for *all* stubs
        Displaced = cms.string(DISPLACED),# "Prompt", "Displaced", "Both"
        L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
        TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
        TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
        TP_minPt = cms.double(2.0),       # only save TPs with pt > X GeV
        TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
        TP_maxZ0 = cms.double(15.0),      # only save TPs with |z0| < X cm
        L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),                                                      # TTTracks, prompt
        L1TrackExtendedInputTag = cms.InputTag("TTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),                                      # TTTracks, extended
        MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"),                                               # MCTruth track, prompt
        MCTruthTrackExtendedInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigisExtended", "Level1TTTracks"),                               # MCTruth track, extended
        L1TrackGTTInputTag = cms.InputTag("L1GTTInputProducer","Level1TTTracksConverted"),                                                      # TTTracks, prompt, GTT converted
        L1TrackExtendedGTTInputTag = cms.InputTag("L1GTTInputProducerExtended","Level1TTTracksExtendedConverted"),                              # TTTracks, extended, GTT converted
        L1TrackSelectedInputTag = cms.InputTag("L1TrackSelectionProducer", "Level1TTTracksSelected"),                                           # TTTracks, prompt, selected
        L1TrackSelectedEmulationInputTag = cms.InputTag("L1TrackSelectionProducer", "Level1TTTracksSelectedEmulation"),                         # TTTracks, prompt, emulation, selected

        L1PosKaonTrackSelectedInputTag = cms.InputTag("L1KaonTrackSelectionProducer", "Level1TTKaonTracksSelectedPositivecharge"),                                           # TTTracks, prompt, selected
        L1PosKaonTrackSelectedEmulationInputTag = cms.InputTag("L1KaonTrackSelectionProducer", "Level1TTKaonTracksSelectedEmulationPositivecharge"),                         # TTTracks, prompt, emulation, selected                               
        L1NegKaonTrackSelectedInputTag = cms.InputTag("L1KaonTrackSelectionProducer", "Level1TTKaonTracksSelectedNegativecharge"),                                           # TTTracks, prompt, selected
        L1NegKaonTrackSelectedEmulationInputTag = cms.InputTag("L1KaonTrackSelectionProducer", "Level1TTKaonTracksSelectedEmulationNegativecharge"),                         # TTTracks, prompt, emulation, selected                               

        SaveTrackPhiCands = cms.bool(True), #includes emulated jets
        TrackPhiCandsInputTag = cms.InputTag("L1PhiMesonSelectionProducer", "Level1TTPhiMesonSelected"),                               

        SaveTrackBsCands = cms.bool(True), #includes emulated jets
        TrackBsCandsInputTag = cms.InputTag("L1BsMesonSelectionProducer", "Level1TTBsMesonSelected"),                               

        L1TrackExtendedSelectedInputTag = cms.InputTag("L1TrackSelectionProducerExtended", "Level1TTTracksExtendedSelected"),                   # TTTracks, extended, selected
        L1TrackExtendedSelectedEmulationInputTag = cms.InputTag("L1TrackSelectionProducerExtended", "Level1TTTracksExtendedSelectedEmulation"), # TTTracks, extended, emulation, selected

        L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
        MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
        MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
        TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
        TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
        GenJetInputTag = cms.InputTag("ak4GenJets", ""),
        ##track jets and track MET
        SaveTrackJets = cms.bool(True), #includes emulated jets
        SaveTrackSums = cms.bool(True), #includes simulated/emulated track MET, MHT, and HT
        TrackFastJetsInputTag = cms.InputTag("L1TrackFastJets","L1TrackFastJets"),
        TrackFastJetsExtendedInputTag = cms.InputTag("L1TrackFastJetsExtended","L1TrackFastJetsExtended"),
        TrackJetsInputTag = cms.InputTag("L1TrackJets", "L1TrackJets"),
        TrackJetsExtendedInputTag=cms.InputTag("L1TrackJetsExtended", "L1TrackJetsExtended"),
        TrackJetsEmuInputTag = cms.InputTag("L1TrackJetsEmulation","L1TrackJets"),
        TrackJetsExtendedEmuInputTag = cms.InputTag("L1TrackJetsExtendedEmulation","L1TrackJetsExtended"),
        TrackMETInputTag = cms.InputTag("L1TrackerEtMiss","L1TrackerEtMiss"),
        TrackMETExtendedInputTag = cms.InputTag("L1TrackerEtMissExtended","L1TrackerExtendedEtMiss"),
        TrackMETEmuInputTag = cms.InputTag("L1TrackerEmuEtMiss","L1TrackerEmuEtMiss"),
        TrackMHTInputTag = cms.InputTag("L1TrackerHTMiss","L1TrackerHTMiss"), #includes HT
        TrackMHTExtendedInputTag = cms.InputTag("L1TrackerHTMissExtended","L1TrackerHTMissExtended"),
        TrackMHTEmuInputTag = cms.InputTag("L1TrackerEmuHTMiss",process.L1TrackerEmuHTMiss.L1MHTCollectionName.value()),
        TrackMHTEmuExtendedInputTag = cms.InputTag("L1TrackerEmuHTMissExtended",process.L1TrackerEmuHTMissExtended.L1MHTCollectionName.value()),
        GenParticleInputTag = cms.InputTag("genParticles",""),
        RecoVertexInputTag=cms.InputTag("L1VertexFinder", "l1vertices"),
        RecoVertexEmuInputTag=cms.InputTag("L1VertexFinderEmulator", "l1verticesEmulation"),
)

process.ntuple = cms.Path(process.L1TrackNtuple)

process.out = cms.OutputModule( "PoolOutputModule",
                                fastCloning = cms.untracked.bool( False ),
                                fileName = cms.untracked.string("test.root" )
		               )
process.pOut = cms.EndPath(process.out)


# use this if you want to re-run the stub making
# process.schedule = cms.Schedule(process.TTClusterStub,process.TTClusterStubTruth,process.TTTracksEmuWithTruth,process.ntuple)

# use this if cluster/stub associators not available
# process.schedule = cms.Schedule(process.TTClusterStubTruth,process.TTTracksEmuWithTruth,process.ntuple)

process.schedule = cms.Schedule(process.TTClusterStub, process.TTClusterStubTruth, process.dtc, process.TTTracksEmuWithTruth, process.pL1GTTInput, process.pPV, process.pPVemu, process.pL1TrackSelection, process.pL1KaonTrackSelection, process.pL1PhiMesonSelection, process.pL1BsMesonSelection, process.pL1TrackJets, process.pL1TrackJetsEmu, process.pL1TrackFastJets, process.pTkMET, process.pTkMETEmu, process.pTkMHT, process.pTkMHTEmulator, process.ntuple)

#process.schedule = cms.Schedule(process.TTTracksEmuWithTruth, process.pL1GTTInput, process.pPV, process.pPVemu, process.pL1TrackSelection, process.pL1KaonTrackSelection, process.pL1PhiMesonSelection, process.pL1BsMesonSelection, process.pL1TrackJets, process.pL1TrackJetsEmu, process.pL1TrackFastJets, process.pTkMET, process.pTkMETEmu, process.pTkMHT, process.pTkMHTEmulator, process.ntuple)
