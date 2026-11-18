# TauT (Tau Transformer) scores for HPS taus

Branch for tau ID studies, based on CMSSW_17_0_0_pre5. Model: Tim's TauT (`RecoTauTag/TrainingFiles/data/TauTransformer/taut_tim-ep59.onnx`).
Scores are between 0 and 1 for HPS taus and -1 for all other taus.

## Install

```bash
cmsrel CMSSW_17_0_0_pre5 && cd CMSSW_17_0_0_pre5/src && cmsenv
git cms-init
git cms-merge-topic pritampalit:taut-nano-study
scram b -j 8
```

## Part 1: TauT on MiniAOD taus

```python
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, originalTauName = "slimmedTaus",
                                          updatedTauName = "slimmedTausWithTauT", toKeep = ["tauT"])
tauIdEmbedder.runTauID()
process.p = cms.Path(process.rerunMvaIsolationSequence * process.slimmedTausWithTauT)
```

The taus in `slimmedTausWithTauT` then have `tauID("byTauTVSjetraw")`, `tauID("byTauTVSeraw")`, `tauID("byTauTVSmuraw")`
and `tauID("byTauTProbTau")`, `tauID("byTauTProbEle")`, `tauID("byTauTProbMu")`, `tauID("byTauTProbJet")`.

## Part 2 (additional): TauT in NanoAOD

Add to the NANO step (MiniAOD to NanoAOD):

```bash
--customise PhysicsTools/NanoAOD/nano_cff.nanoAOD_addTauT
```

New Tau branches: `Tau_rawTauTVSjet`, `Tau_rawTauTVSe`, `Tau_rawTauTVSmu`, `Tau_probTauTTau`, `Tau_probTauTEle`, `Tau_probTauTMu`, `Tau_probTauTJet`.
Without the customisation these branches are -1.
