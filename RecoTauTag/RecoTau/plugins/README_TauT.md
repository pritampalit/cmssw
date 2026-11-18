# TauTransformerIdONNX

Computes the Tau Transformer (TauT) scores for HPS taus. Model: `RecoTauTag/RecoTau/data/taut_tim-ep59.onnx` (Tim's TauT).
Scores are between 0 and 1 for HPS taus and -1 for all other taus.

Outputs: `VSjet`, `VSe`, `VSmu` (tau vs jet / electron / muon).

- **Part 1 (MiniAOD taus):** `toKeep = ["tauT"]` in `RecoTauTag/RecoTau/python/tools/runTauIdMVA.py` adds the scores to the taus
  as `byTauTVSjetraw`, `byTauTVSeraw`, `byTauTVSmuraw`.
- **Part 2 (additional, NanoAOD):** `--customise PhysicsTools/NanoAOD/nano_cff.nanoAOD_addTauT` in the NANO step gives
  `Tau_rawTauTVSjet`, `Tau_rawTauTVSe` and `Tau_rawTauTVSmu`.

Install and examples: `TauT_README.md` at the top of this branch.
