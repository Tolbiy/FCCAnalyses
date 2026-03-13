## Example search for Bs to tau tau decay

### Overview

Analysis set to study sensitivity for Bs to tau tau decay.
Initial commit based on BuBc2TauNu analysis, focusing on tau to 3pi decays.  
To run the `analysis_script.py` described below scripts, use `fccanalysis run analysis_script.py`. Check the [FCCAnalyses tutorial](https://hep-fcc.github.io/fcc-tutorials/main/index.html) for more details on how to write these scripts.

### Scripts
`analysis_stage1.py` Stage 1 sample production for Bs->Tau(->3PiNu)Tau(->3PiNu). Creates the event-level variables columns and Tau->3Pi candidates. Selects two in the signal hemisphere and compute the colinear mass, expected to be used for signal extraction.  

`mumu/analysis_stage1.py` Stage 1 sample production for Bs->Tau(->MuNuNu)Tau(->MuNuNu). Creates the event-levels variables columns and muon collection properties as well as electron truth-matching. Defines a dimuon system, the two best muon candidates to come from the signal decay, and Stage 1 cuts flags.  
`mumu/anaylsis_stage2_bkg` Used to study the MC background that passes stage1 cuts. Reconstruct the MC decay that generated the dimuon system.  

`ee/analysis_stage1.py` Stage 1 sample production for Bs->Tau(->ENuNu)Tau(->ENuNu). Creates the event-levels variables columns and electron collection properties as well as electron truth-matching. Defines a dielectron system, the two best electron candidates to come from the signal decay, and Stage 1 cuts flags.  
`ee/analysis_stage1_TMStudy.py` Very similar to `mumu/anaylsis_stage2_bkg` and used to study a duplication of the MC truth-matched indices.  

All repo have a `functions.h` which includes additional analyzer functions used to find the dilepton systems, perform truth-matching, define Stage 1 cuts flags.  

