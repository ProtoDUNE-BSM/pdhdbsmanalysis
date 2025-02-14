#!/bin/bash
mrb uc
mrbslp
mrbsetenv
mrb i -j10

# /pnfs/dune/scratch/users/chasnip/ProtoDUNEBSM/NeutrinoSim_v2/wnp04/numu/decay0/reco/wnp04_numu_decaymode0_recotriggereddetsim2g4gen_26_1_16477398_0.root
# /pnfs/dune/scratch/users/chasnip/ProtoDUNEBSM/NeutrinoSim_v2/wnp04/numu/decay0/reco/wnp04_numu_decaymode0_recotriggereddetsim2g4gen_84_1_17341409_0.root


# if something goes wrong, exit
if [ $? -ne 0 ]; then
  echo "Error: mrb i -j10 failed"
else
  # Run the analysis
  cd /exp/dune/app/users/dpullia/protodunedm_analysis/srcs/LArSoftAnalysisTools/pdhdbsmanalysis
  # lar -c run_analyzeEvents.fcl -s /exp/dune/data/users/dpullia/reco_files/nu_events/np04hd_raw_run029424_0020_dataflow1_datawriter_0_20241004T180407_reco_stage1_reco_stage2_20241004T211145_keepup.root -n 5
  # lar -c run_analyzeEvents.fcl -s /pnfs/dune/scratch/users/achatter/ProtoDUNEBSM/NeutrinoSim/wnp04/numu/decay0/reco/wnp04_numu_decaymode0_recotriggereddetsim2g4gen_37_1_59137416_0.root 
  lar -c run_analyzeEvents.fcl -s /pnfs/dune/scratch/users/chasnip/ProtoDUNEBSM/NeutrinoSim_v2/wnp04/numu/decay0/reco/wnp04_numu_decaymode0_recotriggereddetsim2g4gen_26_1_16477398_0.root
  # lar -c run_analyzeEvents.fcl -s /pnfs/dune/scratch/users/chasnip/ProtoDUNEBSM/NeutrinoSim_v2/wnp04/numu/decay0/reco/wnp04_numu_decaymode0_recotriggereddetsim2g4gen_84_1_17341409_0.root
  cd /exp/dune/app/users/dpullia/protodunedm_analysis/srcs/LArSoftAnalysisTools
fi

# lar -c pdhdbsmanalysis/eventdump.fcl -s /exp/dune/data/users/dpullia/reco_files/nu_events/np04hd_raw_run029424_0020_dataflow1_datawriter_0_20241004T180407_reco_stage1_reco_stage2_20241004T211145_keepup.root -n 5
# lar -c pdhdbsmanalysis/eventdump.fcl -s /pnfs/dune/scratch/users/chasnip/ProtoDUNEBSM/NeutrinoSim_v2/wnp04/numu/decay0/reco/wnp04_numu_decaymode0_recotriggereddetsim2g4gen_26_1_16477398_0.root -n 5

# Begin processing the 5th record. run: 29424 subRun: 1 event: 54959 at 04-Dec-2024 04:05:23 CST


