#!/bin/bash
OUT_FOLDER=/exp/dune/app/users/dpullia/protodunedm_analysis/out_products
# /exp/dune/data/users/dpullia/reco_files/nu_events
# iterate over all the files in the nu events folder
for file in /exp/dune/data/users/dpullia/reco_files/nu_events/*; do
    cd /exp/dune/app/users/dpullia/protodunedm_analysis/srcs/LArSoftAnalysisTools/pdhdbsmanalysis
    lar -c run_analyzeEvents.fcl -s $file 
    mv analysisOutput.root $OUT_FOLDER/$(basename $file .root)_analysisOutput.root
    cd /exp/dune/app/users/dpullia/protodunedm_analysis/srcs/LArSoftAnalysisTools
    done

