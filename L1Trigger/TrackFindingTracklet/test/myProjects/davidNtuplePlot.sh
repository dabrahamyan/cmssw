#!/bin/bash

filename='retryTest_TTbar_PU200_D88_HYBRID_Comb_LatestDev_2023_10_20'

filedir='/eos/user/d/dabraham/L1NtupleTrackExamples/'

root -l -b -q "davidNtuplePlot.C(\"$filename\", \"$filedir\", \"\", 0)" | tail -n 20 > NtuplePlotOutput/${filename}.out 

echo "donezo"
