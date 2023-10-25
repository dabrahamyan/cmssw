#!/bin/bash

filename='DisplacedMuon_PU0_D88_DISPLACED_Comb_LatestDev_2023_10_24'

filedir='/eos/user/d/dabraham/L1NtupleTrackExamples/'

root -l -b -q "davidNtuplePlot_displaced.C(\"$filename\", \"$filedir\", \"\", 0)" | tail -n 20 > NtuplePlotOutput/${filename}.out 

echo "donezo"
