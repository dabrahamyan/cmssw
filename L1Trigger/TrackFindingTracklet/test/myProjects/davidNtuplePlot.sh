#!/bin/bash

filename='DisplacedMuon_PU0_D88_DISPLACED_Uncomb_BugFix_2023_11_8'

filedir='/eos/user/d/dabraham/DisplacedCombinedBugFix/'

root -l -b -q "davidNtuplePlot_displaced.C(\"$filename\", \"$filedir\", \"\", 0)" > NtuplePlotOutput/${filename}.out 

echo "donezo"
