#!/bin/bash

fileName="combinedDebug_SingleMuon_set1_HYBRID_Aug8commit"
dirName="/eos/user/d/dabraham/L1NtupleTrackExamples/"
numbers='0'

for number in $numbers
do
    root -l -b -q "davidNtuplePlot.C(\"${fileName}\", \"${dirName}\", \"\", ${number})"  | tail -n 19 > NtuplePlotOutput/${fileName}_${number}.out
done

echo "donezo"
