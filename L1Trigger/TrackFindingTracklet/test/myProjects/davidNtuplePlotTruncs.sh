#!/bin/bash

# This code is made to run davidNtuplePlot on my Truncation Examination outputs.
# The goal is to automate running it over all of the ntuples I make for the different 
# truncation settings and of the hybrid and hybrid_newkf algos and all 4 injet settings 
# which are no in jet, regular injet, highpt, vhighpt (0, 1, 2, 3)
#
# By David Abrahamyan Sep 7, 2023

# numbers for the injet settings
numbers='0 1'

# algorithms were testing
algos="HYBRID NEWKF"

# different truncation 
truncs="noMETrunc noTETrunc"

# Where to get the NtupleRoot files from 
dirName="/eos/user/d/dabraham/L1NtupleTrackExamples/"

for algo in $algos
do
    for trunc in $truncs
    do 
        for number in $numbers
        do
            # The root file you'll operate on (w/o ".root")
            fileName="TTbar_PU200_D88_${algo}_${trunc}_assertsOn"
            # Open root, run davidNtuplePlot, quit root
            root -l -b -q "davidNtuplePlot.C(\"${fileName}\", \"${dirName}\", \"\", ${number})"
        done
    done
done

echo "donezo"
