#!/bin/bash

function transfer () {
    cp /eos/user/d/dabraham/code_backups/CMSSW_12_6_0_pre5_2023_10_04/src/L1Trigger/TrackFindingTracklet/test/$1 /afs/cern.ch/user/d/dabraham/private/test/CMSSW_12_6_0_pre5/src/L1Trigger/TrackFindingTracklet/test/$1
}

transfer SingleElectronPU0_all.txt
transfer SingleMuonPU0_all.txt
transfer SingleMuonPU0_set1.txt
transfer TTbarPU200_all.txt
transfer TTbarPU200_set1.txt