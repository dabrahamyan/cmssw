#!/bin/bash

function transfer () {
    cp /eos/user/d/dabraham/code_backups/CMSSW_12_6_0_pre5_2023_10_04/src/L1Trigger/TrackFindingTracklet/test/$1 /afs/cern.ch/user/d/dabraham/private/test/CMSSW_12_6_0_pre5/src/L1Trigger/TrackFindingTracklet/test/$1
}

transfer myProjects/overlayHists_combinedvsNot.C