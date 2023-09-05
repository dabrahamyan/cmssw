#!/bin/bash
export CMSSW_PROJECT_SRC=/afs/cern.ch/user/d/dabraham/private/test/CMSSW_12_6_0_pre5/src
cd $CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
export X509_USER_PROXY=/afs/cern.ch/user/d/dabraham/x509up_u161126
cd /afs/cern.ch/user/d/dabraham/private/test/CMSSW_12_6_0_pre5/src/L1Trigger/TrackFindingTracklet/test
cmsRun L1TrackNtupleMaker_cfg_grid_NEWKF_trkJets.py inputFiles=root://cms-xrd-global.cern.ch///store/mc/CMSSW_12_6_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_125X_mcRun4_realistic_v5_2026D88PU200RV183v2-v1/30000/0976683f-a9a5-4ba6-aa7d-307235d7448c.root outputFile=/eos/user/d/dabraham/L1NtupleTrackExamples/TTbarPU200D88_NEWKF_combined_fullTrunc_assertsOn_2.root maxEvents=10000

