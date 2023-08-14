#!/bin/bash
export CMSSW_PROJECT_SRC=/afs/cern.ch/user/d/dabraham/private/test/CMSSW_12_6_0_pre5/src
cd $CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
export X509_USER_PROXY=/afs/cern.ch/user/d/dabraham/x509up_u161126
cd /afs/cern.ch/user/d/dabraham/private/test/CMSSW_12_6_0_pre5/src/L1Trigger/TrackFindingTracklet/test
cmsRun L1TrackNtupleMaker_cfg_grid_NEWKF.py inputFiles=root://cms-xrd-global.cern.ch///store/relval/CMSSW_12_6_0/RelValSingleMuPt2p0to100p0/GEN-SIM-DIGI-RAW/125X_mcRun4_realistic_v5_2026D88noPURV183-v1/2590000/9948e49d-c40f-43fc-a050-1f74b3a04c7b.root outputFile=/eos/user/d/dabraham/L1NtupleTrackExamples/SingleMuonPU0D88_NEWKF_11.root maxEvents=10000

