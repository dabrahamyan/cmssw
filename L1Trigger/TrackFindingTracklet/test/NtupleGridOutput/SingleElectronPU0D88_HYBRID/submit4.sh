#!/bin/bash
export CMSSW_PROJECT_SRC=/afs/cern.ch/user/d/dabraham/private/test/CMSSW_12_6_0_pre5/src
cd $CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
export X509_USER_PROXY=/afs/cern.ch/user/d/dabraham/x509up_u161126
cd /afs/cern.ch/user/d/dabraham/private/test/CMSSW_12_6_0_pre5/src/L1Trigger/TrackFindingTracklet/test
cmsRun L1TrackNtupleMaker_cfg_grid_HYBRID.py inputFiles=root://cms-xrd-global.cern.ch///store/relval/CMSSW_12_6_0/RelValSingleElPt2p0to100p0/GEN-SIM-DIGI-RAW/125X_mcRun4_realistic_v5_2026D88noPURV183-v1/2590000/13679a5d-c3a0-459a-bc95-c5ba3d65e082.root outputFile=/eos/user/d/dabraham/L1NtupleTrackExamples/SingleElectronPU0D88_HYBRID_4.root maxEvents=10000

