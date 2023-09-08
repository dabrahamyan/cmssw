#!/bin/bash

echo "Starting job on `date`" # date/time of start of job
echo "Running on: `uname -a`" # condor job is running on this node
echo "System software: `cat /etc/redhat-release`" # operating system on that node

source /cvmfs/cms.cern.ch/cmsset_default.sh # to get cmsenv

export X509_USER_PROXY=/afs/cern.ch/user/d/dabraham/.globus/x509up_u161126 # for grid stuff
export EOS_MGM_URL=root://eosuser.cern.ch

cp /path/to/CMSSW_10_6_30_patch1.tgz ./CMSSW_10_6_30_patch1.tgz
tar -xf CMSSW_10_6_30_patch1.tgz
rm CMSSW_10_6_30_patch1.tgz
cd CMSSW_10_6_30_patch1/src/
scramv1 b ProjectRename # this handles linking the already compiled code - do NOT recompile
cmsenv
echo $CMSSW_BASE "is the CMSSW on the local worker node"