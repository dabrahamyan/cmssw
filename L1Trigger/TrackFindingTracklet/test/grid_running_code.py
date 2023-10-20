import os
import sys
from itertools import islice

############## Things to Change ###################
# dataName = 'TTbar_PU200_D88_HYBRID_LatestDev_2023_10_20' # what submit files and output files are named (no ".root")
# filename = 'TTbarPU200_all.txt' ## location of the txt file used to help locate samples
# cfgFile = 'L1TrackNtupleMaker_cfg_grid_HYBRID.py' # cfg file used to run each job

dataName = 'DisplacedMuon_PU200_D88_DISPLACED_LatestDev_2023_10_20' # what submit files and output files are named (no ".root")
filename = 'Muon_Displaced_all.txt' ## location of the txt file used to help locate samples
cfgFile = 'L1TrackNtupleMaker_cfg_grid_DISPLACED.py' # cfg file used to run each job
######################### Change to CMSSW you are using #########################

newCMSSW = 'CMSSW_13_3_0_pre2'  # tempCMSSWs/two  tempCMSSWs/one  CMSSW_12_6_0_pre5

###################################################################################
origPath = '/afs/cern.ch/user/d/dabraham/private/new_2023_10_10/CMSSW_13_3_0_pre2/src/L1Trigger/TrackFindingTracklet/test/'

submit = 'universe = vanilla\n' ## writing .sub file
submit += 'arguments = "$(argument)"\n'
submit += 'output = ' + origPath + 'NtupleGridOutput/' + dataName + '/submit$(ClusterId).out\n' # where to output outputs
submit += 'error = ' + origPath + 'NtupleGridOutput/' + dataName + '/submit$(ClusterId).err\n'
submit += 'log = ' + origPath + 'NtupleGridOutput/' + dataName + '/submit$(ClusterId).log\n'
submit += '+JobFlavour = "tomorrow"\n' 
submit += 'queue\n'
submitName = origPath + 'NtupleGridOutput/' + dataName + '/submit01.sub'
sub1 = open(submitName,'w')
sub1.write(submit+'\n')
sub1.close() ##finish writing .sub file
nfile = 1 ##number of sample file processed per each condor job
with open(filename,'r') as f:
    counter1 = 1 ## the nth job
    while True:
        lines = list(islice(f, nfile))
        if not lines:
            break
        counter2 = 1 ## the nth file in one job
        for line in lines:
            if counter2 == 1:
               input='root://cms-xrd-global.cern.ch//'+line.rstrip() ## Remove  'root:..'+  if you don't have grid certificate. 
            else:
               input =input+',root://cms-xrd-global.cern.ch//'+line.rstrip() ## Remove  ',root:..'+  if you don't have grid certificate.
            counter2+=1
        create = '#!/bin/bash\n' ##writng .sh file
        create += 'export CMSSW_PROJECT_SRC=/afs/cern.ch/user/d/dabraham/private/new_2023_10_10/' + newCMSSW + '/src\n' 
        create += 'cd $CMSSW_PROJECT_SRC\n' ## go to the src directly 
        create += 'eval `scramv1 runtime -sh`\n' ## this is effectively "cmsenv"
        create += 'export X509_USER_PROXY=/afs/cern.ch/user/d/dabraham/x509up_u161126\n' ## exporting the proxy
        create += 'cd /afs/cern.ch/user/d/dabraham/private/new_2023_10_10/' + newCMSSW + '/src/L1Trigger/TrackFindingTracklet/test\n'## go the directly where you have the producer and python config file
        create += 'cmsRun ' + cfgFile + ' inputFiles='+input+' outputFile=/eos/user/d/dabraham/L1NtupleTrackExamples/' + dataName + '_'+str(counter1)+'.root maxEvents=10000\n'
        createName = origPath + 'NtupleGridOutput/' + dataName + '/submit'+str(counter1)+'.sh'
        sub2 = open(createName,'w')
        sub2.write(create+'\n')
        sub2.close() ## finish writing .sh file
        counter1+=1
        os.system('chmod 755 '+createName) ## make .sh file executable
        os.system('condor_submit '+ submitName+' executable='+createName) ## submit the job using condor_submit command
