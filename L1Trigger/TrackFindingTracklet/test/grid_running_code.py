import os
import sys
from itertools import islice

############## Things to Change ###################
# dataName = 'TTbarPU200D88_HYBRID_TCTrunc_assertsOn_oneTrunc' # what submit files and output files are named (no ".root")
# filename = 'TTbarPU200_all.txt' ## location of the txt file used to help locate samples
# cfgFile = 'L1TrackNtupleMaker_cfg_grid_HYBRID_trkJets.py' # cfg file used to run each job

dataName = 'TTbarPU200D88_NEWKF_TCTrunc_assertsOn_oneTrunc' # what submit files and output files are named (no ".root")
filename = 'TTbarPU200_all.txt' ## location of the txt file used to help locate samplefs
cfgFile = 'L1TrackNtupleMaker_cfg_grid_NEWKF_trkJets.py' # cfg file used to run each job
###################################################

submit = 'universe = vanilla\n' ## writing .sub file
submit += 'arguments = "$(argument)"\n'
submit += 'output = NtupleGridOutput/' + dataName + '/submit$(ClusterId).out\n' # where to output outputs
submit += 'error = NtupleGridOutput/' + dataName + '/submit$(ClusterId).err\n'
submit += 'log = NtupleGridOutput/' + dataName + '/submit$(ClusterId).log\n'
submit += '+JobFlavour = "tomorrow"\n' 
submit += 'queue\n'
submitName = 'NtupleGridOutput/' + dataName + '/submit01.sub'
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
        create += 'export CMSSW_PROJECT_SRC=/afs/cern.ch/user/d/dabraham/private/test/CMSSW_12_6_0_pre5/src\n' 
        create += 'cd $CMSSW_PROJECT_SRC\n' ## go to the src directly 
        create += 'eval `scramv1 runtime -sh`\n' ## this is effectively "cmsenv"
        create += 'export X509_USER_PROXY=/afs/cern.ch/user/d/dabraham/x509up_u161126\n' ## exporting the proxy
        create += 'cd /afs/cern.ch/user/d/dabraham/private/test/CMSSW_12_6_0_pre5/src/L1Trigger/TrackFindingTracklet/test\n'## go the directly where you have the producer and python config file
        create += 'cmsRun ' + cfgFile + ' inputFiles='+input+' outputFile=/eos/user/d/dabraham/L1NtupleTrackExamples/' + dataName + '_'+str(counter1)+'.root maxEvents=10000\n'
        createName = 'NtupleGridOutput/' + dataName + '/submit'+str(counter1)+'.sh'
        sub2 = open(createName,'w')
        sub2.write(create+'\n')
        sub2.close() ## finish writing .sh file
        counter1+=1
        os.system('chmod 755 '+createName) ## make .sh file executable
        os.system('condor_submit '+ submitName+' executable='+createName) ## submit the job using condor_submit command
