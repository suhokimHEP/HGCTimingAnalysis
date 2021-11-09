#!/bin/bash 

echo "TEST"
voms-proxy-info --all
ls -l
echo "DONE"
outDir="/eos/uscms/store/user/skim2/"
echo output directory, $outDir
export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

source /cvmfs/cms.cern.ch/cmsset_default.sh
tar -xzf CMSSW_12_1_0_pre3.tar.gz
rm -fv CMSSW_12_1_0_pre3.tar.gz
export SCRAM_ARCH=slc7_amd64_gcc900
cd CMSSW_12_1_0_pre3/
scramv1 b ProjectRename
eval `scram runtime -sh`
cd ../

echo "doing pwd"
pwd
echo "doing ls"
ls
echo "now run"
echo cmsRun
cmsRun runTest.py $1
echo "after run ls"
ls

#for FILE in *.root
#do
# xrdcp -f ${FILE} root://cmseos.fnal.gov//store/user/skim2/
#  rm ${FILE}
#done

