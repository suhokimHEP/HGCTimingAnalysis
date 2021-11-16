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
tar -xzf CMSSW_12_1_0_pre4.tar.gz
rm -fv CMSSW_12_1_0_pre4.tar.gz
export SCRAM_ARCH=slc7_amd64_gcc900
cd HGCTiming_CMSSW_12_1_0_pre4/
scramv1 b ProjectRename
eval `scram runtime -sh`
cd ../
mkdir -p plotdir/cellsPhi
sample=$1"_"$2
echo "now run"
echo $1
echo $2
echo ${sample}
root -l -b 'lookAtMIP_correctRate.C('\""${sample}"\"','\""$1"\"')'
echo "after run ls"
ls
