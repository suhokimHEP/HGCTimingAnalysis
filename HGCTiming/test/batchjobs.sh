#!/bin/bash
doSubmit=true
modes=( \ 
 "eol"       \
# "startup"       \
) 
num=0
upnum=2

makeasubmitdir () {
# write base for submit file
 printf "Making submits for $3/$1\n"

 # go to the directory
 origindir=$(pwd)
 submitdir=$(pwd)/gitignore/$3/$1
 mkdir -p ${submitdir}
 pushd    ${submitdir}  > /dev/null
 printf " The directory is %s\n" $(pwd)

 mkdir -p logs



 printf "Universe = vanilla\n" > submitfile
 printf "Executable = ${origindir}/run_job.sh\n" >> submitfile
 printf "Should_Transfer_Files = YES \n" >> submitfile
 printf "WhenToTransferOutput = ON_EXIT\n" >> submitfile

 printf "notify_user = skim2@cern.ch\n" >> submitfile
 printf "\n" >> submitfile
 printf "Output = logs/SN_\$(Cluster)_\$(Process).stdout\n" >> submitfile
 printf "Error  = logs/SN_\$(Cluster)_\$(Process).stderr\n" >> submitfile
 printf "Log    = logs/SN_\$(Cluster)_\$(Process).log\n" >> submitfile
 printf "\n" >> submitfile
 until [ ${num} -gt ${upnum} ]
 do
 printf "Transfer_Input_Files = ${origindir}/CMSSW_12_1_0_pre4.tar.gz,${origindir}/runTest.py,${origindir}/lists/${mode}${num}.list\n" >> submitfile 
 printf "Arguments = inputFile=${mode}${num}\n" >> submitfile
 printf "Queue\n" >> submitfile
 printf "\n" >> submitfile
 printf "\n" >> submitfile
 printf "\n" >> submitfile
 num=$(( ${num} + 1 ))
 done
 if [ ${doSubmit} = true ]
 then
  condor_submit submitfile
 fi
 popd > /dev/null



}

for mode in ${modes[@]}
do 
 makeasubmitdir ${num} ${upnum} ${mode}
done
