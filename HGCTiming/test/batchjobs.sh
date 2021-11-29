#!/bin/bash
python Setupmake.py
doSubmit=true
modes=( \ 
 "Rand40mu_eol"       \
# "NeutrinoGun_startup_sn2.0"       \
# "Rand40mu_startup_sn2.0"       \
# "40mu_startup_sn2.0"       \
# "old_startup_sn2.0"       \
# "old_startup_sn2.5"       \
# "old_startup_sn3.0"       \
# "old_startup_sn4.0"       \
# "old_startup"       \
# "old_eol"       \
# "old_startup_sn2.0"       \
# "Gun50"       \
# "startup"       \
) 
num=0
upnum=44

makeasubmitdir () {
# write base for submit file
 printf "Making submits for $3\n"

 # go to the directory
 origindir=$(pwd)
 submitdir=$(pwd)/gitignore/$3/
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
 printf "Transfer_Input_Files = ${origindir}/CMSSW_12_1_0_pre4.tar.gz,${origindir}/runTest.py,${origindir}/lists/${mode}_${num}.list\n" >> submitfile 
 printf "Arguments = inputFile=${mode}_${num}\n" >> submitfile
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

