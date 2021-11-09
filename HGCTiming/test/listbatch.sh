#!/bin/bash
doSubmit=true
modes=( \ 
 "eol"       \
# "startup"       \
) 
num=1500
upnum=25
count=0
ind=0

makealistdir () {
# write base for submit file
#until [ ${base} -gt ${num} ]
#do
#submitfile="${mode}${ind}.list"
# until [ ${count} -gt ${upnum} ]
# do
# printf "root://cms-xrd-global.cern.ch//store/group/dpg_hgcal/comm_hgcal/suhokim/RECO/GEN_13Pt10_Vtx0_flatEta_1p5_1p8_26D49_${base}_${mode}_step3.root\n" >> submitfile
# count=$(( ${count} + 1 ))
# base=$(( ${base} + 1 ))
# done
# if [ ${count} = ${upnum} ]
# then
#  count=0
#  ind=$((ind+1))
# fi
#popd > /dev/null
#done

for i in {1..500}
do
 submitfile="${mode}${ind}.list"
 printf 'root://cmsxrootd.fnal.gov//store/user/skim2/GEN_13Pt10_Vtx0_flatEta_1p5_1p8_26D49_'${i}'_'${mode}'_step3.root\n' >> "$submitfile"
 #printf 'root://cms-xrd-global.cern.ch//store/group/dpg_hgcal/comm_hgcal/suhokim/RECO/GEN_13Pt10_Vtx0_flatEta_1p5_1p8_26D49_'${i}'_'${mode}'_step3.root\n' >> "$submitfile"
 count=$(( ${count} + 1 ))
 if [ ${count} = ${upnum} ]
 then
  count=0
  ind=$((ind+1))
 fi
done


}

for mode in ${modes[@]}
do 
 makealistdir ${num} ${upnum} ${mode}
done

