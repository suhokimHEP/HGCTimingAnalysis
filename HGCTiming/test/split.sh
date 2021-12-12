#!/bin/bash
input=text
output=Rand40mu
mode=startup_sn4.0
HGCP=/eos/cms/store/group/dpg_hgcal/comm_hgcal/suhokim/${output}
find ${HGCP}/RECO/*${mode}_st*root > ${input}.list
#sed -i "s/\/eos\/cms/root:\/\/cms-xrd-global.cern.ch\//g" ${input}.list
sed -i "s/\/eos\/cms/root:\/\/cmsxrootd.fnal.gov\//g" ${input}.list
split -l 10 --numeric-suffixes --suffix-length=2 --additional-suffix=.list ${input}.list ${output}_${mode}_
rename ${output}_${mode}_0 ${output}_${mode}_ *
mv ${output}*list lists/
rm ${input}.list
