#!/bin/bash
input=text
output=NeutrinoGun
mode=startup_sn2.0
HGCP=/eos/cms/store/group/dpg_hgcal/comm_hgcal/suhokim/${output}
find ${HGCP}/RECO/*${mode}*root > ${input}.list
#sed -i "s/\/eos\/cms/root:\/\/cms-xrd-global.cern.ch\//g" ${input}.list
sed -i "s/\/eos\/cms/root:\/\/cmsxrootd.fnal.gov\//g" ${input}.list
split -l 10 --numeric-suffixes --suffix-length=2 --additional-suffix=.list ${input}.list ${output}_${mode}_
rename ${output}_${mode}_0 ${output}_${mode}_ *
mv ${output}*list lists/
rm ${input}.list
