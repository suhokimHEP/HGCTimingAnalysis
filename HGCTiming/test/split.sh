#!/bin/bash
input=text
output=40mu
mode=startup_sn2.0
HGCP=/eos/cms/store/group/dpg_hgcal/comm_hgcal/suhokim/${output}
find ${HGCP}/RECO/*${mode}*root > ${input}.list
sed -i "s/\/eos\/cms/root:\/\/cms-xrd-global.cern.ch\//g" ${input}.list
split -l 1 --numeric-suffixes --suffix-length=3 --additional-suffix=.list ${input}.list ${output}_${mode}_
rename ${output}_${mode}_0 ${output}_${mode}_ *
rename ${output}_${mode}_0 ${output}_${mode}_ *
mv ${output}*list lists/
rm ${input}.list
