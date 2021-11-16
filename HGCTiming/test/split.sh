#!/bin/bash
input=text
output=Gun50
HGCP=/eos/cms/store/group/dpg_hgcal/comm_hgcal/suhokim/${output}
find ${HGCP}/RECO/*root > ${input}.list
sed -i "s/\/eos\/cms/root:\/\/cms-xrd-global.cern.ch\//g" ${input}.list
split -l 5 --numeric-suffixes --suffix-length=2 --additional-suffix=.list ${input}.list ${output}_
rename ${output}_0 ${output}_ *
mv ${output}*list lists/
rm ${input}.list
