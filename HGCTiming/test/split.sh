#!/bin/bash
input=text
output=myname
HGCP=/eos/cms/store/group/dpg_hgcal/comm_hgcal/suhokim
find ${HGCP}/RECO > text.list
sed -i "s/\/eos\/cms/root:\/\/cms-xrd-global.cern.ch\//g" ${input}.list
split -l 10 --numeric-suffixes --suffix-length=2 --additional-suffix=.list ${input}.list ${output}
mv ${output}*list lists/
