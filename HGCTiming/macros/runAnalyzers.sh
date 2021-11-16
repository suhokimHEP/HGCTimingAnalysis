#!/bin/bash

samples=( \
##"singleMuon_startup0"   \
"singleMuon_eol0"   \
#"singleMuon_eol1"   \
#"singleMuon_eol2"   \
#"singleMuon_eol3"   \
#"singleMuon_eol4"   \
#"singleMuon_eol5"   \
#"singleMuon_eol6"   \
#"singleMuon_eol7"   \
#"singleMuon_eol8"   \
#"singleMuon_eol9"   \
#"singleMuon_newGun"   \
#"singleMuon_newGun2"   \
#"singleMuon_newGun3"   \
)
for sample in ${samples[@]}
do
  echo ${sample}
  root -l -b 'lookAtMIP_correctRate.C('\""${sample}"\"')'

done
