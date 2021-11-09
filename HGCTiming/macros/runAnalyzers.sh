#!/bin/bash

samples=( \
##"singleMuon_startup0.root"   \
"singleMuon_eol0.root"   \
#"singleMuon_eol1.root"   \
#"singleMuon_eol2.root"   \
#"singleMuon_eol3.root"   \
#"singleMuon_eol4.root"   \
#"singleMuon_eol5.root"   \
#"singleMuon_eol6.root"   \
#"singleMuon_eol7.root"   \
#"singleMuon_eol8.root"   \
#"singleMuon_eol9.root"   \
#"singleMuon_newGun.root"   \
#"singleMuon_newGun2.root"   \
#"singleMuon_newGun3.root"   \
)
for sample in ${samples[@]}
do
  echo ${sample}
  ./lookAtMIP_correctRate ${sample}

done
