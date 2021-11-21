#!/bin/bash

samples=( \
"40mu_startup_sn2.0"
)
for rawsample in ${samples[@]}
do
 for num in {13..13}
  do
   sample=${rawsample}"_"${num}
   echo ${sample}
   root -l -b 'lookAtMIP_correctRate.C('\""${sample}"\"','\""${rawsample}"\"')'
  done
done
