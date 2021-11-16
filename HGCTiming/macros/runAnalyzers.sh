#!/bin/bash

samples=( \
"Gun50"
)
for rawsample in ${samples[@]}
do
 for num in {1..61}
  do
   sample=${rawsample}"_"${num}
   echo ${sample}
   root -l -b 'lookAtMIP_correctRate.C('\""${sample}"\"','\""${rawsample}"\"')'
  done
done
