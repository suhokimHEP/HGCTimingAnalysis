#!/bin/bash

cversions=( \
 "Rand40mu_eol" \
# "Rand40mu_startup_sn2.0" \
)


for cversion in ${cversions[@]}
do
 hadd -f ${cversion}.root gitignore/${cversion}/*${cversion}*.root
done

