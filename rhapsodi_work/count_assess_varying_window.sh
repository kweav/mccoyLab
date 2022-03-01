#!/bin/bash

#Usage: ./count_assess_varying_window.sh numgams numsnps cov

OVLP=2
COUNTER=0

for WIN in 250 500 750 `seq 1000 500 $2`
do
  for AVGR in 0.6 1 3
  do
    for SEQERR in 0.001 0.005 0.05
    do
      for RSD in 42 357 1848
      do
        let COUNTER+=1
      done
    done
  done
done

echo $COUNTER
