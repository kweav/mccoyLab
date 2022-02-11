#!/bin/bash

# ./run_wgets.sh nl_sumstat_wgets.txt

cat $1 | while read line 
do
   $line
done