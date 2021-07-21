#!/bin/bash

IDEAS=/project/vision/Data/IDEAS/ideasVisionV20p8Seg
BED=.bed
for MATCH in 1 2 3 4 5 6 7 8 9 10
do
  cut -f1,2,3 S3V2_IDEAS_mm10_ccre2.cCRE.M.notall0.N.bed_matched.${MATCH}.bed > random_${MATCH}.bed
  date; time bedtools intersect -a random_${MATCH}.bed -b ${IDEAS}*.bed -wa -wb -filenames > random_${MATCH}_int_all.bed
  date; time bedtools intersect -a random_${MATCH}.bed -b ${IDEAS}*.bed -v > random_${MATCH}_v_all.bed
  for CT in B Cd4 Cd8 Cfue Cfum Clp Cmp Er4 Eryad Eryfl G1e Gmp Hpc7 Imk Mk Lsk Mep Mon Neu Nk
  do
    awk -v var=$IDEAS$CT$BED '{if ($4 == var) {print}}' random_${MATCH}_int_all.bed > random_${MATCH}_int_$CT$BED
    awk '{print $1, $2, $3, $6"_"$7"_"$8}' OFS='\t' random_${MATCH}_int_$CT$BED > random_${MATCH}_istate_$CT$BED
    bedtools groupby -i random_${MATCH}_istate_$CT$BED -g 1,2,3 -c 4 -o collapse > random_${MATCH}_collapse_$CT$BED
  done
done
