#!/bin/bash

wget --timestamping 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' -O ~/scr4_rmccoy22/kweave23/resources/hg38toHg19.over.chain.gz
wget --no-check-certificate "https://drive.google.com/uc?export=download&id=1YcrXGJKu8NLk8jMwce11F3ylEjVdI3gZ" -O 'S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.coe_mat.PDmerged.clusterID.txt'

ml anaconda
conda create -c bioconda -n liftover ucsc-liftover -y
conda activate liftover

awk 'BEGIN { OFS = "_" } ; NR==FNR{a[$1, $2, $3]; next} !(($1, $2, $3) in a)' ../IDEAS_sept/S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.hg38.failed.nocomments.bed S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.coe_mat.PDmerged.clusterID.txt > metaclusters_withoutFailedLiftover.txt

for clusternum in {1..25}
do
  head -n 1 metaclusters_withoutFailedLiftover.txt > metacluster${clusternum}_annotation.txt
  awk '{ if ($6 == '${clusternum}') {print}' metaclusters_withoutFailedLiftover.txt >> metacluster${clusternum}_annotation.txt
  awk -v clusternum=$clusternum 'BEGIN { OFS = "\t"} { if ($6 == clusternum) {print $1,$2,$3}}' metaclusters_withoutFailedLiftover.txt > metacluster${clusternum}_annotation.bed
  liftOver metacluster${clusternum}_annotation.bed ~/scr4_rmccoy22/kweave23/resources/hg38toHg19.over.chain.gz metacluster${clusternum}_liftover.bed unMapped_${clusternum}.bed
done

conda deactivate
