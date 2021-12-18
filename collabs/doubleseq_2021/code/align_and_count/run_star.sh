#!/bin/bash

for sample in `ls *.fastq.gz | cut -f 1-2 -d '_' | uniq`
do
  echo ${sample}
  STAR --runMode alignReads \
       --runThreadN 20 \
       --genomeDir ~/workzfs-rmccoy22/rmccoy22/create/feb_2020/hg38_index \
       --readFilesCommand zcat \
       --readFilesIn ${sample}_L001_R1_001.fastq.gz,${sample}_L002_R1_001.fastq.gz ${sample}_L001_R2_001.fastq.gz,${sample}_L002_R2_001.fastq.gz \
       --outFileNamePrefix /scratch/groups/rmccoy22/kweave23/create_outcome_pred/data_202111/bam/${sample} \
       --outSAMtype BAM SortedByCoordinate \
       --quantMode GeneCounts
done