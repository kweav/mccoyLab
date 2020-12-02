#!/bin/bash

#Usage: source run_kallisto.sh /project/vision/Data/mm10_genome/allChromff.fa /project/vision/Data/mm10_genome/gencode.vM4.annotation.gtf /project/vision/Data/mm10.len
# should probably redo this to do the wget, unzip, concat all by itself for the allChromff.fa file, annotation, and mm10.len

#get fastq files

dir=$(cd -P -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P) #set the working directory

wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.96.gtf.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
gunzip Mus_musculus.GRCm38.96.gtf.gz
cat mm10.chrom.sizes | sed 's/^chr//' > mm10.chrom.sizes.new #need to match the chromosome names between the gtf and chrom sizes file

source activate sra_tools_env #activate conda environment
mkdir rnaseq_quant rnaseq_fastq rnaseq_fastq/LaraAstiaso2014 rnaseq_fastq/Heuston2018 #start file structure

# # #fastq-dump single-end reads from Lara-Astiaso 2014 PRJNA257656 SRP045264 GSE60101
cd ${dir}/rnaseq_fastq/LaraAstiaso2014
allSamples=(SRX669482 SRX669483 SRX669484 SRX669485 SRX669486 SRX669487 SRX669488 SRX669489 SRX669490 SRX669491 SRX669492 SRX669493 SRX669494 SRX669495 SRX669496 SRX669497)
allCTs=(LTHSC HSC MPP CLP CMP GMP MF Granulocyte MON B CD4 CD8 NK MEP ERYA ERYB)
for INDEX in ${!allSamples[@]}
do
	SAMPLE=${allSamples[$INDEX]}
	CT=${allCTs[$INDEX]}
	mkdir ${dir}/rnaseq_quant/${SAMPLE}_${CT} #make corresponding directory for later kallisto output
	mkdir ${SAMPLE}_${CT}
	cd ${dir}/rnaseq_fastq/LaraAstiaso2014/${SAMPLE}_${CT}
	touch ${SAMPLE}_${CT}.txt
	mapfile -t SRA_array < ${dir}/toSRA_files/${SAMPLE}_toSRA.txt
	for INDEX_R in ${!SRA_array[@]}
	do
		prefetch ${SRA_array[$INDEX_R]} &>> ${SAMPLE}_${CT}.txt && fastq-dump ${SRA_array[$INDEX_R]} &>> ${SAMPLE}_${CT}.txt
	done
	cd ${dir}/rnaseq_fastq/LaraAstiaso2014
done

# #fastq-dump paired-end reads from Heuston et al 2018 PRJNA63471 SRP013703
cd ${dir}/rnaseq_fastq/Heuston2018
allSamples2=(SRX3010287 SRX3010288 SRX3009974 SRX3009975 SRX3010022 SRX3010023 SRX3010290 SRX3010291 SRX3010176 SRX3010177 SRX3009997 SRX3009998 SRX3010228 SRX3010229 SRX3010199 SRX3010200 SRX7517165 SRX7517166 SRX7517167 SRX7517168 SRX3010044 SRX3010045 SRX3010068 SRX3010069)
allCTs2=(LSK LSK CMP CMP GMP GMP MEP MEP CFUE CFUE ERY ERY CFUMK CFUMK IMK IMK MON MON NEU NEU G1E G1E ER4 ER4)
for INDEX2 in ${!allSamples2[@]}
do
	SAMPLE2=${allSamples2[$INDEX2]}
	CT2=${allCTs2[$INDEX2]}
	mkdir ${dir}/rnaseq_quant/${SAMPLE2}_${CT2} #make corresponding directory for later kallisto output
	mkdir ${SAMPLE2}_${CT2}
	cd ${dir}/rnaseq_fastq/Heuston2018/${SAMPLE2}_${CT2}
	touch ${SAMPLE2}_${CT2}.txt
	mapfile -t SRA_array < ${dir}/toSRA_files/${SAMPLE2}_toSRA.txt
	for INDEX_R in ${!SRA_array[@]}
	do
		prefetch ${SRA_array[$INDEX_R]} &>> ${SAMPLE2}_${CT2}.txt && fastq-dump -I --split-files ${SRA_array[$INDEX_R]} &>> ${SAMPLE2}_${CT2}.txt
	done
	cd ${dir}/rnaseq_fastq/Heuston2018
done

conda deactivate
cd $dir

#run kallisto
source activate kallisto_env

# #generate index
kallisto index -i mm10_index_kallisto.idx Mus_musculus.GRCm38.cdna.all.fa

#quantify the single-end reads
for INDEX in ${!allSamples[@]}
do
	SAMPLE=${allSamples[$INDEX]}
	CT=${allCTs[$INDEX]}
	# #find average and stdev of read length
	# python read_len.py ${dir}/average_len_temp.txt ${dir}/stdev_len_temp.txt ${dir}/rnaseq_fastq/LaraAstiaso2014/${SAMPLE}_${CT}/*.fastq
	# AVGLEN="$(cat ${dir}/average_len_temp.txt)"
	AVGLEN=200
	# STDEVLEN="$(cat ${dir}/stdev_len_temp.txt)"
	STDEVLEN=30
	kallisto quant -i mm10_index_kallisto.idx -b 100 --seed 38 -o ${dir}/rnaseq_quant/${SAMPLE}_${CT} ${dir}/rnaseq_fastq/LaraAstiaso2014/${SAMPLE}_${CT}/*.fastq --single -l $AVGLEN -s $STDEVLEN -t 40 --pseudobam --genomebam -g Mus_musculus.GRCm38.96.gtf -c mm10.chrom.sizes.new
	# rm ${dir}/average_len_temp.txt ${dir}/stdev_len_temp.txt

done


#quantify the paired-end reads
for INDEX2 in ${!allSamples2[@]}
do
	SAMPLE2=${allSamples2[$INDEX2]}
	CT2=${allCTs2[$INDEX2]}
	kallisto quant -i mm10_index_kallisto.idx -b 100 --seed 38 -o ${dir}/rnaseq_quant/${SAMPLE2}_${CT2} ${dir}/rnaseq_fastq/Heuston2018/${SAMPLE2}_${CT2}/*.fastq -t 40 --pseudobam --genomebam -g Mus_musculus.GRCm38.96.gtf -c mm10.chrom.sizes.new
done

conda deactivate

cd $dir

source activate basic

date;time python create_tx2gene_df.py Mus_musculus.GRCm38.96.gtf

conda deactivate

source activate r_env

cd ${dir}/rnaseq_quant
date;time Rscript ../gene_level_quant.R &> ../quant_gene_level.txt

conda deactivate

source activate basic

cd $dir
date; time python look_at_reps.py

conda deactivate
