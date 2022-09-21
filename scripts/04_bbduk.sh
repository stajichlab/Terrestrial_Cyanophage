#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/bbduk.log
#SBATCH -o logs/bbduk.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel,batch
#SBATCH -J virus_qc

IN=fastq
SAMPFILE=fastq/sra_samples.csv
OUT=fastq/filtered

mkdir $OUT

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read PREFIX 
do
	mkdir $OUT/$PREFIX
	
	/rhome/cassande/bigdata/software/bbmap/bbduk.sh in1=$IN/$PREFIX/$PREFIX'_1.fastq.gz' in2=$IN/$PREFIX/$PREFIX'_2.fastq.gz' out1=$OUT/$PREFIX/$PREFIX'_1_filtered.fastq' out2=$OUT/$PREFIX/$PREFIX'_2_filtered.fastq' ref=/rhome/cassande/bigdata/software/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10 maq=10
		
done
