#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/06_bamm.log
#SBATCH -o logs/06_bamm.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel,batch
#SBATCH -J virus_bamm



OUT=fastq/bbmap
DIR=vOTU


module unload miniconda2
module load miniconda3

conda activate bamm

bamm parse -c $DIR/output_file_tpmean_121421.tsv -b $OUT/*.bam -m 'tpmean'
bamm parse -c $DIR/output_file_count_121421.tsv -b $OUT/*.bam -m 'counts'

gzip $OUT/*bam
gzip $OUT/*bai
