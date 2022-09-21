#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -e logs/01_checkv.log
#SBATCH -o logs/01_checkv.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel,batch
#SBATCH -J checkv

CPU=24
IN=vibrant_results/VIBRANT_virsorter1_kbase_all/VIBRANT_phages_virsorter1_kbase_all/virsorter1_kbase_all.phages_combined.fna
SAMPFILE=samples.csv
OUT=vibrant_results/checkv_results

mkdir $OUT

#export CHECKVDB=/rhome/cassande/bigdata/software/checkv-db-v1.0

CHECKVDB=/rhome/cassande/bigdata/software/checkv-db-v1.0


module unload miniconda2
module load miniconda3

conda activate checkv

#checkv to trim
checkv end_to_end $IN $OUT/TC_vibrant_checkv -t $CPU -d $CHECKVDB

