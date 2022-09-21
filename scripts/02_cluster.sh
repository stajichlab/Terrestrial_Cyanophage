#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -e logs/02_cdhit.log
#SBATCH -o logs/02_cdhit.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel
#SBATCH -J cdhit

CPU=24
DIR=vOTU

#module load cd-hit/4.8.1  

/rhome/cassande/bigdata/software/cd-hit-v4.8.1-2019-0228/cd-hit-est -i $DIR/TC_vibrant_checkv_combined.fa -o $DIR/TC_vibrant_checkv_combined_clustered.fa -c 0.95 -aS 0.85 -M 0 -d 0 -T $CPU
