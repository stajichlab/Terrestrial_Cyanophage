#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/vcontact2.prep.log
#SBATCH -o logs/vcontact2.prep.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel
#SBATCH -J vcontact2



DIR=vOTU
PREFIX=TC_vibrant_checkv_combined_clustered
OUTFILE=TC_vibrant_checkv_combined_clustered_prot.fa


module load centos
centos.sh

#vcontact2 on vOTUs

conda activate vContact2

prodigal -i $DIR/$PREFIX.fa -o $DIR/$PREFIX.coords.gbk -a $DIR/$OUTFILE -p meta

vcontact2_gene2genome -p $DIR/$OUTFILE -o $DIR/$PREFIX.gene2genome.csv -s Prodigal-FAA

