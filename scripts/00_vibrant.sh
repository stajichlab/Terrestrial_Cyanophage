#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -e logs/vibrant.log
#SBATCH -o logs/vibrant.log
#SBATCH --mem 100G #memory per node in Gb
#SBATCH -p intel,batch
#SBATCH -J vibrant


module unload miniconda2
module load miniconda3

source activate vibrant-1.2.1

#download database
download-db.sh

#ALLVIR=vOTU/TC_checkv_combined.fa
ALLVIR=data/virsorter1_kbase_all.fa

#run vibrant - we have lots of prophage - here will will see if it can pull them out of the contigs that VirSorter/CheckV did not 
VIBRANT_run.py -i $ALLVIR -folder vibrant_results -virome

