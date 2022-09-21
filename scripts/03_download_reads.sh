#!/usr/bin/bash
#SBATCH -p intel,batch -N 1 -n 2 --mem 4gb
#SBATCH -J virus_download
module load aspera

KEY=/rhome/cassande/private.openssh.pub

/bigdata/stajichlab/shared/bin/sra_download.pl --ascp --id $ASPERAKEY sra.txt
