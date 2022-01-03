#!/bin/bash
#SBATCH --job-name=metaxadbb_com
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=4-00:00:00

#METAXA2 database builder

module purge
module load BLAST+/2.8.1-foss-2018b \
            HMMER/3.2.1-foss-2018b \
            MAFFT/7.471-foss-2018b-with-extensions \
            USEARCH/6.1.544

metaxa2_dbb -i combined_sequences.fasta -o COMBINED -g COMBINED -t combined_taxonomy.fasta --cpu 4 --correct_taxonomy T --auto_rep T --plus T
