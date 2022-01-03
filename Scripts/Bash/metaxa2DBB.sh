#!/bin/bash
#SBATCH --job-name=metaxaDBB
#SBATCH --partition=general
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=3-00:00:00

module load BLAST+/2.8.1-foss-2016b-Python-2.7.13
module load USEARCH/6.1.544
module load HMMER/3.1b2-foss-2016a

metaxa2_dbb -i metaxa2_seqs_zooplanktondb_upper.fasta -o zooplanktonMETAXAdb -t metaxa2_taxonomy_zooplanktondb_clean.fasta -r 1134784017 --cpu 2 --plus T
