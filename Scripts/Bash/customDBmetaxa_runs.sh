#!/bin/bash
#SBATCH --job-name=SPLUS 
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=7-00:00:00

module load BLAST/2.2.22-Linux_x86_64
module load USEARCH/6.1.544
module load HMMER/3.1b2-foss-2016a

folders=("singlendRUN0" "fullpipelineRUN1" "fullpipelineRUN2" "fullpipelineRUN3" "fullpipelineRUN4")
barcodes=("152" "166" "272" "179" "261" "134")

for dir in ${folders[@]}; do

  runid=${dir: -4}
  pipeline_dir="project/scratch_backup/"$dir
  if (( $dir == "fullpipelineRUN2" )); then
    barcodes=("152" "166" "272" "17nine" "261" "134")
  else
    barcodes=("152" "166" "272" "179" "261" "134")
  fi

  for barcode in ${barcodes[@]}; do
    
    manifest="$pipeline_dir"/manifest_"$barcode"
    metadata="$pipeline_dir"/metadata_"$barcode".tsv
    exportpath="$pipeline_dir"/Q2_"$barcode"/export

    #Run Metaxa2 on the whole feature fasta file for each barcode since it is fast

    metaxa2 -i $exportpath/$barcode-sequences-str.fasta -o SPLUSDB_results/METAXA2_plus_$runid$barcode -g SILVAPLUS -T 0,60,60,60,60,65,70,75,85,90,97 --reltax T -R 70
    
    echo METAXA2 SILVAPLUS started for barcode "$barcode"!

  done

done