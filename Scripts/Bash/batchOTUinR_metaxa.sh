#!/bin/bash
#SBATCH --job-name=assignmentR
#SBATCH --partition=general
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=0-04:00:00

runid="RUN1"
pipeline_dir="fullpipelineRUN1"
barcodes=("152" "166" "272" "179" "261" "134")
manifest="$pipeline_dir"/manifest
metadata="$pipeline_dir"/metadata.tsv

module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1

for barcode in ${barcodes[@]}; do

  sed 's/;/\t/g' "$pipeline_dir"/Q2_"$barcode"/"METAXA2_$runid$barcode".taxonomy.txt > "$pipeline_dir"/Q2_"$barcode"/"METAXA2_$runid$barcode".taxonomy-table.tsv
  sed 's/;/\t/g' "$pipeline_dir"/Q2_"$barcode"/"METAXA2_default_$runid$barcode".taxonomy.txt > "$pipeline_dir"/Q2_"$barcode"/"METAXA2_default_$runid$barcode".taxonomy-table.tsv

  pathtoAssignments="$pipeline_dir"/Q2_"$barcode"/"METAXA2_$runid$barcode".taxonomy-table.tsv
  pathtoAssignments_def="$pipeline_dir"/Q2_"$barcode"/"METAXA2_default_$runid$barcode".taxonomy-table.tsv
  pathtoFeatures="$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout-sorted.tsv

  #Create table with OTU IDs and frequency in each sample
  Rscript --vanilla "$pipeline_dir"/OTUwrange_metaxa.R $runid $barcode $pathtoAssignments $pathtoFeatures "$pipeline_dir/"
  Rscript --vanilla "$pipeline_dir"/OTUwrange_metaxa_default.R $runid $barcode $pathtoAssignments_def $pathtoFeatures "$pipeline_dir/"

done