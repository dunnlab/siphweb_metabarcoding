#!/bin/bash
#SBATCH --job-name=silvaTa
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=0-00:50:00


module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1

folders=("fullpipelineRUN1" "fullpipelineRUN3" "fullpipelineRUN4" "fullpipelineRUN5")
barcodes=("152" "166" "272" "179" "261" "134")

for dir in ${folders[@]}; do
  echo $barcodes
  runid=${dir: -4}
  pipeline_dir="project/scratch_backup/"$dir
  for barcode in ${barcodes[@]}; do
    exportpath="$pipeline_dir"/Q2_"$barcode"/export

    sed 's/;/\t/g' "SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy.txt > "SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy-table.tsv
    echo "SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy-table.tsv
    pathtoAssignments_custom="SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy-table.tsv
    pathtoFeatures="$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout-sorted.tsv

    #Create table with OTU IDs and frequency in each sample
    Rscript --vanilla SPLUSDB_results/OTUwrange_metaxa_splus.R $runid $barcode $pathtoAssignments_custom $pathtoFeatures "SPLUSDB_results"

  done

done

runid="RUN0"
dir="RUN0_distal"
pipeline_dir="project/scratch_backup/"$dir
for barcode in ${barcodes[@]}; do
    exportpath="$pipeline_dir"/Q2_"$barcode"/export

    sed 's/;/\t/g' "SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy.txt > "SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy-table.tsv
    echo "SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy-table.tsv
    pathtoAssignments_custom="SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy-table.tsv
    pathtoFeatures="$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout-sorted.tsv

    #Create table with OTU IDs and frequency in each sample
    Rscript --vanilla SPLUSDB_results/OTUwrange_metaxa_splus.R $runid $barcode $pathtoAssignments_custom $pathtoFeatures "SPLUSDB_results"
done

runid="RUN2"
dir="fullpipelineRUN2"
pipeline_dir="project/scratch_backup/"$dir
barcodes=("152" "166" "272" "17nine" "261" "134")
for barcode in ${barcodes[@]}; do
    exportpath="$pipeline_dir"/Q2_"$barcode"/export

    sed 's/;/\t/g' "SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy.txt > "SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy-table.tsv
    echo "SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy-table.tsv
    pathtoAssignments_custom="SPLUSDB_results/METAXA2_plus_$runid$barcode".taxonomy-table.tsv
    pathtoFeatures="$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout-sorted.tsv

    #Create table with OTU IDs and frequency in each sample
    Rscript --vanilla SPLUSDB_results/OTUwrange_metaxa_splus.R $runid $barcode $pathtoAssignments_custom $pathtoFeatures "SPLUSDB_results"
done