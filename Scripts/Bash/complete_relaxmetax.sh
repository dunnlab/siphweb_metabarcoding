#!/bin/bash
#SBATCH --job-name=MetaxaRUN
#SBATCH --partition=general
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=0-01:00:00

runid="RUN1"
pipeline_dir="fullpipelineRUN1"
barcodes=("152" "166" "272" "179" "261" "134")

module load dSQ

for barcode in ${barcodes[@]}; do

  exportpath="$pipeline_dir"/Q2_"$barcode"/export

  #Make DeadSimple jobfile for each fasta subfile
  rm "$pipeline_dir"/Q2_"$barcode"/METAXA2"$runid""$barcode"_jobfile.txt
  rm "$pipeline_dir"/Q2_"$barcode"/METAXA2_custom_"$runid""$barcode"_jobfile.txt

  echo "module load BLAST/2.2.22-Linux_x86_64; module load USEARCH/6.1.544; module load HMMER/3.1b2-foss-2016a; metaxa2 -i $exportpath/$barcode-sequences-str.fasta -o $pipeline_dir/Q2_$barcode/METAXA2_$runid$barcode" -g SSU_SILVA123.1 -T 0,60,60,60,60,65,70,75,85,90,97 --reltax T -R 70 | cat >> "$pipeline_dir"/Q2_"$barcode"/METAXA2"$runid""$barcode"_jobfile.txt
  echo "module load BLAST/2.2.22-Linux_x86_64; module load USEARCH/6.1.544; module load HMMER/3.1b2-foss-2016a; metaxa2 -i $exportpath/$barcode-sequences-str.fasta -o $pipeline_dir/Q2_$barcode/METAXA2_custom_$runid$barcode" -g COMBINED -T 0,60,60,60,60,65,70,75,85,90,97 --reltax T -R 70 | cat >> "$pipeline_dir"/Q2_"$barcode"/METAXA2_custom_"$runid""$barcode"_jobfile.txt
  dSQ --jobfile "$pipeline_dir"/Q2_"$barcode"/METAXA2"$runid""$barcode"_jobfile.txt --ntasks 4 --cpus-per-task=4 --mem-per-cpu=4g -t 7-10:00:00 > "$pipeline_dir"/Q2_"$barcode"/METAXA2_"$runid""$barcode".sh
  dSQ --jobfile "$pipeline_dir"/Q2_"$barcode"/METAXA2_custom_"$runid""$barcode"_jobfile.txt --ntasks 4 --cpus-per-task=4 --mem-per-cpu=4g -t 7-10:00:00 > "$pipeline_dir"/Q2_"$barcode"/METAXA2_custom_"$runid""$barcode".sh
  sbatch dsq-METAXA2"$runid$barcode"_jobfile-20$(date +%y-%m-%d).sh
  sbatch dsq-METAXA2_custom_"$runid$barcode"_jobfile-20$(date +%y-%m-%d).sh

done
