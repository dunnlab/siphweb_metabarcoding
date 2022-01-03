#!/bin/bash

folders=("singlendRUN0" "fullpipelineRUN1" "fullpipelineRUN2" "fullpipelineRUN3" "fullpipelineRUN4" "fullpipelineRUN5")

for dir in ${folders[@]}; do

  runid=${dir: -4}
  pipeline_dir="project/scratch_backup/"$dir
  barcodes=("152" "166" "272" "179" "261" "134")
  if [[ $runid == "RUN2" ]]
  then
  barcodes=("152" "166" "272" "17nine" "261" "134")
  fi
  for barcode in ${barcodes[@]}; do

    exportpath="$pipeline_dir"/Q2_"$barcode"/export

    cat $exportpath/$barcode"-sequences-str.fasta" >> allsequences.fasta

  done

done