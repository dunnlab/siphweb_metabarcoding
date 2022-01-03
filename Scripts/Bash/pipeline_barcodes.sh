#!/bin/bash
#SBATCH --job-name=Q2_pipe  
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=7-00:00:00

### USAGE ###

# Make sure the manifest files and metadata files are formatted adequately
# mkdir pipeline_directory
# Have a copy of the manifest and metadata files in the pipeline_directory
# Manifests paths should be $PWD/$bb_output/filename for each filename
# Metadata_barcode.tsv for QIIME must have unpaired sample IDs
# OTUwrangle.R script should be in the pipeline_directory

#############

#Load packages
module load FastQC/0.11.5-Java-1.8.0_121
module load BBMap/36.62-foss-2016b-Java-1.8.0_121

echo "Modules and environment loaded!"

#Declare temporary variables
runid="RUN5"
pipeline_dir="fullpipelineRUN5"
#pathtoTARGZs="/SAY/standard2/cwd7-CC0522-FASEEB/data/sequences/illumina/SIPHWEB_MiSeq_Lane_4" must cp *.gz from /SAY/ folder
pathtoTARGZs=$pipeline_dir"/targzs"
pathtoFASTQs=$pipeline_dir"/fastq_files"
pathtoHTMLS=$pipeline_dir"/fastqc_HTMLs"
barcodes=("152" "166" "272" "179" "261" "134")
queries=("^TGACGGAAGGGCACCACCAG|^TCCACCAACTAAGAACGGCC|^CTGGTGGTGCCCTTCCGTCA|^GGCCGTTCTTAGTTGGTGGA" "^AACGGCTACCACATCCAAGG|^CACCAGACTTGCCCTCCAAT|^CCTTGGATGTGGTAGCCGTT|^ATTGGAGGGCAAGTCTGGTG" "^AAACGATGCCGACTAGCGAT|^TCCACCAACTAAGAACGGCC|^ATCGCTAGTCGGCATCGTTT|^GGCCGTTCTTAGTTGGTGGA" "^GGCCGTTCTTAGTTGGTGGA|^TGCGGCCCAGAACATCTAAG|^TCCACCAACTAAGAACGGCC|^CTTAGATGTTCTGGGCCGCA" "^AACAGGTCTGTGATGCCCTT|^TGTGTACAAAGGGCAGGGAC|^AAGGGCATCACAGACCTGTT|^GTCCCTGCCCTTTGTACACA" "^CTTTGTACACACCGCCCGTC|^CCTTGTTACGACTTTTACTTCCTCT|^GACGGGCGGTGTGTACAAAG|^AGAGGAAGTAAAAGTCGTAACAAGG")
grep_output=$pipeline_dir/"grep_output"
bb_output=$pipeline_dir"/bb_output"
bb_orphans=$pipeline_dir"/bb_orphans"
cutadapt_o=$pipeline_dir"/cutadapt_o"
manifest="$pipeline_dir"/manifest

mkdir $pipeline_dir
mkdir $pathtoFASTQs
mkdir $pathtoHTMLS
mkdir $grep_output
mkdir $bb_output
mkdir $bb_orphans
mkdir $cutadapt_o

echo "Directories and paths declared!"

#Extract GZ balls
gzfiles=($(ls "$pathtoTARGZs"))
for f in ${gzfiles[@]}; do gunzip -c "$pathtoTARGZs"/"$f" > "$pathtoFASTQs"/"${f%.*}" ; done

echo ".gz files extracted!"

#Quality HTMLS
fastqc "$pathtoFASTQs"/* -o $pathtoHTMLS

echo "Fastqc quality html reports generated!"

#Segragate by Barcode/Amplicon
filenames=($(ls "$pathtoFASTQs")) #for the grep by amplicon stage
for file in ${filenames[@]}; do
	for bc in 0 1 2 3 4 5; do
		grep -E "${queries[$bc]}" -C 1 -A 2 "$pathtoFASTQs"/"$file"  | sed 's\^--$\\g' | sed '/^$/d' >> "$grep_output"/"${barcodes[$bc]}"_"$file"
	done
done

echo "Reads sorted by PCR primers into the 6 different amplicons!"

#Remove PCR primers
files=($(ls "$grep_output")) #for the cutadapt PCR primer removal
for file in ${files[@]}; do
	cutadapt -m 50 -b TGACGGAAGGGCACCACCAG -b TCCACCAACTAAGAACGGCC -b CTGGTGGTGCCCTTCCGTCA -b GGCCGTTCTTAGTTGGTGGA  -b AACGGCTACCACATCCAAGG -b CACCAGACTTGCCCTCCAAT -b CCTTGGATGTGGTAGCCGTT -b ATTGGAGGGCAAGTCTGGTG  -b AAACGATGCCGACTAGCGAT -b TCCACCAACTAAGAACGGCC -b ATCGCTAGTCGGCATCGTTT -b GGCCGTTCTTAGTTGGTGGA  -b GGCCGTTCTTAGTTGGTGGA -b TGCGGCCCAGAACATCTAAG -b TCCACCAACTAAGAACGGCC -b CTTAGATGTTCTGGGCCGCA  -b AACAGGTCTGTGATGCCCTT -b TGTGTACAAAGGGCAGGGAC -b AAGGGCATCACAGACCTGTT -b GTCCCTGCCCTTTGTACACA  -b CTTTGTACACACCGCCCGTC -b CCTTGTTACGACTTTTACTTCCTCT -b GACGGGCGGTGTGTACAAAG -b AGAGGAAGTAAAAGTCGTAACAAGG -o "$cutadapt_o"/"$file" "$grep_output"/"$file"
done

echo "PCR primer sequences removed from the 5' end!"

#Make general manifest
names=($(ls $cutadapt_o | grep -oE "[0-9]\w+_L001"))
let i=1
for name in ${names[@]}; do 
 if (( $i % 2 )); then
  echo $name,"\$PWD/\$bb_output/"$name"_R1_001.fastq,forward" >> $manifest
 else
  echo $name,"\$PWD/\$bb_output/"$name"_R2_001.fastq,reverse" >> $manifest
 fi
  let i=i+1
done

#Repair non-matching read pairs due to differential sequencing errors in R1 and R2
for barcode in ${barcodes[@]}; do
  samplenames=($(grep -oE "\w+_L001" $manifest))
  for file in ${samplenames[@]}; do
    repair.sh in="$cutadapt_o"/"$file"_R1_001.fastq in2="$cutadapt_o"/"$file"_R2_001.fastq out1="$bb_output"/"$file"_R1_001.fastq out2="$bb_output"/"$file"_R2_001.fastq outsingle="$bb_orphans"/"$barcode"_"$file"_orphans.fastq
  done
done

echo "R1 and R2 reads paired up and ready to import!"

# *** IMPORTANT *** MANIFEST MUST CONTAIN THE PATHS TO BB_OUTPUT!
module purge
source activate qiime2-2018.11
module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1
module load dSQ

for barcode in ${barcodes[@]}; do


  echo Starting QIIME2 for barcode_"$barcode"

  manifest="$pipeline_dir"/manifest_"$barcode"
  metadata="$pipeline_dir"/metadata_"$barcode".tsv
  exportpath="$pipeline_dir"/Q2_"$barcode"/export

  #Make manifest files
  names=($(ls $bb_output | grep -oE ^$barcode"\w+_L001"))
  echo "sample-id,absolute-filepath,direction" > $manifest
  let i=1
  for name in ${names[@]}; do 
     if (( $i % 2 )); then
          echo $name",/home/ad2258/scratch60/"$bb_output"/"$name"_R1_001.fastq,forward" >> $manifest
     else
          echo $name",/home/ad2258/scratch60/"$bb_output"/"$name"_R2_001.fastq,reverse" >> $manifest
     fi
    let i=i+1
  done

  #Make metadata files
  uniqnames=($(ls $bb_output | grep -oE ^$barcode"\w+_L001" | uniq))
  echo -e "SampleID\\tBarcodeSequence\\tLinkerPrimerSequence\\tlane\\ttemplate" > $metadata
  for uname in ${uniqnames[@]}; do 
    echo -e $uname"\\tGCACCACCAA\\tCAGCACGGAG\\t1\\t"$barcode >> $metadata
  done

  #Remove pre-existing dada outputs
  rm -rf "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"_dada2_out
  rm -rf "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results

  mkdir "$pipeline_dir"/Q2_"$barcode"
  mkdir $exportpath

  #Load data into QIIME2
  qiime metadata tabulate --m-input-file $metadata --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"_metadata.qzv

  qiime tools import \
     --type 'SampleData[PairedEndSequencesWithQuality]' \
     --input-path $manifest \
     --output-path "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qza \
     --input-format PairedEndFastqManifestPhred33

  qiime demux summarize \
    --i-data "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qza \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qzv

  echo "Paired reads imported into QIIME2!"

  #DADA2 Denoise and dereplicate
  #let trunclenF=$barcode-($barcode/5)
  #let trunclenR=$barcode-($barcode/5)
  qiime dada2 denoise-paired \
          --i-demultiplexed-seqs "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qza \
          --p-trunc-len-f 0 \
          --p-trunc-len-r 0 \
          --p-trunc-q 28 \
          --o-representative-sequences "$pipeline_dir"/Q2_"$barcode"/rep-"$runid""$barcode".qza \
          --o-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-table.qza \
          --p-n-threads 24 \
          --output-dir "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"_dada2_out

  qiime feature-table summarize \
    --i-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-table.qza \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-table.qzv \
    --m-sample-metadata-file $metadata

  qiime feature-table tabulate-seqs \
    --i-data "$pipeline_dir"/Q2_"$barcode"/rep-"$runid""$barcode".qza \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/rep-"$runid""$barcode".qzv

  qiime vsearch cluster-features-de-novo \
    --i-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-table.qza \
    --i-sequences "$pipeline_dir"/Q2_"$barcode"/rep-"$runid""$barcode".qza \
    --p-perc-identity 0.95 \
    --o-clustered-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-095-table.qza \
    --o-clustered-sequences "$pipeline_dir"/Q2_"$barcode"/rep-095-"$runid""$barcode".qza

  qiime tools export \
    --input-path "$pipeline_dir"/Q2_"$barcode"/rep-095-"$runid""$barcode".qza \
    --output-path $exportpath

  qiime tools export \
  --input-path "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-095-table.qza \
  --output-path $exportpath

  echo "Data denoised, merged, and demultiplexed with DADA2 and clustered with VSEARCH!"

  #Multiple sequence alignment using Mafft
   qiime alignment mafft \
    --i-sequences "$pipeline_dir"/Q2_"$barcode"/rep-095-"$runid""$barcode".qza \
    --o-alignment "$pipeline_dir"/Q2_"$barcode"/aligned-rep-095-"$runid""$barcode".qza

  #Mask the alignment to remove positions that are highly variable.
  qiime alignment mask \
    --i-alignment "$pipeline_dir"/Q2_"$barcode"/aligned-rep-095-"$runid""$barcode".qza \
    --o-masked-alignment "$pipeline_dir"/Q2_"$barcode"/masked-aligned-rep-095-"$runid""$barcode".qza

  echo "Sequences aligned!"

  #Create the tree using the Fasttree program
  qiime phylogeny fasttree \
    --i-alignment "$pipeline_dir"/Q2_"$barcode"/masked-aligned-rep-095-"$runid""$barcode".qza \
    --o-tree "$pipeline_dir"/Q2_"$barcode"/unrooted-"$runid""$barcode".qza

  #Root the tree using the longest root
  qiime phylogeny midpoint-root \
    --i-tree "$pipeline_dir"/Q2_"$barcode"/unrooted-"$runid""$barcode".qza \
    --o-rooted-tree "$pipeline_dir"/Q2_"$barcode"/rooted-"$runid""$barcode".qza

  echo "Phylogenetic tree generated!"

  #Alpha rarefaction
  qiime diversity alpha-rarefaction \
    --i-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-095-table.qza \
    --i-phylogeny "$pipeline_dir"/Q2_"$barcode"/rooted-"$runid""$barcode".qza \
    --p-max-depth 4000 \
    --m-metadata-file $metadata \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-rarefaction.qzv

  #Alpha diversity metrics
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "$pipeline_dir"/Q2_"$barcode"/rooted-"$runid""$barcode".qza \
    --i-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-095-table.qza \
    --p-sampling-depth 1100 \
    --m-metadata-file $metadata \
    --output-dir "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results

  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/faith_pd_vector.qza \
    --m-metadata-file $metadata \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/"$runid""$barcode"faith-pd-group-significance.qzv

  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/evenness_vector.qza \
    --m-metadata-file $metadata \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/"$runid""$barcode"evenness-group-significance.qzv

  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/shannon_vector.qza \
    --m-metadata-file $metadata \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/"$runid""$barcode"shannon_group-significance.qzv

  echo Diversity metrics computed and plotted. End of the QIIME process for barcode "$barcode"

  #Convert feature table to TSV
  biom convert --to-tsv -i "$exportpath"/feature-table.biom -o "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout.tsv
  #Remove the header 2 rows and Sort feature table by abundance row sum
  tail -n +2 "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout.tsv | awk '{for(x=2;x<=NF;x++)s+=$x;print s,$0;s=0}' | sort -nr |sed 's/^\S* //' > "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout-1.tsv
  #Get sample name header line
  sed '2q;d' "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout.tsv | sed 's/^#OTU //' > "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-header.txt
  #Concatenate header line with sorted rows
  cat "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-header.txt "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout-1.tsv > "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout-sorted.tsv
  #Slice the top 100 features
  head -101 "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout-sorted.tsv > "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-feature-095-table.tsv
  cat "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-feature-095-table.tsv | awk '{print $1}' > "$pipeline_dir"/Q2_"$barcode"/top100"$barcode".txt
  tail -n +1 "$pipeline_dir"/Q2_"$barcode"/top100"$barcode".txt > "$pipeline_dir"/Q2_"$barcode"/top100"$barcode"_neat.txt
  awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' "$exportpath"/dna-sequences.fasta > "$exportpath"/"$barcode"-sequences-str.fasta
  grep -f "$pipeline_dir"/Q2_"$barcode"/top100"$barcode"_neat.txt -A 1 "$exportpath"/"$barcode"-sequences-str.fasta | sed 's\^--$\\g' | sed '/^$/d' >> "$pipeline_dir"/Q2_"$barcode"/top100"$barcode".fasta

  #Split sequences in sets of 10 for SAP
  #mkdir "$pipeline_dir"/Q2_"$barcode"/splitout"$runid""$barcode"
  #split -l 10 "$pipeline_dir"/Q2_"$barcode"/top100"$barcode".fasta "$pipeline_dir"/Q2_"$barcode"/splitout"$runid""$barcode"/splitout
  #splitfiles=($(ls "$pipeline_dir"/Q2_"$barcode"/splitout"$runid""$barcode"))
  #mkdir "$pipeline_dir"/Q2_"$barcode"/splitSAP"$runid""$barcode"
  #rm "$pipeline_dir"/Q2_"$barcode"/SAP"$runid""$barcode"_jobfile.txt
  #for i in ${splitfiles[@]}; do
  #  echo "module load BLAST+/2.8.1-foss-2016b-Python-2.7.13; module load ClustalW2/2.1-foss-2016b; ~/miniconda2/bin/sap --database /home/ad2258/project/scratchfullpipelineRUN1_resources/zooplanktonSAP_db.fasta --project $pipeline_dir/Q2_$barcode/splitSAP$runid$barcode/$i --email alejandro.damianserrano@yale.edu $pipeline_dir/Q2_$barcode/splitout$runid$barcode/$i" | cat >> "$pipeline_dir"/Q2_"$barcode"/SAP"$runid""$barcode"_jobfile.txt
    # --minsignificance 0.9 -s 0.1 --minidentity 0.1 -n 10 --ppcutoff 60 --svg # optional arguments for SAP
  #done
  #dSQ --jobfile "$pipeline_dir"/Q2_"$barcode"/SAP"$runid""$barcode"_jobfile.txt --ntasks 4 --cpus-per-task=4 --mem-per-cpu=4g -t 7-10:00:00 > "$pipeline_dir"/Q2_"$barcode"/runSAP_"$runid""$barcode".sh
  
done

echo QIIME2 finished

## METAXA2 ##

module purge
module load BLAST/2.2.22-Linux_x86_64
module load USEARCH/6.1.544
module load HMMER/3.1b2-foss-2016a

for barcode in ${barcodes[@]}; do
  echo METAXA2 started for barcode "$barcode"!
  #Run Metaxa2 on the whole feature fasta file for each barcode since it is fast
  metaxa2 -i $exportpath/$barcode-sequences-str.fasta -o $pipeline_dir/Q2_$barcode/METAXA2_$runid$barcode -g SSU_SILVA123.1 -T 0,60,60,60,60,65,70,75,85,90,97 --reltax T -R 70
  metaxa2 -i $exportpath/$barcode-sequences-str.fasta -o $pipeline_dir/Q2_$barcode/METAXA2_custom_$runid$barcode -g SILVAPLUS -T 0,60,60,60,60,65,70,75,85,90,97 --reltax T -R 70 
  echo METAXA2 finished for barcode "$barcode"!
done

#Continue to OTUwrange batch script and linked OTU_wrange.R file (in pipeline_dir) for subsequen steps

module purge
module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1

for barcode in ${barcodes[@]}; do

  sed 's/;/\t/g' "$pipeline_dir"/Q2_"$barcode"/"METAXA2_$runid$barcode".taxonomy.txt > "$pipeline_dir"/Q2_"$barcode"/"METAXA2_$runid$barcode".taxonomy-table.tsv
  sed 's/;/\t/g' "$pipeline_dir"/Q2_"$barcode"/"METAXA2_custom_$runid$barcode".taxonomy.txt > "$pipeline_dir"/Q2_"$barcode"/"METAXA2_custom_$runid$barcode".taxonomy-table.tsv

  pathtoAssignments="$pipeline_dir"/Q2_"$barcode"/"METAXA2_$runid$barcode".taxonomy-table.tsv
  pathtoAssignments_def="$pipeline_dir"/Q2_"$barcode"/"METAXA2_custom_$runid$barcode".taxonomy-table.tsv
  pathtoFeatures="$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-biomout-sorted.tsv

  #Create table with OTU IDs and frequency in each sample
  Rscript --vanilla "$pipeline_dir"/OTUwrange_metaxa.R $runid $barcode $pathtoAssignments $pathtoFeatures "$pipeline_dir/"
  Rscript --vanilla "$pipeline_dir"/OTUwrange_metaxa_default.R $runid $barcode $pathtoAssignments_def $pathtoFeatures "$pipeline_dir/"

done
