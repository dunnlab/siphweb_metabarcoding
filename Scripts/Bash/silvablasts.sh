#! /bin/bash

blastn -query 18S272insert_pruned.fasta -db silva_metazoa_db -max_target_seqs 1 -out silvablast_272.out
echo 272_done!
blastn -query 18S261insert.fasta -db silva_metazoa_db -max_target_seqs 1 -out silvablast_261.out
echo 261_done!
blastn -query 18S179insert.fasta -db silva_metazoa_db -max_target_seqs 1 -out silvablast_179.out
echo 179_done!
blastn -query 18S166insert_pruned.fasta -db silva_metazoa_db -max_target_seqs 1 -out silvablast_166.out
echo 166_done!
blastn -query 18S153insert_pruned.fasta -db silva_metazoa_db -max_target_seqs 1 -out silvablast_153.out
echo 153_done!
blastn -query 18S134insert_pruned.fasta -db silva_metazoa_db -max_target_seqs 1 -out silvablast_134.out
echo 134_done!
