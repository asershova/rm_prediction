#!/bin/bash

for i in `cat $1`; do mkdir ../ab/acinetobacter_genomes/${i}_rm_pred_specific; done
echo "mkdir done"
 
for i in `cat $1`; do python3 comparison_domain_types.py ../ab/acinetobacter_genomes/${i}_rm_prediction/rps_blast_results/ ../ab/acinetobacter_genomes/${i}_rm_pred_specific/ family 0; done 
 
echo "domain_comparison done"

for i in `cat $1`; do python3 rm_predictor_types.py ../ab/acinetobacter_genomes/${i}_rm_pred_specific/ ../ab/acinetobacter_genomes/${i}_rm_pred_specific/domain_comparison/${i}* ../ab/acinetobacter_genomes/${i}_rm_prediction/blast_best_hits/${i}* family 0 thresholds; done 

echo "rm_predictor_types done"

for i in `cat $1`; do python3 add_RE_from_blast.py ../ab/acinetobacter_genomes/${i}_rm_pred_specific/predicted_rm_proteins/${i}* ../ab/acinetobacter_genomes/ab_refseq_5_21_genomic_gff/${i}*  ../ab/acinetobacter_genomes/${i}_rm_prediction/blast_best_hits/${i}* ../ab/acinetobacter_genomes/${i}_rm_pred_specific/; done 

echo "add_RE_from_blast done"

for i in `cat $1`; do python3 gene_distance_gff.py ../ab/acinetobacter_genomes/${i}_rm_pred_specific/rm_system_with_RE/${i}* ../ab/acinetobacter_genomes/ab_refseq_5_21_genomic_gff/${i}* ../ab/acinetobacter_genomes/${i}_rm_pred_specific/; done 

echo "gene_distance_gff done"
