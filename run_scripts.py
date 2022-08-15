#!/usr/bin/python3
import os, sys, argparse, subprocess
import run_blast, get_blast_best_hit, make_fasta_from_best, run_rps_blast,\
        comparison_domain_types, rm_predictor_types, gene_distance_gff

parser = argparse.ArgumentParser(description='Run R-M system protein prediction')
parser.add_argument('genome_faa', help="The path to the genome faa file")
parser.add_argument('project_dir', help="The name of the project result directory")
parser.add_argument('genome_annotation', help = "The path to the genome annotation file, gff or feature_table format")
parser.add_argument('mode', type=str, choices=["family", "superfamily"],
                    help="compare domains by family or superfamily")
parser.add_argument('m', type=int, choices=[0, 1],
                    help="1 to add several domains in gold_standard data")
parser.add_argument('-gff', action="store_true", help="If an annotation file is in gff format")

args = parser.parse_args()
try:
    os.mkdir(args.project_dir)
    print ("Path is created:"+args.project_dir)
except Exception:
    print ("Path exists:"+args.project_dir)

faa_name = args.genome_faa.split("/")[-1]

blast_outdir = run_blast.create_dir(args.project_dir)
run_blast.run_blastp(args.genome_faa, blast_outdir)
print ("Blast calculated")
best_blast_outdir = get_blast_best_hit.create_dir(args.project_dir)
get_blast_best_hit.best_blast_hits(blast_outdir, best_blast_outdir)
print ("Best blast hits are found")
best_blast_fasta = make_fasta_from_best.create_dir(args.project_dir)
make_fasta_from_best.best_blast_hits_fasta(args.genome_faa, best_blast_outdir, best_blast_fasta)
print ("Best blast hits fasta are obtained")

rps_blast_dir = run_rps_blast.create_dir(args.project_dir)
full_best_faa = best_blast_fasta+faa_name+".0.01.blast.best_hits.fasta"
run_rps_blast.rps_blast(full_best_faa, rps_blast_dir)

domain_comparison_dir = comparison_domain_types.create_dir(args.project_dir)
comparison_domain_types.compare_domains(rps_blast_dir, args.mode, args.m, domain_comparison_dir)

predict_rm_type_dir = rm_predictor_types.create_dir(args.project_dir)
cdd_file = domain_comparison_dir+faa_name+".0.01.blast.best_hits.fasta.cdd.asn.dom.out"+"."+str(args.m)+"."+args.mode
blast_file = best_blast_outdir+faa_name+".0.01.blast.best_hits"
rm_predictor_types.write_results(cdd_file, blast_file, str(args.m), args.mode, predict_rm_type_dir)
if args.gff:
    rm_anno_dir = gene_distance_gff.create_dir(args.project_dir)
    for rm_proteins in os.scandir(predict_rm_type_dir):
        rm_file = predict_rm_type_dir+rm_proteins.name
        gene_distance_gff.rm_annotate(rm_file,args.genome_annotation, rm_anno_dir)
print (f'Done\nThe R-M protein prediction list is in the {rm_anno_dir}\n')

#python3 rm_predictor_types.py test_data/cdd_result/GCF_000196795.1_ASM19679v1_protein.faa.0.01.blast.best_hits.cdd.asn.dom.out.family test_data/blast_best_hits/GCF_000196795.1_ASM19679v1_protein.faa.0.01.blast.best_hits
#python3 comparison_domain.types.py test_data/cdd_result/GCF_000196795.1_ASM19679v1_protein.faa.0.01.blast.best_hits.cdd.asn.dom.out test_data/ family 0
#python3 run_rps_blast.py test_data/seqret_results/GCF_000196795.1_ASM19679v1_protein.faa.0.01.blast.best_hits.fasta /mnt/md0/anna/seq_db/blastdb/cdd_profiles/Cdd test_data/cdd_python



#python3 run_blast.py test_data/proteins_faa/GCF_000196795.1_ASM19679v1_protein.faa /mnt/md0/anna/rm_systems/rebase_gold_standard/gold_standard.txt.fasta /mnt/md0/anna/rm_systems/scripts/test_data/blast_python

#working directory
#directory with files
#file list

