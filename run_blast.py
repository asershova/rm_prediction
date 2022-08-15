from sys import argv
import os
import subprocess

def run_blastp(fasta, outpath):
    db = "rebase_gold_standard/gold_standard.txt.fasta"
    #db = "rebase_007/rebase_007"
    evalue = 0.01
    fasta_name = fasta.split("/")[-1]
    subprocess.run(["blastp", "-query", fasta, "-db", db, "-outfmt", "6", "-evalue", str(evalue), \
        "-num_threads", "10", "-out", f'{outpath}{fasta_name}.{str(evalue)}.blast'])
def create_dir(working_dir):    
    try:
        outpath = working_dir+"blastp_results/" #output directory
        os.mkdir(outpath)
        print ("Path is created:"+outpath)
    except Exception:
        print ("Path exists:"+outpath)
    return outpath
def main(working_dir, faa):
    outpath = create_dir(working_dir)
    print ("See results in "+outpath)
    run_blastp(faa, outpath)
    print ("blast results are ready")

if __name__ == "__main__":
    fasta = argv[1]
    working_dir = argv[2]
    main(working_dir, fasta)

#fasta = argv[1] #fasta file (faa)
#db = argv[2]
#fasta_name = fasta.split("/")[-1]
#blastp -query $i -db /mnt/md0/anna/rm_systems/rebase_gold_standard/gold_standard.txt.fasta -outfmt 6 -evalue 0.01 -num_threads 10 -out blast_result/${i}.0.01.blast
#subprocess.run(["blastp", "-query", fasta, "-db", db, "-outfmt", "6", "-evalue", "0.01", \
#        "-num_threads", "10", "-out", f'{outpath}/{fasta_name}.0.01.blast'])

