from sys import argv
import os
import subprocess

def create_dir(working_dir):
    try:
        outpath = working_dir+"rps_blast_results/" #output directory
        os.mkdir(outpath)
        print ("Path is created:"+outpath)
    except Exception:
        print ("Path exists:"+outpath)
    return outpath

def rps_blast(fasta, outdir):
    db = "/mnt/md0/anna/seq_db/blastdb/cdd_profiles/Cdd"
    rpsbproc_d = "/mnt/md0/anna/distrib/blast/ncbi-blast-2.10.1+-src/c++/ReleaseMT/bin/data" 
    fasta_name = fasta.split("/")[-1]
    evalue = 0.01
    subprocess.run(["rpsblast", "-db", db, "-query", fasta, "-outfmt", "11", "-evalue", str(evalue), \
        "-num_threads", "10", "-out", f'{outdir}/{fasta_name}.cdd.asn'])
    subprocess.run(["rpsbproc", "-i", f'{outdir}/{fasta_name}.cdd.asn', "-o", f'{outdir}/{fasta_name}.cdd.asn.dom.out',\
        "-e", "0.01", "-m", "std", "-t", "doms", "-d", rpsbproc_d])
def main(working_dir, faa):
    outpath = create_dir(working_dir)
    print ("See results in "+outpath)
    rps_blast(faa, outpath)
    print ("rpsblast results are ready")

if __name__ == "__main__":
    fasta = argv[1]
    working_dir = argv[2]
    main(working_dir, fasta)


#rpsblast -db ~/seq_db/blastdb/cdd_profiles/Cdd -query gold_standard.txt.fasta -out gold_standard.cdd.table -evalue 0.01 -outfmt 11 -num_threads 10
#rpsbproc -i blast_rebase.cdd.asn -o blast_rebase.cdd.asn.dom.out -e 0.01 -m std -t doms -d ~/distrib/blast/ncbi-blast-2.10.1+-src/c++/ReleaseMT/bin/data


