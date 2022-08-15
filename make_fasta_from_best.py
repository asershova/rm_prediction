from sys import argv
import subprocess
import sys
import os
import tempfile

def create_dir(working_dir):
    directory = working_dir+"blast_best_hits_fasta/"

    try:
       os.mkdir(directory)
    except:
       pass
    return directory
def best_blast_hits_fasta(protein_faa, best_blast_dir, outdir):
    entries = os.scandir(best_blast_dir)
    for entry in entries:
        besthits_file = open(best_blast_dir+entry.name, "r")#best_hits filename
        name = entry.name
        fp = tempfile.NamedTemporaryFile(mode = 'w', buffering = 1)
        out = outdir+"/"+name + ".fasta"
        try:
           os.unlink(out)
        except Exception:
           pass

        #read line by line, strip and split lines
        for line in besthits_file:
            if line.startswith("#"):
                continue
            sline = line.strip().split()
            protein_id = sline[2]
            fp.write(f'{protein_faa}:{protein_id}\n')
        subprocess.run(["/mnt/md0/kroegerlab/bin/bin/seqret", f'@{fp.name}', "-outseq", out])
def main(protein_faa, best_blast_dir, working_dir):
    outpath = create_dir(working_dir)
    print ("See results in "+outpath)
    best_blast_hits_fasta(protein_dir, best_blast_dir, outpath)
    print ("best hits fasta are ready")

if __name__ == "__main__":
    protein_faa = argv[1]
    best_blast_dir = argv[2]
    working_dir = argv[3]
    main(protein_dir, best_blast_dir, working_dir)

