import os
import sys
import getopt
from Bio import SeqIO

def write_batch_file_forACDCNN(SEQFILE,PROFILE,PDBFILE='',CHAIN='',OUTFILE="batch_input_forACDCNN.txt"):
    out = open(OUTFILE,'w')
    seq = SeqIO.read(SEQFILE, "fasta").seq
    aminoacids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    for res in range(len(seq)):
        for aa in aminoacids:
            if seq[res]!=aa:
                out.write(f'{seq[res]}{res+1}{aa}')
                if(PROFILE):
                    out.write(f'\t{PROFILE}')
                if(PDBFILE):
                    out.write(f'\t{PDBFILE}')
                    if(CHAIN):
                        out.write(f'\t{CHAIN}')
                out.write("\n")
    out.close()

if __name__ == "__main__":
    
    opt_list, args = getopt.getopt(sys.argv[1:], 's:n:' ,["seqfile=", "profile=", "pdbfile=", "chain=", "outfile="])
    opt_list = dict(opt_list)
    SEQFILE = os.path.abspath(opt_list['--seqfile'])
    PROFILE = os.path.abspath(opt_list['--profile'])
    try:
        PDBFILE = os.path.abspath(opt_list['--pdbfile'])
    except KeyError:
        PDBFILE = ""
    try:
        CHAIN = opt_list['--chain']
    except KeyError:
        CHAIN = ""
    try:
        OUTFILE = os.path.abspath(opt_list['--outfile'])
    except KeyError:
        OUTFILE = "batch_input_forACDCNN.txt"
        
    write_batch_file_forACDCNN(SEQFILE,PROFILE,PDBFILE,CHAIN,OUTFILE)
    

    