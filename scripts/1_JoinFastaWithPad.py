__author__ = 'jgwall'

import argparse
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import math


debug = False


def main():
    args = parse_args()
    print("Joinging sequences in",args.infile,"and padding with",args.n_padding,"Ns between each")
    
    # Open files
    IN = get_filehandle(args.infile, "rt")
    record_iter = SeqIO.parse(IN,"fasta")
    KEY = get_filehandle(args.keyfile, "wt")
    KEY.write("start\tstop\tscaffold\n")
    
    # Go through and join
    i=0
    joined=Seq("", generic_dna)
    pad=Seq("N" * args.n_padding, generic_dna)
    for record in record_iter:
        i+=1
        if debug and i>=10: break
        if i % 1000 ==0 : print("Joined",i,"sequences")
        
        # Start and stop for the key
        
        
        if len(joined)==0:
            start = 0
            joined += record.seq    # No need for padding the first time
        else:
            start=len(joined) + len(pad)
            joined = joined + pad + record.seq
        stop=start + len(record) - 1
            
        # Write key record
        KEY.write(str(start) + "\t" + str(stop) + "\t" + record.description + "\n")
    print("\tJoined",i,"total sequences")
    IN.close()    
    
    # Write out
    OUT = get_filehandle(args.outfile, "wt")
    #OUT.write(">" + args.newname + "\n")
    newseq = SeqRecord(joined)
    newseq.id = args.newname
    SeqIO.write(newseq, OUT, "fasta")
    OUT.close()
    

    

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile", help="Output file of joined FASTA sequences")
    parser.add_argument("-k", "--keyfile", help="Output keyfile showing location of each original scaffold")
    parser.add_argument("--newname", default="joined_sequence", help="Name for the joined sequence")
    parser.add_argument("-n", "--n-padding", type=int, default=1000, help="Number of N's to stick in between each sequence")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def get_filehandle(file, mode):
    if file.endswith(".gz"):
        return gzip.open(file, mode)
    return open(file, mode)




if __name__ == '__main__': main()