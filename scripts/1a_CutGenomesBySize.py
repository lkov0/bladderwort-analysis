__author__ = 'jgwall'

# Script to cut genome into chunks of user specified size

import argparse
import gzip
from Bio import SeqIO
import math


debug = False


def main():
    args = parse_args()
    IN = get_filehandle(args.infile, "rt")
    record_iter = SeqIO.parse(IN,"fasta")
    print("Splitting sequences in",args.infile,"into pieces of size",args.winsize)

    # Go over the records in each file, split, and print out
    i=0
    for record in record_iter:
        seq = record.seq
        numwins = math.ceil(len(seq) / args.winsize)
        print("\tProcessing",record.id,"of length", len(seq), "into",numwins, "windows")

        #Split into windows and write to separate files
        start, stop = 0, args.winsize
        for offset in range(numwins):
            # Get the current start and stop coordinates of the subsequence
            start = offset * args.winsize
            stop = start + args.winsize
            if stop > len(seq): stop = len(seq) # Make sure don't go over

            if stop - start < args.min:
                print("\tSkipped due to being less than",args.min,"bp")
                i+=1
                continue

            # Make output data
            subseq = record[start:stop]
            # newdescription = record.description.replace(" ","_") + "_window" + str(i) + "_" \
            #                  + str(start+1) + "_" + str(stop+1) # Add 1 because nucleotides usually start at 1, not 0
            outfile = args.outprefix + ".win" + str(i) + ".fa"
            if args.rename:
                subseq.id=args.rename + "." + "window" + str(i) + "_" + str(start+1) + "_" + str(stop+1)
                # record.description=newdescription


            OUT = get_filehandle(outfile, "wt")

            # Write out
            SeqIO.write(subseq, OUT, "fasta")

            # Cleanup
            OUT.close()
            i+=1
            if debug and i>=10: break
        if debug and i>=10: break

    IN.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outprefix", help="Prefix for output files")
    parser.add_argument("-r", "--rename", help="Change fasta names to this followed by the window")
    parser.add_argument("-w", "--winsize", type=int, default=1000000, help="Window size to cut each sequence into")
    parser.add_argument("-m", "--min", type=int, default=0, help="Minimum window size (windows less than this are skipped)")
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