import sys
from Bio import SeqIO

'''Finds the longest sequence in a FASTA file and writes
it to stdout.

'''

def find_max(infile):
    max_length = None
    for rec in SeqIO.parse(infile, 'fasta'):
        # print rec.id, len(rec.seq)
        if not max_length:
            max_length = rec
        else:
            if len(max_length.seq) < len(rec.seq):
                max_length = rec

    # print "The longest sequence is", max_length.id, len(max_length.seq)
    SeqIO.write(max_length, sys.stdout, 'fasta')


if __name__=='__main__':
    for infile in sys.argv[1:]:
        find_max(infile)
