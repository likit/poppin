import sys
from Bio import SeqIO


def find_max(infile):
    max_length = None
    for rec in SeqIO.parse(infile, 'fasta'):
        print rec.id, len(rec.seq)
        if not max_length:
            max_length = rec
        else:
            if len(max_length.seq) < len(rec.seq):
                max_length = rec

    # print "The longest sequence is", max_length.id, len(max_length.seq)
    SeqIO.write(max_length, infile + '.longest', 'fasta')


if __name__=='__main__':
    find_max(sys.argv[1])
