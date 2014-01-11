import glob
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

'''Concatenate sequences from all input files to create
a supercontig/gene.

Usage: concat-seq.py <input files> <seqid>

'''

try:
    seqid = sys.argv[2]
except IndexError:
    seqid = 'supercontig'

superseq = SeqRecord('')  # empty sequence
superseq.id = seqid
numseq = 0
for infile in glob.glob(sys.argv[1]):
    for rec in SeqIO.parse(infile, 'fasta'):
        superseq += rec
        numseq += 1

SeqIO.write(superseq, sys.stdout, 'fasta')

print >> sys.stderr, 'Total sequences: %d' % numseq
print >> sys.stderr, 'Total length: %d bp' % len(superseq)
