import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

'''Concatenate sequences from all input files to create
a supercontig/gene.

Usage: concat-seq.py <seqid> <input files>

'''
seqid = sys.argv[1]
superseq = SeqRecord(seq='')  # empty sequence
numseq = 0
for infile in sys.argv[2:]:
    print >> sys.stderr, 'working on %s' % infile
    for rec in SeqIO.parse(infile, 'fasta'):
        superseq += rec
        numseq += 1

superseq.id = seqid
superseq.description = ''
SeqIO.write(superseq, sys.stdout, 'fasta')

print >> sys.stderr, 'Total sequences: %d' % numseq
print >> sys.stderr, 'Total length: %d bp' % len(superseq)
