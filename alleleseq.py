'''Find allele sequences from assembly to be used as

an input for EcMLST website. Output sequences are of
the same length as provided reference sequences.

'''

import sys
from Bio.Blast import NCBIXML
from Bio import SeqIO
from collections import defaultdict

geneSet = sys.argv[1]
blastOut = sys.argv[2]
refSeqLen = {}

print >> sys.stderr, 'Reading reference sequences'.center(50, '=')
for rec in SeqIO.parse(geneSet, 'fasta'):
    print >> sys.stderr, rec.id + ', %d bp' % len(rec)
    refSeqLen[rec.id] = len(rec)

print >> sys.stderr, 'Getting sequences'.center(50, '=')
handle = open(blastOut)
for blastRec in NCBIXML.parse(handle):
    for alignment in blastRec.alignments:
        subject = alignment.hit_id.split('|')[-1]
        print '>' + subject, blastRec.query
        hps  = alignment.hsps[0]
        n = hps.sbjct_start
        i = 0
        seq = defaultdict(str)
        while n <= hps.sbjct_end:
            seq[n] = hps.query[i]
            i += 1
            n += 1
        alleleSeq = ''
        for i in range(refSeqLen[subject]):
            base = seq[i + 1]
            alleleSeq += '-' if base == '' else base
        print alleleSeq, len(alleleSeq)
        break  # ignore other alignments
