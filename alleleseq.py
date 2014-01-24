'''Find allele sequences from assembly to be used as

an input for EcMLST website. Output sequences are of
the same length as provided reference sequences.

'''

import sys
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict

refseq = sys.argv[1]
blastinput = sys.argv[2]
refseq_len = {}

print >> sys.stderr, 'Reading reference sequences'.center(50, '=')
for rec in SeqIO.parse(refseq, 'fasta'):
    print >> sys.stderr, rec.id + ', %d bp' % len(rec)
    refseq_len[rec.id] = len(rec)

print >> sys.stderr, 'Getting sequences'.center(50, '=')
handle = open(blastinput)
for blastrec in NCBIXML.parse(handle):
    for alignment in blastrec.alignments:
        subject = alignment.hit_id.split('|')[-1]
        hps  = alignment.hsps[0]
        n = hps.sbjct_start
        i = 0
        seq = defaultdict(str)
        while n <= hps.sbjct_end:
            seq[n] = hps.query[i]
            i += 1
            n += 1
        alleleseq = ''
        for i in range(refseq_len[subject]):
            base = seq[i + 1]
            alleleseq += '-' if base == '' else base

        seqrec = SeqRecord(Seq(alleleseq))
        seqrec.id = subject
        seqrec.description = ''
        SeqIO.write(seqrec, sys.stdout, 'fasta')
        break  # ignore other alignments
