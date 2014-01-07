import sys
import subprocess
import glob
from collections import defaultdict

def pslx_to_fasta(output):
    queries = defaultdict(int)
    '''Reads pslx file and convert it to fasta format.

    Every alignment is saved as an individual fasta records.
    Sequence ID is query ID + the number of alignments.
    '''
    curr_query = None
    for line in open(output):
        cols = line.split('\t')
        new_query = cols[9]
        sequences = cols[22].split(',')[:-1]

        if curr_query == None:
            curr_query = new_query
            op = open('%s_%s.fa' % (output, curr_query), 'w')

        if curr_query != new_query:
            op.close()
            curr_query = new_query
            op = open('%s_%s.fa' % (output, curr_query), 'w')

        for seq in sequences:
            seqid = '>%s-%d' % (curr_query, queries[curr_query])
            queries[curr_query] += 1
            print >> op, '%s\n%s' % (seqid, seq)

def run_blat(contig, query, output):
    '''run BLAT in shell'''

    return subprocess.call('blat -noHead -out=pslx %s %s %s' % (contig, query, output), shell=True)


def run_cap3():
    for infile in glob.glob('*pslx*.fa'):
        ret = subprocess.call('cap3 %s' % infile,
                                shell=True,
                                stdout=open('/dev/null', 'w'),
                            )
        if ret != 0:
            return ret
    return 0


def main():
    if len(sys.argv) < 2:
        print >> sys.stderr, 'Not enough arguments.'
        print >> sys.stderr, 'Usage: poppin.py <contig.fa> <query.fa>'
        sys.exit(1)

    contig = sys.argv[1]
    query = sys.argv[2]
    output = contig + '.pslx'

    print >> sys.stderr, 'Step 1: align query sequences to contig(s)'
    print >> sys.stderr, '=' * 50
    ret = run_blat(contig, query, output)  # if return code not equal to 0
    if ret != 0:
        print 'BLAT failed. Please check the input files.'
        sys.exit(ret)

    print >> sys.stderr, 'Step 2: convert alignments to FASTA'
    print >> sys.stderr, '=' * 50
    pslx_to_fasta(output)

    print >> sys.stderr, 'Step 3: finding consensus sequences using CAP3'
    print >> sys.stderr, '=' * 50
    ret = run_cap3()
    if ret != 0:
        print 'CAP3 failed. Please check the input files.'
        sys.exit(ret)

if __name__=='__main__':
    main()
