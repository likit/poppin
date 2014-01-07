import sys
import subprocess
import glob
import os
from collections import defaultdict

def pslx_to_fasta(blat_output, outdir):
    '''Reads pslx file and convert it to fasta format.

    Every alignment is saved as an individual fasta record.
    Sequence ID is query ID + an alignment number.
    '''

    currdir = os.path.abspath(os.curdir)
    os.chdir(outdir)

    queries = defaultdict(int)
    curr_query = None
    for line in open(blat_output):
        cols = line.split('\t')
        new_query = cols[9]
        sequences = cols[22].split(',')[:-1]

        if curr_query == None:
            curr_query = new_query
            op = open('%s_%s.fa' % (blat_output, curr_query), 'w')

        if curr_query != new_query:
            op.close()
            curr_query = new_query
            op = open('%s_%s.fa' % (blat_output, curr_query), 'w')

        for seq in sequences:
            seqid = '>%s-%d' % (curr_query, queries[curr_query])
            queries[curr_query] += 1
            print >> op, '%s\n%s' % (seqid, seq)
    os.chdir(currdir)

def run_blat(contig, query, blat_output, outdir):
    '''run BLAT in shell'''
    blat_output = os.path.join(outdir, blat_output)

    return subprocess.call('blat -noHead -out=pslx %s %s %s' % (contig, query, blat_output), shell=True)


def run_cap3(outdir):
    '''run CAP3 in shell on all files'''
    currdir = os.path.abspath(os.curdir)
    os.chdir(outdir)
    num_files = 0
    for infile in glob.glob('*pslx*.fa'):
        num_files += 1
        ret = subprocess.call('cap3 %s' % infile,
                                shell=True,
                                stdout=open('%s.cap.align' % infile, 'w'),
                            )
        if ret != 0:
            return ret
    os.chdir(currdir)
    if num_files == 0:  # files not found
        return 1
    else:
        return 0


def main():
    if len(sys.argv) < 3:
        print >> sys.stderr, 'Not enough arguments.'
        print >> sys.stderr, 'Usage: poppin.py <contig.fa> <query.fa> <outdir>'
        sys.exit(1)

    contig = sys.argv[1]
    query = sys.argv[2]
    blat_output = os.path.split(contig)[-1] + '.pslx'
    outdir = os.path.abspath(sys.argv[3])

    if os.path.exists(outdir):
        print >> sys.stderr, '%s exists.' % outdir
        sys.exit(1)
    else:
        os.mkdir(outdir)

    print >> sys.stderr, 'Step 1: align query sequences to contig(s)'
    print >> sys.stderr, '=' * 50
    ret = run_blat(contig, query, blat_output, outdir)  # if return code not equal to 0
    if ret != 0:
        print 'BLAT failed. Please check the input files.'
        sys.exit(ret)

    print >> sys.stderr, 'Step 2: convert alignments to FASTA'
    print >> sys.stderr, '=' * 50
    pslx_to_fasta(blat_output, outdir)

    print >> sys.stderr, 'Step 3: build consensus sequences using CAP3'
    print >> sys.stderr, '=' * 50
    ret = run_cap3(outdir)
    if ret != 0:
        print 'CAP3 failed. Please check the input files.'
        sys.exit(ret)

    print >> sys.stderr, 'Done.'

if __name__=='__main__':
    main()
