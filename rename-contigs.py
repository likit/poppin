import sys
import glob
import os

for f in glob.glob(os.path.join(sys.argv[1], '*.fa.cap.contigs')):
    tmp = '%s.tmp' % f
    tmpfile = open(tmp, 'w')
    seqid = os.path.split(f)[-1]
    for line in open(f):
        if line.startswith('>'):
            line = line.strip('\n').replace('>', '>' + seqid + '_')
            print >> tmpfile, line
        else:
            print >> tmpfile, line.strip('\n')
    tmpfile.close()
    os.rename(tmp, f)
