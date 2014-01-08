import sys
import glob
import os

for f in glob.glob(os.path.join(sys.argv[1], '*.fa.cap.contigs')):
    geneid = f.split('.')[2].split('_')[-1]
    tmp = '%s.tmp' % f
    tmp = open(tmp, 'w')
    for line in open(f):
        if line.startswith('>'):
            line.strip('\n').replace('Contig', geneid + '-')
            print >> tmp, line
        else:
            print >> tmp, line.strip('\n')
    tmp.close()
