import sys
import glob
import os

for f in glob.glob(os.path.join(sys.argv[1], '*.fa.cap.contigs')):
    geneid = f.split('.')[2].split('_')[-1]
    tmp = '%s.tmp' % f
    tmpfile = open(tmp, 'w')
    for line in open(f):
        if line.startswith('>'):
            line = line.strip('\n').replace('>', '>' + geneid + '_')
            print >> tmpfile, line
        else:
            print >> tmpfile, line.strip('\n')
    tmpfile.close()
    os.rename(tmp, f)
