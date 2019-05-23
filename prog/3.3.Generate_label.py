import sys

usage='''
Usage: python 3.3.Generate_label.py  [peaks] [selected cell types] [output]
-- a program for generating label matrix (L x C)
OPTIONS:
    [peaks] -- union peak bed file in 3.0
    [selected cell types] -- selected cell types defined in 3.0
    [output]  --  label matrix as a plain text
'''

if len(sys.argv) != 4:
    print usage
    exit(0)

celltypes = [line.rstrip() for line in open(sys.argv[2]).readlines()]

peaks = {}
for c in celltypes:
    peaks[c] = {line.split()[0]+':'+line.split()[1]+'-'+line.rstrip().split()[2]:1 for line in open('../data/processed_RNA_DNase/%s.peak.bins.bed' % (c)).readlines()}

fd = open(sys.argv[1])
gd = open(sys.argv[3],'w')
gd.write('#region\t' + '\t'.join(celltypes) + '\n')
line = fd.readline()
while not line == '':
    tag = line.split()[0]+':'+line.split()[1]+'-'+line.rstrip().split()[2]
    _labels = []
    for c in celltypes:
        if peaks[c].has_key(tag): _labels.append('1')
        else: _labels.append('0')
    gd.write(tag + '\t' + '\t'.join(_labels) + '\n')
    line = fd.readline()
gd.close()
fd.close()
