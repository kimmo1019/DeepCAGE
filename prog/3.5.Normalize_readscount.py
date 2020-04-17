#generate a table to store the normalized readscount, refer to BIRD paper for normalization
import sys
import numpy as np
import pandas as pd

usage='''
Usage: python 3.5.Normalize_readscount.py  [selected cell types] [output]
-- a program for generating readscount matrix (L x C)
OPTIONS:
    [selected cell types] -- selected cell types defined in 3.0
    [output]  --  readscount matrix as a plain text
'''
if len(sys.argv) != 3:
    print usage
    exit(0)
    
celltypes = [line.rstrip() for line in open(sys.argv[1]).readlines()]

peaks = {}
total_reads={}
for c in celltypes:
    total_reads[c] = int(open('../data/processed_RNA_DNase/%s.total.readscount'%c).readline())
min_reads = min(total_reads.values())   
#f_out = open('../data/processed_RNA_DNase/union.peaks.readscount','w')
label_file='../data/processed_RNA_DNase/union.peaks.labels'
label = pd.read_csv(label_file,sep='\t',header=0,index_col=[0])
rownames = label.index
reads_mat = np.zeros(label.shape)
cellidx=0
for c in celltypes:
    #print c
    regions_readscount = open('../data/processed_RNA_DNase/%s.union.peaks.readscount'%c).readlines()
    if len(regions_readscount)!=label.shape[0]:
        print 'reads count error!'
        sys.exit(1)
    reads_vec = np.array([int(item.strip().split('\t')[-1]) for item in regions_readscount])
    reads_mat[:,cellidx] = reads_vec*min_reads*1.0/total_reads[c]
    cellidx+=1

#formulate pandas dataframe
df = pd.DataFrame(data=reads_mat,columns=celltypes,index=rownames)
df.to_csv(sys.argv[2],sep='\t')
