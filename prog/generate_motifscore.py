#from motifscan outcome by homer to motif scores, note that multiple tfs may share the same motif, so the size motif score matrix is nb_peaks x nb_tfs, the highest motif score will assign to all correpsonding tfs
import numpy as np 
import pandas as pd 
import sys,os
import scipy.io as scio
import time

motif_file = '../data/motif_db/all_motif_rmdup.motif'
motif_tf_file = '../data/motif_db/MotifMatch_human_rmdup.mat'
peak_file = '../data/encode/preprocessed/union.peaks.pad1k.bed'
motifscan_file = '../data/encode/preprocessed/union.pad1k.motif.txt'


motif_tf_match = scio.loadmat(motif_tf_file)
motif2tf = {item[0][0]:[] for item in motif_tf_match['motifName']}
for each in motif_tf_match['Match2']:
    motif = each[0][0]
    tf = each[1][0]
    motif2tf[motif]+=[tf]

motif_set = []
for line in open(motif_file).readlines():
    if line[0]=='>':
        motif_set.append(line.split('\t')[1])
final_motif_set = [item for item in motif2tf.keys() if item in motif_set]
ind_motif={}
for each in motif_set:
    if each in final_motif_set:
        ind_motif[each] = 1
    else:
        ind_motif[each] = 0
len(motif2tf.keys())
len(motif_set)
len(final_motif_set)

final_tf_set = []
for each in final_motif_set:
    final_tf_set += motif2tf[each]
final_tf_set = list(set(final_tf_set))
final_tf_set.sort()
len(final_tf_set)

peakid2position = {item.split('\t')[3]:item.split('\t')[0]+':'+item.split('\t')[1]+'-'+item.split('\t')[2] for item in open(peak_file).readlines()}
nb_peaks = len(peakid2position.keys())
peaks_info = [item.split('\t')[0]+':'+item.split('\t')[1]+'-'+item.split('\t')[2] for item in open(peak_file).readlines()]

tf_score_mat = np.zeros((nb_peaks,len(final_tf_set)))
#df = pd.DataFrame(tf_score_mat,columns=final_tf_set,index=peaks_info)
col2idx = {item:i for i,item in enumerate(final_tf_set)}
f = open(motifscan_file,'r')
line = f.readline()
line = f.readline()
count=0
t=time.time()
while line!='':
    if count%500000==0:
        print(count)
    contents = line.split('\t')
    seq_id = contents[0]
    motif = contents[3]
    score = float(contents[-1].rstrip())
    if ind_motif[motif]==1:
       # print 'contents',time.time()-t
       # t=time.time()
        #print seq_id,motif,score
        tfs = motif2tf[motif]
        for each in tfs:
            col_idx = col2idx[each]
            #df[each].loc[peakid2position[seq_id]] = max(score,df[each].loc[peakid2position[seq_id]])
            tf_score_mat[int(seq_id)][col_idx] = max(score,tf_score_mat[int(seq_id)][col_idx])
        #print 'load score',time.time()-t
        #t=time.time()
    line = f.readline()
    #print 'read the next line',time.time()-t
    #t=time.time()
    count+=1
f.close()  
df = pd.DataFrame(tf_score_mat,columns=final_tf_set,index=peaks_info)
df.to_csv('../data/encode/preprocessed/tf_motif_score.csv',sep='\t')
