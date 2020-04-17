import numpy as np 
import pandas as pd
import hickle as hkl
import os ,sys

nb_peaks = 15181508
nb_samples = 129
dic_total = {}
df = pd.DataFrame(np.zeros((nb_peaks,nb_samples)))
idx=0
for each in open('../data/encode/bams_info_final.txt').readlines():

    rseq_exp = each.split('\t')[0]
    dseq_exp = each.split('\t')[1]
    print idx, dseq_exp
    file_peaks = '../data/encode/preprocessed/%s.peaks.readscount'%dseq_exp
    file_total = '../data/encode/preprocessed/%s.total.readscount'%dseq_exp
    cmd = ''' cat %s '''%file_total
    nb_readcounts = int(os.popen(cmd).read().strip())
    dic_total[dseq_exp] = nb_readcounts
    #cmd = ''' wc -l %s'''%file_peaks
    #nb_peaks = os.popen(cmd).read().strip().strip()
    df.index = [item.split('\t')[0]+':'+item.split('\t')[1]+'-'+item.split('\t')[2] for item in open(file_peaks).readlines()]
    readscount_vec = np.array([float(item.split('\t')[-1].strip()) for item in open(file_peaks).readlines()])
    df[idx] = readscount_vec
    df.rename(columns={idx:dseq_exp},inplace=True)
    idx+=1
df.to_csv('../data/encode/readscount_raw.csv',sep='\t')
hkl.dump(dic_total,'../data/encode/readscount_totol_dic.hkl')
