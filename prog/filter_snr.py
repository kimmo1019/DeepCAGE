import numpy as np 
import hickle as hkl 
import os
import sys
import pandas as pd 
import math
from scipy import signal


df = pd.read_csv('../data/encode/readscount_raw.csv',sep='\t',header=0,index_col=[0])

dic_total = hkl.load('../data/encode/readscount_totol_dic.hkl')

readscount_max = np.max(dic_total.values())
readscount_min = np.min(dic_total.values())


for each in df.columns:
    df[each] = df[each].values*readscount_min*1.0/dic_total[each]

#df.to_csv('../data/encode/readscount_normalized.csv',sep='\t')
peak_max_count = df.max(axis=1)
is_peak1 = np.zeros(df.shape[0])# 0 for discard and 1 for contain, thred fitler
is_peak2 = np.zeros(df.shape[0])# 0 for discard and 1 for contain, snr filter
idx = 0
peak_list=[]
for each in df.index:
    if peak_max_count.loc[each] > 10 and peak_max_count.loc[each] < 10000:
        is_peak1[idx] = 1
        peak_list.append(each)
    idx+=1
print('%d peaks were found using lower and upper threds'%len(peak_list))

#df = df.loc[peak_list]
dic_chrom_len = {item.split('\t')[0]:(int(item.split('\t')[-1].strip())-int(item.split('\t')[-1].strip())%200) for item in open('/home/liuqiao/openness_pre/chromosome.txt').readlines()}
print('Data loading done!')

def FilterSNR(x,res=200,nb_bg = 500):#x is a column of df
    idx = 0
    eps = 1e-10
    nb = len(df.index)
    snr = np.zeros(nb)
    for each in df.index:
        chrom, start, end  = each.split(':')[0], int(each.split(':')[1].split('-')[0]),int(each.split(':')[1].split('-')[1].strip())
        chrom_length = dic_chrom_len[chrom]
        if start < res*nb_bg/2:
            nb_before = start/res
            average_DH = (np.sum(x[(idx-nb_before):idx])+np.sum(x[(idx+1):(idx+1+nb_bg/2)]))*1.0/(nb_before+nb_bg/2)
            snr[idx] = x[idx]*1.0/(average_DH+eps)
        elif end > dic_chrom_len[chrom]-res*nb_bg/2:
            nb_after = (dic_chrom_len[chrom]-end)/res
            average_DH = (np.sum(x[(idx-nb_bg/2):idx])+np.sum(x[(idx+1):(idx+1+nb_after)]))*1.0/(nb_bg/2+nb_after)
            snr[idx] = x[idx]*1.0/(average_DH+eps)
        else:
            average_DH = (np.sum(x[idx-nb_bg/2:idx])+np.sum(x[idx+1:idx+1+nb_bg/2]))*1.0/nb_bg
            snr[idx] = x[idx]*1.0/(average_DH+eps)
        idx+=1
    return snr

if False:
    eps = 1e-10
    window = np.ones((500,1))/500.0
    average_signal = signal.convolve2d(df.values, window, mode='same', boundary='fill')
    snr_mat = df.values/(average_signal+eps)
    np.save('../data/encode/snr.npy',snr_mat)
snr_mat = np.load('../data/encode/snr.npy')
print('snr loading finished')
#df_snr = df.apply(FilterSNR,axis=0)
snr_max = snr_mat.max(axis=1)
idx=0
for each in df.index:
    if snr_max[idx]>8:
        is_peak2[idx] = 1
    idx+=1
is_peak = is_peak1*is_peak2
final_peaks = [item for i,item in enumerate(df.index) if is_peak[i]==1]
print('finally, %d peaks were found'%(len(final_peaks)))
hkl.dump(final_peaks,'../data/encode/peak_final_thred10-10000-8.hkl')

