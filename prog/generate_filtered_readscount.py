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
#final_peaks=hkl.load('../data/encode/peak_final_thred8-10000-3.5.hkl')
final_peaks=hkl.load('../data/encode/peak_final_thred10-10000-8.hkl')
print('finally, %d peaks were found'%(len(final_peaks)))
df_filtered = df.loc[final_peaks]
df_filtered.to_csv('../data/encode/readscount_normalized_filtered.csv',sep='\t')
