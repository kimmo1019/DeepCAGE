# tf binding properties (binding number and coverage rate) in open chromain regions and random regions
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import hickle as hkl
import os,sys


cellline = 'K562'
exp = 'ENCSR000AEL'
motif_db = 'HOMER'

gene2expr = {}
for line in open('../data/%s/rseq/%s/0000'%(cellline,exp)).readlines():
    gene2expr[line.strip().split('\t')[1]] = float(line.strip().split('\t')[2])
tf_expressed = []

motifs2id = {name:idx for idx, name in enumerate(open('../data/motif_db/%s/motifs_list'%motif_db).readlines())}
for each in motif2id:
    if each.split('(')[0].upper() in gene2expr:
        tf_expressed.append(each.split('(')[0].upper())
#list to store TFs with different expression level
tf_expr_zero = []
tf_expr_low=[]
tf_expr_high=[]
ratio = 0.5
tf_expr_sort=sorted([[item,gene2expr[item]] for item in tf_expressed],key=lambda a:a[1])
tf_expr_zero = [item[0] for item in tf_expr_sort if item[1]==0]
num_zero = len(tf_expr_zero)
num_low = int((len(tf_expressed)-len(tf_expr_zero))*ratio)
tf_expr_low = [item[0] for item in tf_expr_sort][num_zero:(num_zero+num_low)]
tf_expr_high = [item[0] for item in tf_expr_sort][(num_zero+num_low):]
#print len(tf_expr_zero),len(tf_expr_low),len(tf_expr_high),len(tf_expressed)
cmd = 'wc -l  ../data/%s/dseq/positive_indexd.bed'%cellline
P_num = int(os.popen(cmd).readlines()[0].split()[0])
cmd = 'wc -l  ../data/%s/dseq/negative_indexd.bed'%cellline
N_num = int(os.popen(cmd).readlines()[0].split()[0])

#vector to store the number of TFBS on each positive/negative region
P_binding_vec_zero = np.zeros(P_num)
P_binding_vec_low = np.zeros(P_num)
P_binding_vec_high = np.zeros(P_num)
N_binding_vec_zero = np.zeros(N_num)
N_binding_vec_low = np.zeros(N_num)
N_binding_vec_high = np.zeros(N_num)
for each in open('../data/%s/dseq/%s_overlap_filter_sorted_pos.bed'%(cellline,motif_db)).readlines():
    index = int(each.strip().split('\t')[3])
    tf_name = each.strip().split('\t')[8].split('(')[0].upper()
    if tf_name in tf_expr_zero:
        P_binding_vec_zero[index-1]+=1
    if tf_name in tf_expr_low:
        P_binding_vec_low[index-1]+=1
    if tf_name in tf_expr_high:
        P_binding_vec_high[index-1]+=1
for each in open('../data/%s/dseq/%s_overlap_filter_sorted_neg.bed'%(cellline,motif_db)).readlines():
    index = int(each.strip().split('\t')[3])
    tf_name = each.strip().split('\t')[8].split('(')[0].upper()
    if tf_name in tf_expr_zero:
        N_binding_vec_zero[index-1]+=1
    if tf_name in tf_expr_low:
        N_binding_vec_low[index-1]+=1
    if tf_name in tf_expr_high:
        N_binding_vec_high[index-1]+=1
#draw boxplot 
dic={}
dic['TF binding number'] = np.concatenate([P_binding_vec_zero,P_binding_vec_low,P_binding_vec_high,N_binding_vec_zero,N_binding_vec_low,N_binding_vec_high],axis=0)
dic['region type'] = ['DNase peak']*3*P_num+['Random region']*3*N_num
dic['TF expression level'] = ['Zero']*P_num+['Low']*P_num+['High']*P_num + ['Zero']*P_num+['Low']*P_num+['High']*P_num
df_tfbs = pd.DataFrame(dic)
#store data
hkl.dump(df_tfbs,'../data/%s/tfbs_num_dis.hkl'%cellline,mode='w',compression='gzip')
ax = sns.boxplot(x='TF expression level', y="TF binding number", hue='region type',
                    data=df_tfbs, palette="Set3")
#plt.ylim([0,90])
plt.save('../data/%s/tfbs_num_dis.pdf')
#plt.show()
     
    
    
#vector to store the coverage rate of TFBS on each positive/negative region
P_cover_vec_zero = np.zeros(P_num)
P_cover_vec_low = np.zeros(P_num)
P_cover_vec_high = np.zeros(P_num)
N_cover_vec_zero = np.zeros(N_num)
N_cover_vec_low = np.zeros(N_num)
N_cover_vec_high = np.zeros(N_num)

region_previous_idx = 0
isFirstSample = 1
for each in open('../data/%s/dseq/%s_overlap_filter_sorted_pos.bed'%(cellline,motif_db)).readlines():
    index = int(each.strip().split('\t')[3])
    tf_name = each.strip().split('\t')[8].split('(')[0].upper()
    region_start = int(each.strip().split('\t')[1])
    region_end = int(each.strip().split('\t')[2])
    if index != region_previous_idx:
        if isFirstSample != 1:
            P_cover_vec_zero[region_previous_idx-1] = np.mean(cover_vec_current_zero)
            P_cover_vec_low[region_previous_idx-1] = np.mean(cover_vec_current_low)
            P_cover_vec_high[region_previous_idx-1] = np.mean(cover_vec_current_high)
        cover_vec_current_zero = np.zeros(region_end-region_start+1)
        cover_vec_current_low = np.zeros(region_end-region_start+1)
        cover_vec_current_high = np.zeros(region_end-region_start+1)
    tf_start = int(each.strip().split('\t')[6])
    tf_end = int(each.strip().split('\t')[7])
    overlap_start = max(tf_start,region_start)
    overlap_end = min(tf_end,region_end)
    if tf_name in tf_expr_zero:
        #note the cooridate of a vector should start at 0
        cover_vec_current_zero[(overlap_start-region_start):(overlap_end+1-region_start)] = 1
    if tf_name in tf_expr_low:
        cover_vec_current_low[(overlap_start-region_start):(overlap_end+1-region_start)] = 1
    if tf_name in tf_expr_high:
        cover_vec_current_high[(overlap_start-region_start):(overlap_end+1-region_start)] = 1
    isFirstSample=0
    region_previous_idx = index
# for the last region
P_cover_vec_zero[region_previous_idx-1] = np.mean(cover_vec_current_zero)
P_cover_vec_low[region_previous_idx-1] = np.mean(cover_vec_current_low)
P_cover_vec_high[region_previous_idx-1] = np.mean(cover_vec_current_high)

region_previous_idx = 0
isFirstSample = 1
for each in open('../data/%s/dseq/%s_overlap_filter_sorted_neg.bed'%(cellline,motif_db)).readlines():
    index = int(each.strip().split('\t')[3])
    tf_name = each.strip().split('\t')[8].split('(')[0].upper()
    region_start = int(each.strip().split('\t')[1])
    region_end = int(each.strip().split('\t')[2])
    if index != region_previous_idx:
        if isFirstSample != 1:
            N_cover_vec_zero[region_previous_idx-1] = np.mean(cover_vec_current_zero)
            N_cover_vec_low[region_previous_idx-1] = np.mean(cover_vec_current_low)
            N_cover_vec_high[region_previous_idx-1] = np.mean(cover_vec_current_high)
        cover_vec_current_zero = np.zeros(region_end-region_start+1)
        cover_vec_current_low = np.zeros(region_end-region_start+1)
        cover_vec_current_high = np.zeros(region_end-region_start+1)
    tf_start = int(each.strip().split('\t')[6])
    tf_end = int(each.strip().split('\t')[7])
    overlap_start = max(tf_start,region_start)
    overlap_end = min(tf_end,region_end)
    if tf_name in tf_expr_zero:
        cover_vec_current_zero[(overlap_start-region_start):(overlap_end+1-region_start)] = 1
    if tf_name in tf_expr_low:
        cover_vec_current_low[(overlap_start-region_start):(overlap_end+1-region_start)] = 1
    if tf_name in tf_expr_high:
        cover_vec_current_high[(overlap_start-region_start):(overlap_end+1-region_start)] = 1
    isFirstSample=0
    region_previous_idx = index
# for the last region
N_cover_vec_zero[region_previous_idx-1] = np.mean(cover_vec_current_zero)
N_cover_vec_low[region_previous_idx-1] = np.mean(cover_vec_current_low)
N_cover_vec_high[region_previous_idx-1] = np.mean(cover_vec_current_high)

dic={}
dic['TF coverage rate'] = np.concatenate([P_cover_vec_zero,P_cover_vec_low,P_cover_vec_high,N_cover_vec_zero,N_cover_vec_low,N_cover_vec_high],axis=0)
dic['region type'] = ['DNase peak']*3*P_num+['Random region']*3*N_num
dic['TF expression level'] = ['Zero']*P_num+['Low']*P_num+['High']*P_num + ['Zero']*P_num+['Low']*P_num+['High']*P_num

df_tfbs = pd.DataFrame(dic)
#store data
hkl.dump(df_tfbs,'../data/%s/tfbs_cover_dis.hkl'%cellline,mode='w',compression='gzip')
ax = sns.boxplot(x='TF expression level', y="TF coverage rate", hue='region type',
                    data=df_tfbs, palette="Set3")
#plt.ylim([0,1.4])
plt.save('../data/%s/tfbs_cover_dis.pdf')
#plt.show()
