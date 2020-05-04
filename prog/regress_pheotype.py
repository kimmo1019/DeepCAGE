import numpy as np 
import pandas as pd 
import os, sys
import gzip
import hickle as hkl
from pyfasta import Fasta
from sklearn.linear_model import LinearRegression, Ridge, Lasso
import scipy.stats as ss
from sklearn.metrics import r2_score
from sklearn.model_selection import cross_validate

def get_causal_score(rs):
    genotype_head_file = '%s/GTEx/GTEx_WGS_genotype_first500.txt' %DPATH
    muscle_info = '%s/GTEx/GTEx_rnaseq_mucsle.txt'%DPATH
    first500 = open(genotype_head_file).readlines()
    line = first500[1].strip().split('\t')
    donor_muscle_ids = [item.split('\t')[0] for item in open(muscle_info).readlines()]
    donor_ids = ['-'.join(item.rstrip().split('-')[:2]) for item in donor_muscle_ids]
    muscle_tissue_ids=[donor_muscle_ids[donor_ids.index(item)] for item in line[9:] if item in donor_ids]

    Y_pred_ref,Y_pred_alt = hkl.load('%s/height_gwas/neighbor_var/%s.score.hkl'%(DPATH,rs))
    var_info = '%s/height_gwas/neighbor_var/%s.hkl'%(DPATH,rs)
    vars_list = hkl.load(var_info)
    rs2info = {item.split('\t')[0]:item.split('\t')[1:5] for item in open(height_GWAS).readlines()[1:]}
    pos_wgas = rs2info[rs][1]
    score_mat = np.zeros((len(muscle_tissue_ids),len(vars_list)))
    idx=0
    epsilon=1e-10
    columns = ['_'.join([item[0],str(item[1]),item[2],item[3]]) for item in vars_list]
    for var in vars_list:
        chrom, pos, ref, alt, related_donor_idx = var[0], var[1], var[2], var[3], var[4]
        if len(ref)==1 and len(alt)==1:
            for i in related_donor_idx:
                score_mat[i][vars_list.index(var)] = abs(np.log2(abs(Y_pred_ref[idx,0]*1.0/Y_pred_alt[idx,0])+epsilon))
                idx+=1
    pd_data = pd.DataFrame(score_mat,index=muscle_tissue_ids,columns=columns)
    pd_data.to_csv('%s/height_gwas/neighbor_var/%s.causal_score.csv'%(DPATH,rs),sep='\t')

def get_height_mat(rs):
    pos_wgas = rs2info[rs][1]
    score_file = '%s/height_gwas/neighbor_var/%s.causal_score_thred1.csv'%(DPATH,rs)
    causal_score = pd.read_csv(score_file ,sep='\t',header=0,index_col=[0])
    return causal_score.values, causal_score.columns #may duplicate




if __name__=="__main__":
    DPATH='/home/liuqiao/software/DeepCAGE/data'
    height_GWAS='%s/height_gwas/Meta-analysis_Wood_et_al+UKBiobank_2018_top_3290_from_COJO_analysis.txt'%DPATH
    phenotype_file ='%s/GTEx/phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt'%DPATH 
    GWAS_variants = '%s/height_gwas/neighbor_var'%DPATH
    rs2info = {item.split('\t')[0]:item.split('\t')[1:5] for item in open(height_GWAS).readlines()[1:]}
    dic_height = {item.split('\t')[1]:float(item.split('\t')[7]) for item in open(phenotype_file).readlines()[11:]}

    rs_list = [item.split('\t')[0] for item in open(height_GWAS).readlines()[1:]]
    
    genotype_head_file = '%s/GTEx/GTEx_WGS_genotype_first500.txt' %DPATH
    muscle_info = '%s/GTEx/GTEx_rnaseq_mucsle.txt'%DPATH
    first500 = open(genotype_head_file).readlines()
    line = first500[1].strip().split('\t')
    donor_muscle_ids = [item.split('\t')[0] for item in open(muscle_info).readlines()]
    donor_ids = ['-'.join(item.rstrip().split('-')[:2]) for item in donor_muscle_ids]
    muscle_tissue_ids=[donor_muscle_ids[donor_ids.index(item)] for item in line[9:] if item in donor_ids]
    height_vec = [dic_height['-'.join(item.rstrip().split('-')[:2])] for item in muscle_tissue_ids]
    #min-max norm
    Y_height = (height_vec - np.min(height_vec))/(np.max(height_vec)-np.min(height_vec))
    # results = map(get_height_mat,rs_list[135:150])
    # all_score_mat = np.hstack([item[0] for item in results])
    # col_all = []
    # for each in [item[1] for item in results]:
    #     col_all += list(each)
    # pd_data_all = pd.DataFrame(all_score_mat,index=muscle_tissue_ids,columns=col_all)
    # pd_data_all = pd_data_all.T.drop_duplicates().T
    pd_data_all = pd.read_csv('%s/height_gwas/vars_score.csv'%DPATH,sep='\t',header=0,index_col=[0])
    var2score={}
    for each in pd_data_all.columns:
        col_data = pd_data_all[each]
        if len(col_data.shape)>1:
            col_data = col_data.values[:,0]
        nb_donor_has_var = len([value for value in col_data if value>0])
        var2score[each] = np.sum(col_data)*1.0/nb_donor_has_var

    
    #pd_data_all.to_csv('%s/height_gwas/vars_score.csv'%DPATH,sep='\t')
    hkl.dump(var2score,'%s/height_gwas/var2score.hkl'%DPATH)
    score_info = [[var2score[item], item] for item in pd_data_all.columns]
    score_info_sort = sorted(score_info, key=lambda a: a[0]) 

    col_top = [item[1] for item in score_info_sort]
    X_score = pd_data_all[col_top].values
    print X_score.shape, Y_height.shape
    #reg_all = LinearRegression()
    #reg_all = Ridge(alpha=0.1)
    reg_all = Lasso(alpha=0.1)
    cv_results = cross_validate(reg_all, X_score, Y_height, scoring='r2', cv=5)
    print cv_results['test_score'],np.mean(cv_results['test_score'])
    reg_all = Ridge(alpha=0.1)
    cv_results = cross_validate(reg_all, X_score, Y_height, scoring='r2', cv=5)
    print cv_results['test_score'],np.mean(cv_results['test_score'])

    col_top = [item[1] for item in score_info_sort[-100:]]
    X_score = pd_data_all[col_top].values
    print X_score.shape, Y_height.shape
    #reg_all = LinearRegression()
    #reg_all = Ridge(alpha=0.1)
    reg_all = Lasso(alpha=0.1)
    cv_results = cross_validate(reg_all, X_score, Y_height, scoring='r2', cv=5)
    print cv_results['test_score'],np.mean(cv_results['test_score'])
    reg_all = Ridge(alpha=0.1)
    cv_results = cross_validate(reg_all, X_score, Y_height, scoring='r2', cv=5)
    print cv_results['test_score'],np.mean(cv_results['test_score'])

    col_top = [item[1] for item in score_info_sort[:100]]
    X_score = pd_data_all[col_top].values
    print X_score.shape, Y_height.shape
    #reg_all = LinearRegression()
    #reg_all = Ridge(alpha=0.1)
    reg_all = Lasso(alpha=0.1)
    cv_results = cross_validate(reg_all, X_score, Y_height, scoring='r2', cv=5)
    print cv_results['test_score'],np.mean(cv_results['test_score'])
    reg_all = Ridge(alpha=0.1)
    cv_results = cross_validate(reg_all, X_score, Y_height, scoring='r2', cv=5)
    print cv_results['test_score'],np.mean(cv_results['test_score'])



    

