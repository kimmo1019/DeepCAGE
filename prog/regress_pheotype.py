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
from sklearn.linear_model import LassoCV
from sklearn.model_selection import KFold

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
    
    genotype_head_file = '%s/GTEx/GTEx_WGS_genotype_first500.txt'%DPATH
    muscle_info = '%s/GTEx/GTEx_rnaseq_mucsle.txt'%DPATH
    first500 = open(genotype_head_file).readlines()
    line = first500[1].strip().split('\t')
    donor_muscle_ids = [item.split('\t')[0] for item in open(muscle_info).readlines()]
    donor_ids = ['-'.join(item.rstrip().split('-')[:2]) for item in donor_muscle_ids]
    muscle_tissue_ids=[donor_muscle_ids[donor_ids.index(item)] for item in line[9:] if item in donor_ids]
    height_vec = [dic_height['-'.join(item.rstrip().split('-')[:2])] for item in muscle_tissue_ids]

    #min-max norm
    #print height_vec
    #Y_height = (height_vec - np.min(height_vec))/(np.max(height_vec)-np.min(height_vec))
    #print Y_height
    #Y_height = np.array(height_vec,dtype='float64')/np.max(height_vec)
    Y_height = np.array(height_vec,dtype='float64')
    #print Y_height




    rs_list.remove('rs13429408')
        
    results = map(get_height_mat,rs_list)
    print 'load done'
    # all_mats_list = [item[0] for item in results]
    # average_list = [np.mean(item) for item in all_mats_list]
    # merge_list = sorted(zip(average_list,np.arange(len(average_list))),key=lambda a:a[0],reverse=True)
    # top_100_idx_list = [item[1] for item in merge_list[:100]]
    # top_100_score_mat = []
    # for i in top_100_idx_list:
    #     top_100_score_mat.append(all_mats_list[i])
    # top_100_score_mat = np.hstack(top_100_score_mat)


    # all_score_mat = np.hstack(all_mats_list)
    # col_all = []
    # for each in [item[1] for item in results]:
    #     col_all += list(each)
    # pd_data_all = pd.DataFrame(all_score_mat,index=muscle_tissue_ids,columns=col_all)
    # pd_data_all = pd_data_all.T.drop_duplicates().T
    # pd_data_all = pd_data_all.to_csv('%s/height_gwas/vars_score.csv'%DPATH,sep='\t')

    pd_data_all = pd.read_csv('%s/height_gwas/vars_score.csv'%DPATH,sep='\t',header=0,index_col=[0])
    print 'pd_data_all shape:', pd_data_all.shape
    # var2score={}
    # for each in pd_data_all.columns:
    #     col_data = pd_data_all[each]
    #     if len(col_data.shape)>1:
    #         col_data = col_data.values[:,0]
    #     nb_donor_has_var = len([value for value in col_data if value>0])
    #     var2score[each] = np.sum(col_data)*1.0/nb_donor_has_var
    # hkl.dump(var2score,'%s/height_gwas/var2score.hkl'%DPATH)
    var2score = hkl.load('%s/height_gwas/var2score.hkl'%DPATH)


    # #all vars in all loci
    # score_info = [[var2score[item], item] for item in pd_data_all.columns]
    # score_info_sort = sorted(score_info, key=lambda a: a[0]) 
    
    # col_top = [item[1] for item in score_info_sort]
    # #col_top = [item[1] for item in score_info_sort if sum(pd_data_all[item[1]]!=0)>10 and sum(pd_data_all[item[1]]!=0)<481 ]
    # X_score = pd_data_all[col_top].values
    # np.save('%s/height_gwas/vars_score_all.npy'%DPATH,X_score)
    # print 'total shape',X_score.shape
    # kf = KFold(n_splits=10)
    # r2_train_list = []
    # r2_test_list = []
    # for train_index, test_index in kf.split(X_score):
    #     X_train, X_test = X_score[train_index], X_score[test_index]
    #     y_train, y_test = Y_height[train_index], Y_height[test_index]
    #     #reg = LassoCV(cv=5, random_state=0,n_jobs=-1).fit(X_train, y_train)
    #     reg = Lasso(0.5).fit(X_train, y_train)
    #     r2_train = reg.score(X_train, y_train)
    #     r2_test = reg.score(X_test, y_test)
    #     print r2_train,r2_test#,reg.alpha_
    #     r2_test_list.append(r2_test)
    #     r2_train_list.append(r2_train)
    # print 'all vars in all loci,all-average r2 training: %.4f, test: %.4f'%(np.mean(r2_test_list),np.mean(r2_train_list))

    # #all vars in top 100 loci
    # print 'top100 loci shape',top_100_score_mat.shape
    # kf = KFold(n_splits=10)
    # r2_train_list = []
    # r2_test_list = []
    # for train_index, test_index in kf.split(top_100_score_mat):
    #     X_train, X_test = top_100_score_mat[train_index], top_100_score_mat[test_index]
    #     y_train, y_test = Y_height[train_index], Y_height[test_index]
    #     #reg = LassoCV(cv=5, random_state=0,n_jobs=-1).fit(X_train, y_train)
    #     reg = Lasso(0.5).fit(X_train, y_train)
    #     r2_train = reg.score(X_train, y_train)
    #     r2_test = reg.score(X_test, y_test)
    #     print r2_train,r2_test#,reg.alpha_
    #     r2_test_list.append(r2_test)
    #     r2_train_list.append(r2_train)
    # print 'all vars in top 100 loci all-average r2 training: %.4f, test: %.4f'%(np.mean(r2_test_list),np.mean(r2_train_list))
    
    #top 20% vars in all loci
    #ratio = 0.2
    
    all_col_list = [item[1] for item in results]
    def get_r_square_with_top(ratio):
        top100_col = []
        for each in all_col_list:
            each_filter =[item for item in list(each) if sum(pd_data_all[item]!=0)>5] 
            #score_list = [var2score[item] for item in list(each)]
            #sorted_score_list = sorted(zip(score_list,list(each)),key=lambda a:a[0],reverse=True)
            score_list = [var2score[item] for item in each_filter]
            sorted_score_list = sorted(zip(score_list,each_filter),key=lambda a:a[0],reverse=True)
            top100_col += [item[1] for item in sorted_score_list[:int(ratio*len(sorted_score_list))]]
        X_score = pd_data_all[top100_col].values
        print 'variants shape',X_score.shape,ratio
        X_score = pd_data_all[top100_col].T.drop_duplicates().T.values
        print 'variants shape',X_score.shape,ratio
        np.save('%s/height_gwas/vars_score_top_%.2f_percent.npy'%(DPATH,ratio),X_score)
        kf = KFold(n_splits=10)
        r2_train_list = []
        r2_test_list = []
        for train_index, test_index in kf.split(X_score):
            X_train, X_test = X_score[train_index], X_score[test_index]
            y_train, y_test = Y_height[train_index], Y_height[test_index]
            #reg = LassoCV(cv=5, random_state=0,n_jobs=-1).fit(X_train, y_train)
            reg = Lasso(0.5).fit(X_train, y_train)
            r2_train = reg.score(X_train, y_train)
            r2_test = reg.score(X_test, y_test)
            #print r2_train,r2_test#,reg.alpha_
            r2_test_list.append(r2_test)
            r2_train_list.append(r2_train)
        print np.mean(r2_train_list),np.mean(r2_test_list),'top ratio%.2f'%ratio
    # get_r_square_with_top(0.01)
    # get_r_square_with_top(0.001)
    # get_r_square_with_top(0.01)
    # get_r_square_with_top(0.05)
    # get_r_square_with_top(0.1)
    # get_r_square_with_top(0.2)
    # get_r_square_with_top(0.3)
    # get_r_square_with_top(0.4)
    # get_r_square_with_top(0.5)
    #get_r_square_with_top(0.6)
    #get_r_square_with_top(0.7)
    #get_r_square_with_top(0.8)
    get_r_square_with_top(0.01)
    get_r_square_with_top(0.02)
    get_r_square_with_top(0.03)
    get_r_square_with_top(0.04)
    get_r_square_with_top(0.05)
    get_r_square_with_top(0.06)
    get_r_square_with_top(0.07)
    get_r_square_with_top(0.08)
    get_r_square_with_top(0.09)
    get_r_square_with_top(0.1)
    #bottom 20% vars in all loci

    def get_r_square_with_bottom(ratio):
        all_col_list = [item[1] for item in results]
        bottom100_col = []
        for each in all_col_list:
            each_filter =[item for item in list(each) if sum(pd_data_all[item]!=0)>5] 
            #score_list = [var2score[item] for item in list(each)]
            #sorted_score_list = sorted(zip(score_list,list(each)),key=lambda a:a[0],reverse=True)
            score_list = [var2score[item] for item in each_filter]
            sorted_score_list = sorted(zip(score_list,each_filter),key=lambda a:a[0],reverse=True)
            bottom100_col += [item[1] for item in sorted_score_list[-int(ratio*len(sorted_score_list)):]]
        X_score = pd_data_all[bottom100_col].values
        print 'variants shape',X_score.shape, ratio
        X_score = pd_data_all[bottom100_col].T.drop_duplicates().T.values
        print 'variants shape',X_score.shape,ratio
        np.save('%s/height_gwas/vars_score_bottom_%.2f_percent.npy'%(DPATH,ratio),X_score)
        kf = KFold(n_splits=10)
        r2_train_list = []
        r2_test_list = []
        for train_index, test_index in kf.split(X_score):
            X_train, X_test = X_score[train_index], X_score[test_index]
            y_train, y_test = Y_height[train_index], Y_height[test_index]
            #reg = LassoCV(cv=5, random_state=0,n_jobs=-1).fit(X_train, y_train)
            reg = Lasso(0.5).fit(X_train, y_train)
            r2_train = reg.score(X_train, y_train)
            r2_test = reg.score(X_test, y_test)
            #print r2_train,r2_test#,reg.alpha_
            r2_test_list.append(r2_test)
            r2_train_list.append(r2_train)
        print np.mean(r2_train_list),np.mean(r2_test_list),'bottom ratio%.2f'%ratio
    #get_r_square_with_bottom(0.1)
    #get_r_square_with_bottom(0.2)
    #get_r_square_with_bottom(0.3)
    #get_r_square_with_bottom(0.4)
    #get_r_square_with_bottom(0.5)
    #get_r_square_with_bottom(0.6)
    #get_r_square_with_bottom(0.7)
    #get_r_square_with_bottom(0.8)
    #get_r_square_with_bottom(0.9)
    #get_r_square_with_bottom(0.9)
    get_r_square_with_bottom(0.91)
    get_r_square_with_bottom(0.92)
    get_r_square_with_bottom(0.93)    
    get_r_square_with_bottom(0.94)
    get_r_square_with_bottom(0.95)  
    get_r_square_with_bottom(0.96)  
    get_r_square_with_bottom(0.97)  
    get_r_square_with_bottom(0.98)  
    get_r_square_with_bottom(0.99)    
    sys.exit()


    #top 20% highest vars
    col_top = [item[1] for item in score_info_sort[-int(0.1*len(pd_data_all.shape[1])):]]
    X_score = pd_data_all[col_top].values
    print X_score.shape, Y_height.shape
    kf = KFold(n_splits=10)
    r2_train_list = []
    r2_test_list = []
    for train_index, test_index in kf.split(X_score):
        X_train, X_test = X_score[train_index], X_score[test_index]
        y_train, y_test = Y_height[train_index], Y_height[test_index]
        #reg = LassoCV(cv=5, random_state=0,n_jobs=-1).fit(X_train, y_train)
        reg = Lasso(0.5).fit(X_train, y_train)
        r2_train = reg.score(X_train, y_train)
        r2_test = reg.score(X_test, y_test)
        print r2_train,r2_test
        r2_test_list.append(r2_test)
        r2_train_list.append(r2_train)
    print 'top100 vars in all loci average r2 training: %.4f, test: %.4f'%(np.mean(r2_test_list),np.mean(r2_train_list))

    #bottom 10% vars
    col_top = [item[1] for item in score_info_sort[:int(0.1*len(pd_data_all.shape[1]))]]
    X_score = pd_data_all[col_top].values
    print X_score.shape, Y_height.shape
    kf = KFold(n_splits=10)
    r2_train_list = []
    r2_test_list = []
    for train_index, test_index in kf.split(X_score):
        X_train, X_test = X_score[train_index], X_score[test_index]
        y_train, y_test = Y_height[train_index], Y_height[test_index]
        #reg = LassoCV(cv=5, random_state=0,n_jobs=-1).fit(X_train, y_train)
        reg = Lasso(0.5).fit(X_train, y_train)
        r2_train = reg.score(X_train, y_train)
        r2_test = reg.score(X_test, y_test)
        print r2_train,r2_test
        r2_test_list.append(r2_test)
        r2_train_list.append(r2_train)
    print 'Bottom-100 average r2 training: %.4f, test: %.4f'%(np.mean(r2_test_list),np.mean(r2_train_list))




    

