import numpy as np 
import pandas as pd 
import os, sys
import gzip
from pyfasta import Fasta
import hickle as hkl
import pickle
from multiprocessing.dummy import Pool as ThreadPool

def quantile_norm(matrix):
    rank_mean = matrix.stack().groupby(matrix.rank(method='first').stack().astype(int)).mean()
    return matrix.rank(method='min').stack().astype(int).map(rank_mean).unstack()

def extract_tf_expr():
    tf_gexp = pd.read_csv(tf_gexp_file ,sep='\t',header=0,index_col=[0])
    #tf_gexp_log = np.log(tf_gexp+1)
    #tf_gexp_log = pd.DataFrame.transpose(quantile_norm(pd.DataFrame.transpose(tf_gexp_log)))
    #1672 tf from fantom database
    tf_names = list(tf_gexp.columns)
    #491 gtex sample id from muscle tissue
    sample_ids = [item.split('\t')[0] for item in open(GTEx_sample_info).readlines() if item.strip().split('\t')[1]=='Muscle']
    data_gtex = np.zeros((len(sample_ids),len(tf_names)))
    f_rseq = gzip.open(GTEx_rseq_file)
    line = f_rseq.readline()
    line = f_rseq.readline()
    line = f_rseq.readline()
    all_sample_ids = line.strip().split('\t')[2:]
    idx_sample = [all_sample_ids.index(item) for item in sample_ids]
    tf_count=0
    line_count = 3
    line = f_rseq.readline()
    while line != "":
        line_count+=1
        if len(line.split('\t'))==1:
            print line,line_count
            sys.exit(0)
        tf = line.split('\t')[1]
        if tf in tf_names:
            tf_count+=1
            col_idx = tf_names.index(tf)
            tf_exprs = [float(item) for item in line.strip().split('\t')[2:]]
            for i in range(len(sample_ids)):
                data_gtex[i][col_idx] = tf_exprs[idx_sample[i]]
        line = f_rseq.readline()

    print tf_count
    #1643 out of 1672 tfs can be found expression
    pd_gtex_rseq = pd.DataFrame(data_gtex,index=sample_ids,columns=tf_names)
    pd_gtex_rseq.to_csv('%s/GTEx/GTEx_rseq_tpm.csv'%DPATH,sep='\t')

def find_wgs_vars(snp):
    #chrom, pos, ref, alt = snp[0], int(snp[1]), snp[2], snp[3]
    f_wgs = open(GTEx_genotype)
    line = f_wgs.readline()
    flag = 0
    records=[]
    while line != "":
        if line[:len(snp[0])+1] == snp[0]+' ':
            elements = line.strip().split(' ')
            chrom, pos, ref, alt = elements[0], int(elements[1]), elements[2], elements[3]
            if abs(pos-int(snp[1])) <= Max_range:
                ind_var = [i for i, item in enumerate(elements[4:]) if item=='1']
                if len(ind_var)>=1:
                    records.append([chrom, pos, ref, alt, ind_var])
            flag = 1
        else:
            if flag==1:
                break
        line = f_wgs.readline()
    return records


def process_WGS():
    #remove chrX and chrY
    GWAS_snps = [item.split('\t')[1:5] for item in open(height_GWAS).readlines()[1:] ]
    results = map(find_wgs_vars,GWAS_snps[0:1])
    pickle.dump(results,open('var.hkl','wb'))
    pool = ThreadPool()
    #results = pool.map(find_wgs_vars,GWAS_snps[0:2])
    pool.close()
    pool.join()

if __name__=="__main__":
    DPATH='/home/liuqiao/software/DeepCAGE/data'
    genome = Fasta('/home/liuqiao/openness_pre/genome.fa')
    tf_gexp_file = '%s/encode/preprocessed/fantom_tf_gexp.csv'%DPATH
    GTEx_rseq_file='%s/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz'%DPATH
    GTEx_sample_info='%s/GTEx/GTEx_rnaseq_names.txt'%DPATH
    GTEx_genotype='%s/GTEx/GTEx_WGS_vcf.txt'%DPATH
    height_GWAS='%s/height_gwas/Meta-analysis_Wood_et_al+UKBiobank_2018_top_3290_from_COJO_analysis.txt'%DPATH
    Max_range = 100000 #select variants within 100k
    rs, chrom, pos, ref, alt = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    records = find_wgs_vars([chrom,int(pos), ref, alt])
    hkl.dump(records,'%s/height_gwas/neighbor_var/%s.hkl'%(DPATH,rs))
    #extract_tf_expr()
    #process_WGS()
    