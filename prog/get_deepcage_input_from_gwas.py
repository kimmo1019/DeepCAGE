import numpy as np 
import pandas as pd 
import os, sys
import gzip
import hickle as hkl
from pyfasta import Fasta

SEQ_LEN=1000

def quantile_norm(matrix):
    rank_mean = matrix.stack().groupby(matrix.rank(method='first').stack().astype(int)).mean()
    return matrix.rank(method='min').stack().astype(int).map(rank_mean).unstack()

def one_hot_seq(seq):
    d = {'a':0, 'A':0, 'c':1, 'C':1, 'g':2, 'G':2, 't':3, 'T':3, 'N':4, 'n':4}
    mat = np.zeros((5,SEQ_LEN,1))  
    for i in range(len(seq)):
        mat[d[seq[i]],i,0] = 1
    mat = mat[:4,:,:]
    return mat


def get_data(use_indel):
    f_gwas = open('%s/height_gwas/neighbor_var/%s.log'%(DPATH,rs),'w')
    f_out = open('%s/height_gwas/neighbor_var/%s.fa'%(DPATH,rs),'w')
    #snp direction, tested allele and other allele, sometimes the ref is the former and sometimes the later
    assert genome['chr'+chrom_gwas][int(pos_gwas)-1].upper() == test_allele or genome['chr'+chrom_gwas][int(pos_gwas)-1].upper() == other_allele
    seq_tested = genome['chr'+chrom_gwas][int(pos_gwas)-SEQ_LEN/2:int(pos_gwas)-1]+test_allele+genome['chr'+chrom_gwas][int(pos_gwas):int(pos_gwas)+SEQ_LEN/2]
    seq_other = genome['chr'+chrom_gwas][int(pos_gwas)-SEQ_LEN/2:int(pos_gwas)-1]+other_allele+genome['chr'+chrom_gwas][int(pos_gwas):int(pos_gwas)+SEQ_LEN/2]
    assert len(seq_tested)==SEQ_LEN
    assert len(seq_other)==SEQ_LEN
    f_out.write('>gwas_p_chr%s_%s_%s_%s\n'%(chrom_gwas,pos_gwas,test_allele,other_allele))
    f_out.write(seq_tested)
    f_out.write('\n')
    f_out.write('>gwas_n_chr%s_%s_%s_%s\n'%(chrom_gwas,pos_gwas,test_allele,other_allele))
    f_out.write(seq_other)
    f_out.write('\n')
    vars_list = hkl.load('%s/height_gwas/neighbor_var/%s.hkl'%(DPATH,rs))
    for var in vars_list:
        chrom, pos, ref, alt = var[0], var[1], var[2], var[3]
        if str(pos) == pos_gwas:
            f_gwas.write('%s\t%s\t%s\t%s\n'%(chrom, pos, ref, alt))
        if not use_indel:
            if len(ref)>1 or len(alt)>1:
                continue
        seq_ref = genome['chr'+chrom][int(pos)-SEQ_LEN/2:int(pos)-1]+ref+genome['chr'+chrom][int(pos):int(pos)+SEQ_LEN/2]
        seq_alt = genome['chr'+chrom][int(pos)-SEQ_LEN/2:int(pos)-1]+alt+genome['chr'+chrom][int(pos):int(pos)+SEQ_LEN/2]
        assert len(seq_ref)>=SEQ_LEN
        assert len(seq_alt)>=SEQ_LEN
        f_out.write('>wgs_ref_chr%s_%s_%s_%s\n'%(chrom,pos,ref,alt))
        f_out.write(seq_ref)
        f_out.write('\n')
        f_out.write('>wgs_alt_chr%s_%s_%s_%s\n'%(chrom,pos,ref,alt))
        f_out.write(seq_alt)
        f_out.write('\n')
    f_out.close()
    f_gwas.close()
  
def get_all_data():
    genotype_head_file = '%s/GTEx/GTEx_WGS_genotype_first500.txt' %DPATH
    muscle_info = '%s/GTEx/GTEx_rnaseq_mucsle.txt'%DPATH
    first500 = open(genotype_head_file).readlines()
    line = first500[1].strip().split('\t')
    donor_muscle_ids = [item.split('\t')[0] for item in open(muscle_info).readlines()]
    donor_ids = ['-'.join(item.rstrip().split('-')[:2]) for item in donor_muscle_ids]
    muscle_tissue_ids=[donor_muscle_ids[donor_ids.index(item)] for item in line[9:] if item in donor_ids]

    tf_gexp_file = '%s/GTEx/GTEx_rseq_tpm_peca.csv'%DPATH
    tf_gexp = pd.read_csv(tf_gexp_file ,sep='\t',header=0,index_col=[0])
    tf_gexp_log = np.log(tf_gexp+1)
    tf_gexp_log = pd.DataFrame.transpose(pd.DataFrame.transpose(tf_gexp_log).rank(method='min').stack().astype(int).map(rank_mean).unstack())
    #tf_gexp_log = pd.DataFrame.transpose(quantile_norm(pd.DataFrame.transpose(tf_gexp_log)))


    fasta_file = '%s/height_gwas/neighbor_var/%s.fa'%(DPATH,rs)
    genome = Fasta(fasta_file)
    motifscore_file = '%s/height_gwas/neighbor_var/%s.motifscore.txt'%(DPATH,rs)
    tf_motifscore = pd.read_csv(motifscore_file,sep='\t',header=0,index_col=[0])
    var_info = '%s/height_gwas/neighbor_var/%s.hkl'%(DPATH,rs)
    match_file = '%s/height_gwas/neighbor_var/%s.log'%(DPATH,rs)
    match_line = open(match_file).readlines()
    if len(match_line)>0:
        chrom_gwas, pos_gwas = match_line[0].strip().split('\t')[0],match_line[0].strip().split('\t')[1]
    epsilon = 1e-10
    data_ref_all,data_alt_all=[],[]
    vars_list = hkl.load(var_info)
    for var in vars_list:
        chrom, pos, ref, alt, related_donor_idx = var[0], var[1], var[2], var[3], var[4]
        if len(ref)==1 and len(alt)==1:
            seq_ref_id = 'wgs_ref_chr%s_%s_%s_%s'%(chrom,pos,ref,alt)
            seq_alt_id = 'wgs_alt_chr%s_%s_%s_%s'%(chrom,pos,ref,alt)
            seq_mat_ref = one_hot_seq(genome[seq_ref_id])[np.newaxis,:]
            seq_mat_alt = one_hot_seq(genome[seq_alt_id])[np.newaxis,:]
            motifscore_ref = tf_motifscore.loc[seq_ref_id]
            motifscore_alt = tf_motifscore.loc[seq_alt_id]
            gexp_data_ref,gexp_data_alt =[],[]
            motif_data_ref,motif_data_alt=[],[]
            for i in related_donor_idx:
                gexp = tf_gexp_log.loc[muscle_tissue_ids[i]]
                motif_data_ref.append(motifscore_ref)
                motif_data_alt.append(motifscore_alt)
                gexp_data_ref.append(gexp)
                gexp_data_alt.append(gexp)
            motif_data_ref=np.stack(motif_data_ref)
            motif_data_alt=np.stack(motif_data_alt)
            gexp_data_ref = np.stack(gexp_data_ref)
            gexp_data_alt = np.stack(gexp_data_alt)
            data_ref_all.append([seq_mat_ref,motif_data_ref,gexp_data_ref])
            data_alt_all.append([seq_mat_alt,motif_data_alt,gexp_data_alt])
    hkl.dump([data_ref_all,data_alt_all],'%s/height_gwas/neighbor_var/%s.ref_alt.hkl'%(DPATH,rs))

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
    #pd_data.to_csv('%s/height_gwas/neighbor_var/%s.causal_score.csv'%(DPATH,rs),sep='\t')

    columns_filter = [item for item in columns if np.max(pd_data[item])>1]
    if len(columns_filter)>0:
        pd_data_filter = pd_data[columns_filter]
        pd_data_filter.to_csv('%s/height_gwas/neighbor_var/%s.causal_score_thred1.csv'%(DPATH,rs),sep='\t')

if __name__=="__main__":
    DPATH='/home/liuqiao/software/DeepCAGE/data'
    height_GWAS='%s/height_gwas/Meta-analysis_Wood_et_al+UKBiobank_2018_top_3290_from_COJO_analysis.txt'%DPATH
    GWAS_variants = '%s/height_gwas/neighbor_var'%DPATH
    tf_gexp_train_file = '%s/encode/preprocessed/tf_gexp.csv'%DPATH
    tf_gexp_train = pd.read_csv(tf_gexp_train_file ,sep='\t',header=0,index_col=[0])
    tf_gexp_train_log = np.log(tf_gexp_train+1)
    rank_mean = pd.DataFrame.transpose(tf_gexp_train_log).stack().groupby(pd.DataFrame.transpose(tf_gexp_train_log).rank(method='first').stack().astype(int)).mean()
    #tf_gexp_train_log = pd.DataFrame.transpose(quantile_norm())

    rs, chrom_gwas, pos_gwas, test_allele, other_allele = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    #get_data(use_indel=False)#generate fasta
    get_all_data()#seq_mat, gexp, and motifscore
    #get_causal_score(rs)



    

