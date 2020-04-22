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
    GWAS_snps = [ item.split('\t')[:5] for item in open(height_GWAS).readlines()[1:] ]
    count=0
    for each in GWAS_snps:
        count+=1
        if count%100==0:
            print count
        rs, chrom, pos, test_allele, other_allele = each
        f_out = open('%s/height_gwas/neighbor_var/%s.fa'%(DPATH,rs),'w')
        #snp direction, tested allele and other allele, sometimes the ref is the former and sometimes the later
        assert genome['chr'+chrom][int(pos)-1].upper() == test_allele or genome['chr'+chrom][int(pos)-1].upper() == other_allele
        seq_tested = genome['chr'+chrom][int(pos)-SEQ_LEN/2:int(pos)-1]+test_allele+genome['chr'+chrom][int(pos):int(pos)+SEQ_LEN/2]
        seq_other = genome['chr'+chrom][int(pos)-SEQ_LEN/2:int(pos)-1]+other_allele+genome['chr'+chrom][int(pos):int(pos)+SEQ_LEN/2]
        assert len(seq_tested)==SEQ_LEN
        assert len(seq_other)==SEQ_LEN
        f_out.write('>gwas_p_chr%s_%s\n'%(chrom,pos))
        f_out.write(seq_tested)
        f_out.write('\n')
        f_out.write('>gwas_n_chr%s_%s\n'%(chrom,pos))
        f_out.write(seq_other)
        f_out.write('\n')
        vars_list = hkl.load('%s/height_gwas/neighbor_var/%s.hkl'%(DPATH,rs))
        for var in vars_list:
            chrom, pos, ref, alt = var[0], var[1], var[2], var[3]
            if not use_indel:
                if len(ref)>1 or len(alt)>1:
                    continue
            seq_ref = genome['chr'+chrom][int(pos)-SEQ_LEN/2:int(pos)-1]+ref+genome['chr'+chrom][int(pos):int(pos)+SEQ_LEN/2]
            seq_alt = genome['chr'+chrom][int(pos)-SEQ_LEN/2:int(pos)-1]+alt+genome['chr'+chrom][int(pos):int(pos)+SEQ_LEN/2]
            assert len(seq_ref)>=SEQ_LEN
            assert len(seq_alt)>=SEQ_LEN
            f_out.write('>wgs_ref_chr%s_%s\n'%(chrom,pos))
            f_out.write(seq_ref)
            f_out.write('\n')
            f_out.write('>wgs_alt_chr%s_%s\n'%(chrom,pos))
            f_out.write(seq_alt)
            f_out.write('\n')
        f_out.close()
            
            

            

        

        



        
        

        
        

if __name__=="__main__":
    DPATH='/home/liuqiao/software/DeepCAGE/data'
    genome = Fasta('/home/liuqiao/openness_pre/genome.fa')
    height_GWAS='%s/height_gwas/Meta-analysis_Wood_et_al+UKBiobank_2018_top_3290_from_COJO_analysis.txt'%DPATH
    GWAS_variants = '%s/height_gwas/neighbor_var'%DPATH
    get_data(use_indel=False)

    #deepcage input
    # X_test_mat = np.zeros((nb_regions,4,SEQ_LEN,1))
    # X_test_gexp = np.zeros((nb_regions,NUM_TF))
    # X_test_motif = np.zeros((nb_regions,NUM_TF))