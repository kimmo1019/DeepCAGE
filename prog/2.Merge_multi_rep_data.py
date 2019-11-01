import sys
import os

usage='''
Usage: python 2.Merge_multi_rep_data.py [CELLID]
-- a program for merging multiple replicates data from RNA-seq and DNase-seq data
[CELLID] : pre-defined cell ID (from 1 to 55)
'''
if len(sys.argv)!=2:
    print usage
    sys.exit(1)


DPATH='../data'
#download motif data (769 motifs) from http://hocomoco11.autosome.ru/HUMAN/mono?full=true
motif_db = '%s/motif_db/HOCOMOCO/HUMAN_mono_motifs_full.tsv'

#download refgene data from https://www.gencodegenes.org/human/release_19.html
id2name={line.split('\t')[0]:line.strip().split('\t')[1] for line in open('%s/geneid2name.txt'%DPATH).readlines()}

cell_id=sys.argv[1]
#rna multiple replicates merge
numRNAs=0
mergeRNA={}
for exp in os.listdir('%s/all/%s/rseq'%(DPATH,cell_id)):
    for rna_file in os.listdir('%s/all/%s/rseq/%s'%(DPATH,cell_id,exp)):
        if rna_file[-3:]=='tsv':
            for line in open('%s/all/%s/rseq/%s/%s'%(DPATH,cell_id,exp,rna_file)).readlines()[1:]:
                gene_id = line.split('\t')[0].split('.')[0]#ENG...
                gene_tpm = float(line.split('\t')[5])#0.25
                if id2name.has_key(gene_id):
                    mergeRNA[id2name[gene_id]] = mergeRNA.get(id2name[gene_id],0.0)+gene_tpm
                else:
                    print gene_id
            #sys.exit()
            numRNAs+=1
cmd="mkdir -p %s/processed_RNA_DNase"%DPATH
os.system(cmd) 
f_out = open('%s/processed_RNA_DNase/%s.TPM.tsv'%(DPATH,cell_id),'w')
for gene in  mergeRNA.keys():
    mergeRNA[gene] /= numRNAs
    f_out.write(gene+'\t'+str(mergeRNA[gene])+'\n')
f_out.close()

#dnase mutiple replicates merge
numDNase=0
cmd  = 'bedtools intersect -wa -f 0.5 -c -a %s/whole_genome_200bp.bed -b'%DPATH
for exp in os.listdir('%s/all/%s/dseq'%(DPATH,cell_id)):
    for dnase_file in os.listdir('%s/all/%s/dseq/%s'%(DPATH,cell_id,exp)):
        if dnase_file[-2:]=='gz':
            cmd += ' %s/all/%s/dseq/%s/%s'%(DPATH,cell_id,exp,dnase_file)
            numDNase+=1
cmd += " | awk '$4>=%d{print $1\"\\t\"$2\"\\t\"$3}' > %s/processed_RNA_DNase/%s.peak.bins.bed" % ((numDNase+1)/2, DPATH ,cell_id)
os.system(cmd)  
