import sys


usage='''
Usage: python 3.1.Generate_tf_exp.py  [selected cell types] [output]
-- a program for generating transcription factors expression matrix(N x C)
OPTIONS:
    [selected cell types] -- selected cell types defined in 3.0
    [output]  --  transcription factors expression matrix as a plain text
'''

if len(sys.argv) != 3:
    print usage
    exit(0)

DPATH='../data'

#download motif data (769 motifs) from http://hocomoco11.autosome.ru/HUMAN/mono?full=true
motif_db = '%s/motif_db/HOCOMOCO/HUMAN_mono_motifs_full.tsv'%DPATH
motif_file = '%s/motif_db/HOCOMOCO/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.001.motif'%DPATH
motif2tf={item.split()[0]:item.split()[1].rstrip(';') for item in open(motif_db).readlines()[1:]}
print '%d motif models and %d TFs found'%(len(motif2tf.keys()),len(set(motif2tf.values())))
#load motifs
motifs=[]
f_motif = open(motif_file)
line = f_motif.readline()
while not line=='':
    if line[0]=='>':
        motifs.append(line.split()[1])
    line = f_motif.readline()
f_motif.close()

tfs = [motif2tf[item] for item in motifs]
#print tfs[:10]


#load the selected cell types
celltypes = [line.rstrip() for line in open(sys.argv[1]).readlines()]

gexps = {}
for c in celltypes:
    gexps[c] = {line.split()[0]:line.rstrip().split()[1] for line in open('../data/processed_RNA_DNase/%s.TPM.tsv'%(c)).readlines()}

fd = open(sys.argv[2],'w')
fd.write('#cell\t'+'\t'.join(tfs)+'\n')
for c in celltypes:
    #fd.write(c + '\t' + '\t'.join([gexps[c].get(tf,'0') for tf in tfs])+'\n')
    fd.write(c + '\t' + '\t'.join([gexps[c][tf] for tf in tfs])+'\n')
fd.close()
