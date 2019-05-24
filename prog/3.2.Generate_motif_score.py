import sys



usage='''Usage: python 3.2.Generate_motif_score.py  [peak file] [motif file] [output]
-- a program for generating motif score matrix(L x N)
OPTIONS:
    [peak file] -- union peak with padding in 3.0
    [motif file]  --  HOCOMOCO motif file (769 motifs)
'''
if len(sys.argv) != 4:
    print usage
    exit(0)
    

peaks=[item.split()[0]+':'+item.split()[1]+'-'+item.rstrip().split()[2] for item in open(sys.argv[1]).readlines()]
#print peaks[:10]
#load motif
motifs=[]
f_motif = open(sys.argv[2])
line = f_motif.readline()
while not line=='':
    if line[0]=='>':
        motifs.append(line.split()[1])
    line = f_motif.readline()
f_motif.close()

print "%d motifs were loaded.\n"%len(motifs)
#print(motifs[:10])

#load motifscan outcome from motifmatchr,f1 denotes motifs 1-250, f2 denotes motifs 251-500, f3 denotes motifs 501-769 
f1 = open('../data/motif_db/HOCOMOCO/motif_score_mat_part1.txt')
f2 = open('../data/motif_db/HOCOMOCO/motif_score_mat_part2.txt')
f3 = open('../data/motif_db/HOCOMOCO/motif_score_mat_part3.txt')
#f_out = open('../data/motif_db/HOCOMOCO/motif_score_mat.txt','w')
f_out = open(sys.argv[3],'w')
f_out.write('#region\t'+'\t'.join(motifs)+'\n')
#skip frist line
line1 = f1.readline()
line2 = f2.readline()
line3 = f3.readline()

line1 = f1.readline()
line2 = f2.readline()
line3 = f3.readline()
idx=0
while line1!='':
    f_out.write(peaks[idx]+'\t')
    score1_list = line1.rstrip().split(' ')[1:]
    score2_list = line2.rstrip().split(' ')[1:]
    score3_list = line3.rstrip().split(' ')[1:]
    score_all_list = score1_list+score2_list+score3_list
    f_out.write('\t'.join(score_all_list)+'\n')
    line1 = f1.readline()
    line2 = f2.readline()
    line3 = f3.readline()
    idx+=1
f1.close()
f2.close()
f3.close()
f_out.close()
    
    
    
