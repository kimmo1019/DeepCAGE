#!/bin/bash
#Select candidate/potential regulatory locus
export LC_ALL=C
PPATH=`readlink -f $BASH_SOURCE | sed -r 's/^\.\///g'`
PPATH=`dirname $PPATH`
DPATH=`echo $PPATH | sed -r 's/prog/data/g'`

cd $DPATH/processed_RNA_DNase


# select cell types with peaks 50k ~ 200k
wc -l *peak.bins.bed | awk '$1>=50000&&$$1<=200000{print $2}'|grep peak.bins.bed| awk '{split($1,a,".");print a[1]}' | sort -k1 -n > celltypes.filtered.txt

# take union of peaks that appear at >= 2 cell types remove chrX chrY and chrM
wc -l *peak.bins.bed | awk '$1>=50000&&$$1<=200000{print $2}'|grep peak.bins.bed| xargs cat | sort | uniq -c | awk '$1>1&&$2!="chrX"&&$2!="chrY"&&$2!="chrM"{print $2"\t"$3"\t"$4}' |bedtools sort -i | > union.peaks.bed

# generate padded (up and down 400bp, 1kb in total) bed and fasta
awk '{print $1"\t"$2-400"\t"$3+400}' union.peaks.bed > union.peaks.pad1k.bed
#outrange chr17	81194400	81195400,delete it manuly
#also delete chr17 81194800 81195000 in union.peaks.bed
#finally got 1750033 peaks
bedtools getfasta -fi ../male.hg19.fa -bed union.peaks.pad1k.bed -fo union.peaks.pad1k.fa
