#!/bin/bash
export LC_ALL=C
PPATH=`readlink -f $BASH_SOURCE | sed -r 's/^\.\///g'`
PPATH=`dirname $PPATH`
DPATH=`echo $PPATH | sed -r 's/prog/data/g'`
echo $DPATH
zcat $DPATH/20180113.txt.gz | awk -F '\t' '$5=="RNA-seq"'|awk -F '\t' '$3=="gene quantifications"'|awk -F '\t' '$46=="released"'|grep hg19|grep -P 'tsv'|sort -k 7 -k 5  -t $'\t'|cut -f 1-4,7,26 > $DPATH/rseq.txt
zcat $DPATH/20180113.txt.gz | grep hg19|awk -F '\t' '$5=="DNase-seq"'|awk -F '\t' '$2=="bed narrowPeak"'|awk -F '\t' '$46=="released"'|sort -k 7 -k 5  -t $'\t'|cut -f 1-4,7,26  > $DPATH/dseq.txt
#remove eGFP(荧光蛋白) target and hisone markers 
zcat $DPATH/20180113.txt.gz | grep hg19|awk -F '\t' '$5=="ChIP-seq"' |awk -F '\t' '$2=="bed narrowPeak"'|awk -F '\t' '$46=="released"'|grep -v -P "eGFP"|grep -v -P "\sH\d+\w" |sort -k 7 -k 8  -t $'\t'|cut -f 1-4,7,17,26 > $DPATH/cseq.txt
