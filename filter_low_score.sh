#!/bin/bash
#filter low score instances after overlapping with cseq data
export LC_ALL=C
#current path
PPATH=`readlink -f $BASH_SOURCE | sed -r 's/^\.\///g'`
PPATH=`dirname $PPATH`
DPATH=`echo $PPATH | sed -r 's/prog/data/g'`

echo "PPATH=$PPATH"
echo "DPATH=$DPATH"

declare -A GNMAP
GNMAP=([hg19]=hg19 [hg38]=GRCh38)
gname=${GNMAP['hg19']}

CELL="GM12878"
motif_db="HOMER"
target_abbr="CTCF"
overlap_file="$DPATH/motif_scan/$motif_db/${CELL}_overlap_with_cseq/output_${target_abbr}_overlap.bed"
nonoverlap_file="$DPATH/motif_scan/$motif_db/${CELL}_overlap_with_cseq/output_${target_abbr}_nonoverlap.bed"
echo $overlap_file
echo $nonoverlap_file

awk -F '\t' '$5>10.207' $overlap_file > $DPATH/motif_scan/$motif_db/${CELL}_overlap_with_cseq/output_${target_abbr}_overlap_thred_0.25.bed
awk -F '\t' '$5>10.207' $nonoverlap_file > $DPATH/motif_scan/$motif_db/${CELL}_overlap_with_cseq/output_${target_abbr}_nonoverlap_thred_0.25.bed
