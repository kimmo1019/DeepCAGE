#!/bin/bash
#check the motif scan with cseq data(gold standard)
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


#meta file path
mpath="$DPATH/ucsc/TfbsUniPk/files.txt"

CELL="K562"
motif_db="HOMER"

#motif files
motif_file="$DPATH/motif_db/$motif_db/motifs_list"

cat $mpath	|
while read line; do
	item=(${line//;/ })
	name=${item[0]}
	cell=$(echo ${item[6]}|cut -d '=' -f 2)
	target=$(echo ${item[8]}|cut -d '=' -f 2)
	msum=$(echo ${item[19]}|cut -d '=' -f 2)
	if [ $cell = $CELL ] ; then
		#remove "_()" in target name e.g. ELF1_(SC-631)
		target_abbr=$(echo $target|cut -d '_' -f 1) 
		num=$(cat $motif_file|grep -i -P "^$target_abbr\(" | wc -l)
		if [ $num -gt 0 ];then
			matched_target=$(cat $motif_file|grep -i -P "^$target_abbr\(")
			printf "overlap %s with target %s\n" $name $matched_target
			#to do: filter output.bed with a score threshold
			#note that if use perl regular expression grep 'CTCF(Zf)' should be 'CTCF\(Zf\)'
			grep  "$matched_target" $DPATH/motif_scan/$motif_db/output.bed > $DPATH/motif_scan/$motif_db/output_$target_abbr.bed
			zcat $DPATH/ucsc/TfbsUniPk/$cell/$name > $DPATH/ucsc/TfbsUniPk/$cell/Tfbs_peak.bed
			cseq_file="$DPATH/ucsc/TfbsUniPk/$cell/Tfbs_peak.bed"
			scan_file="$DPATH/motif_scan/$motif_db/output_$target_abbr.bed"
			#todo calculate the coverage rate
			num_cseq=$(wc -l "$cseq_file"|awk '{print $1}')
			num_scan=$(wc -l "$scan_file"|awk '{print $1}')
			bedtools intersect -wa -a $scan_file -b $cseq_file|sort|uniq  > $DPATH/motif_scan/$motif_db/output_${target_abbr}_overlap.bed
			bedtools subtract -a $DPATH/motif_scan/$motif_db/output_$target_abbr.bed -b $DPATH/motif_scan/$motif_db/output_${target_abbr}_overlap.bed > $DPATH/motif_scan/$motif_db/output_${target_abbr}_nonoverlap.bed
			TP1=$(bedtools intersect -wa -a $scan_file -b $cseq_file|sort|uniq|wc -l)
			TP2=$(bedtools intersect -wa -a $cseq_file -b $scan_file|sort|uniq|wc -l)
			printf "cseq peaks: %s\toverlapped peaks: %s\tmotifs scan sites: %s\toverlapped sites: %s\n" $num_cseq $TP2 $num_scan $TP1
			printf "%s\t%s\t%s\t%s\t%s\n" $matched_target $num_cseq $TP2 $num_scan $TP1 >> $DPATH/$CELL/cseq/Tfbs_overlap_with_cseq.txt
			rm $DPATH/motif_scan/$motif_db/output_$target_abbr.bed
			rm $DPATH/ucsc/TfbsUniPk/$cell/Tfbs_peak.bed
		fi
	fi

done


