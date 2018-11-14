#!/bin/bash
#generate negative set (random peak) for dnase peak, calculate the overlap regions with motifs regions(from scan),then sort them
export LC_ALL=C
PPATH=`readlink -f $BASH_SOURCE | sed -r 's/^\.\///g'`
PPATH=`dirname $PPATH`
DPATH=`echo $PPATH | sed -r 's/prog/data/g'`
echo $DPATH
CELL=$1
Database='HOMER'
#generate random peaks against DNase peaks
#file="$DPATH/$CELL/dseq/wgEncodeAwgDnaseUwdukeGm12878UniPk.narrowPeak.gz"
for exp in `ls $DPATH/all/$CELL/dseq`
do
	for file in `ls $DPATH/all/$CELL/dseq/$exp`
	do
		if [ "${file##*.}" = "gz" ]; then
			echo $DPATH/all/$CELL/dseq/$exp/$file
			zcat $DPATH/all/$CELL/dseq/$exp/$file > $DPATH/all/$CELL/dseq/$exp/positive.bed
			bedtools shuffle -chrom -noOverlapping -excl $DPATH/all/$CELL/dseq/$exp/positive.bed -i $DPATH/all/$CELL/dseq/$exp/positive.bed -g /home/liuqiao/bedtools2/genomes/human.hg19.genome > $DPATH/all/$CELL/dseq/$exp/negative.bed
			#index each region
			awk -F '\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,NR,$5,$6,$7}' $DPATH/all/$CELL/dseq/$exp/positive.bed > $DPATH/all/$CELL/dseq/$exp/positive_indexd.bed
			awk -F '\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,NR,$5,$6,$7}' $DPATH/all/$CELL/dseq/$exp/negative.bed > $DPATH/all/$CELL/dseq/$exp/negative_indexd.bed
			#find overlapping regions with motifs binding sites (at least 1 bp)
			bedtools intersect  -wa -a $DPATH/motif_scan/$Database/output.bed -wb -b $DPATH/all/$CELL/dseq/$exp/positive_indexd.bed > $DPATH/all/$CELL/dseq/$exp/${Database}_overlap_pos.bed
			bedtools intersect  -wa -a $DPATH/motif_scan/$Database/output.bed -wb -b $DPATH/all/$CELL/dseq/$exp/negative_indexd.bed > $DPATH/all/$CELL/dseq/$exp/${Database}_overlap_neg.bed
			#filter the regions in chrY, sort the regions by region index
			awk -F '\t' '$7!="chrY" {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$7,$8,$9,$10,$13,$1,$2,$3,$4,$5,$6}' $DPATH/all/$CELL/dseq/$exp/${Database}_overlap_pos.bed | sort -n -k 4 -o  $DPATH/all/$CELL/dseq/$exp/${Database}_overlap_filter_sorted_pos.bed
			awk -F '\t' '$7!="chrY" {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$7,$8,$9,$10,$13,$1,$2,$3,$4,$5,$6}' $DPATH/all/$CELL/dseq/$exp/${Database}_overlap_neg.bed | sort -n -k 4 -o  $DPATH/all/$CELL/dseq/$exp/${Database}_overlap_filter_sorted_neg.bed
			rm -f $DPATH/all/$CELL/dseq/$exp/${Database}_overlap_pos.bed
			rm -f $DPATH/all/$CELL/dseq/$exp/${Database}_overlap_neg.bed
			rm -f $DPATH/all/$CELL/dseq/$exp/positive.bed
			rm -f $DPATH/all/$CELL/dseq/$exp/negative.bed
		fi
	done
done


