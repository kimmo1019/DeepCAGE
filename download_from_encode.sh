#!/bin/bash
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

mkdir -p "$DPATH/$CELL"
Rseq="$DPATH/$CELL/rseq"
Cseq="$DPATH/$CELL/cseq"
Dseq="$DPATH/$CELL/dseq"
#meta file path
mpath="$DPATH/20171106.txt.gz"

mkdir -p $Rseq
mkdir -p $Cseq
mkdir -p $Dseq

#downloading RNA-seq data
zcat -f $mpath		|
tail -n +2			|
grep -P "$gname"    |
awk -F '\t' '$5=="RNA-seq"' |
awk -F '\t' '$3=="gene quantifications"' |
#remove control and histone marker (chip-seq)
#grep -v -P "control" |
#grep -v -P "\sH\d+"  |
grep -P "tsv"		|
grep -P "$CELL"     |
#sort -d	-k1,1	|
awk  -F "\t"   	'
{
	print $1 "\t" $42 "\t" $40 "\t" $4
}' |
while read line; do
	item=(${line//\\t/})
	facc=${item[0]}
	furl=${item[1]}
	msum=${item[2]}
	exp=${item[3]}
	fbed=$Rseq/$exp/$facc.tsv
	fmd5=$Rseq/$exp/$facc.md5
	printf "[Downloading] %s\t%s\t%s\t%s\n" $facc $furl $msum $name
	`
	for ((i = 0; i < 5; i++)); do
		if [[ ! -e $fmd5 || \`cat $fmd5\` != $msum ]]; then
			mkdir -p \`dirname $fbed\`
			curl  -o $fbed -L $furl -C - -s
			if [  -e $fbed ]; then
				md5sum $fbed | cut -d ' ' -f 1 > $fmd5
			fi
		fi
	done
	` &
	sleep 1s; while [ `ps -T | grep -P "\s+curl$" | wc -l` -ge 5 ]; do sleep 1s; done
done
echo RNA-seq done

#downloading Dnase-seq data
zcat -f $mpath		|
tail -n +2			|
grep -P "$gname"	|
awk -F '\t' '$5=="DNase-seq"' |
awk -F '\t' '$2=="bed narrowPeak"' |
#remove control and histone marker (chip-seq)
#grep -v -P "control" |
#grep -v -P "\sH\d+"  |
grep -P "$CELL"		  |
#sort -d	-k1,1	  |
awk  -F "\t"   	'
{
	print $1 "\t" $42 "\t" $40 "\t" $4
}' |
while read line; do
	item=(${line//\\t/})
	facc=${item[0]}
	furl=${item[1]}
	msum=${item[2]}
	exp=${item[3]}
	fbed=$Dseq/$exp/$facc.bed.gz
	fmd5=$Dseq/$exp/$facc.md5
	printf "[Downloading] %s\t%s\t%s\t%s\n" $facc $furl $msum $exp
	`
	for ((i = 0; i < 5; i++)); do
		if [[ ! -e $fmd5 || \`cat $fmd5\` != $msum ]]; then
			mkdir -p \`dirname $fbed\`
			printf "[Downloading] %s\t%s\t%s\t%s\n" $facc $furl $msum $name
			curl  -o $fbed -L $furl -C - -s
			if [  -e $fbed ]; then
				md5sum $fbed | cut -d ' ' -f 1 > $fmd5
			fi
		fi
	done
	` &  #backend (multi-process)
	sleep 1s; while [ `ps -T | grep -P "\s+curl$" | wc -l` -ge 5 ]; do sleep 1s; done
done
echo DNase-seq done

#downloading ChIP-seq data
zcat -f $mpath		|
tail -n +2			|
grep -P "$gname"	|
awk -F '\t' '$5=="ChIP-seq"' |
awk -F '\t' '$2=="bed narrowPeak"' |
#remove control and histone marker (chip-seq)
grep -v -P "eGFP"     |
grep -v -P "\sH\d+\w" |
grep -P "$CELL"		  |
#sort -d	-k1,1	  |
awk  -F "\t"   	'
{
	print $1 "\t" $42 "\t" $40 "\t" $4 "\t" $17
}' |
while read line; do
	item=(${line//\\t/})
	facc=${item[0]}
	furl=${item[1]}
	msum=${item[2]}
	exp=${item[3]}
	target=${item[4]}
	fbed=$Cseq/$target/$exp/$facc.bed.gz
	fmd5=$Cseq/$target/$exp/$facc.md5
	printf "[Downloading] %s\t%s\t%s\t%s\n" $facc $furl $msum $exp
	`
	for ((i = 0; i < 5; i++)); do
		if [[ ! -e $fmd5 || \`cat $fmd5\` != $msum ]]; then
			mkdir -p \`dirname $fbed\`
			printf "[Downloading] %s\t%s\t%s\t%s\n" $facc $furl $msum $name
			curl  -o $fbed -L $furl -C - -s
			if [  -e $fbed ]; then
				md5sum $fbed | cut -d ' ' -f 1 > $fmd5
			fi
		fi
	done
	` &  #backend (multi-process)
	sleep 1s; while [ `ps -T | grep -P "\s+curl$" | wc -l` -ge 10 ]; do sleep 1s; done
done
echo ChIP-seq done




:<<\COMMENT
1	File accession
2	File format
3	Output type
4	Experiment accession
5	Assay
6	Biosample term id
7	Biosample term name
8	Biosample type
9	Biosample life stage
10	Biosample sex
11	Biosample Age
12	Biosample organism
13	Biosample treatments
14	Biosample subcellular fraction term name
15	Biosample phase
16	Biosample synchronization stage
17	Experiment target
18	Antibody accession
19	Library made from
20	Library depleted in
21	Library extraction method
22	Library lysis method
23	Library crosslinking method
24	Library strand specific
25	Experiment date released
26	Project
27	RBNS protein concentration
28	Library fragmentation method
29	Library size range
30	Biological replicate(s)
31	Technical replicate
32	Read length
33	Mapped read length
34	Run type
35	Paired end
36	Paired with
37	Derived from
38	Size
39	Lab
40	md5sum
41	dbxrefs
42	File download URL
43	Assembly
44	Platform
45	Controlled by
46	File Status
47	Audit WARNING
48	Audit INTERNAL_ACTION
49	Audit NOT_COMPLIANT
50	Audit ERROR
COMMENT
