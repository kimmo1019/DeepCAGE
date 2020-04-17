#!/bin/bash
#generate chip-seq, dnase-seq, rna-seq from the same cell line
#usage: bash generate_data.sh cellline_id (e.g. bash generate_data.sh 1)
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
Rseq_exp=$1
Dseq_exp=$2

Rseq="$DPATH/encode/rseq/$1"
Dseq="$DPATH/encode/dseq/$2"
#meta file path
mpath="$DPATH/encode/metadata.tsv"

mkdir -p $Rseq
mkdir -p $Dseq

#downloading RNA-seq data
    cat  $mpath		|
    tail -n +2			|
    grep -P "$gname"    |
    awk -F '\t' '$4=="'$Rseq_exp'"' |
    awk -F '\t' '$5=="RNA-seq"' |
    awk -F '\t' '$3=="gene quantifications"' |
    awk -F '\t' '$48=="released"' |
    #remove control and histone marker (chip-seq)
    #grep -v -P "control" |
    #grep -v -P "\sH\d+"  |
    grep -P "tsv"		|
    #sort -d	-k1,1	|
    awk  -F "\t"   	'
    {
        print $1 "\t" $43 "\t" $41 "\t" $4
    }' |
    while read line; do
        item=(${line//\\t/})
        facc=${item[0]}
        furl=${item[1]}
        msum=${item[2]}
        exp=${item[3]}
        fbed=$Rseq/$facc.tsv
        fmd5=$Rseq/$facc.md5
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
    echo RNA-seq .tsv download finished


    #downloading Dnase-seq data(peak)
    cat  $mpath		|
    tail -n +2			|
    grep -P "$gname"	|
    awk -F '\t' '$4=="'$Dseq_exp'"' |
    awk -F '\t' '$5=="DNase-seq"' |
    awk -F '\t' '$2=="bed narrowPeak"' |
    awk -F '\t' '$48=="released"' |
    #remove control and histone marker (chip-seq)
    #grep -v -P "control" |
    #grep -v -P "\sH\d+"  |
    #sort -d	-k1,1	  |
    awk  -F "\t"   	'
    {
        print $1 "\t" $43 "\t" $41 "\t" $4
    }' |
    while read line; do
        item=(${line//\\t/})
        facc=${item[0]}
        furl=${item[1]}
        msum=${item[2]}
        exp=${item[3]}
        fbed=$Dseq/$facc.bed.gz
        fmd5=$Dseq/$facc.md5
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
    echo DNase-seq .narrowPeak download finished


cat  $mpath		|
tail -n +2		|
grep -P "$gname"	|
grep DNase-seq |
grep bam |
awk -F '\t' '$4=="'$Dseq_exp'"' |
awk -F '\t' '$3=="alignments"' |
#awk -F '\t' '$30=="1" || $30=="2"' |
#awk -F '\t' '$30=="1"' |
awk -F '\t' '$48=="released"' |
#sort -d	-k1,1	  |
awk  -F "\t"   	'
{
	print $1 "\t" $43 "\t" $41 "\t" $4
}' |
while read line; do
	item=(${line//\\t/})
	facc=${item[0]}
	furl=${item[1]}
	msum=${item[2]}
	exp=${item[3]}
	fbed=$Dseq/$facc.bam
	fmd5=$Dseq/$facc.md5
	printf "[Downloading] %s\t%s\t%s\t%s\n" $facc $furl $msum $exp
	for ((i = 0; i < 1; i++)); do
     cat $fmd5
     echo $msum $fbed $furl
		if [[ ! -e $fmd5 || `cat $fmd5` != $msum ]]; then
			mkdir -p `dirname $fbed`
			printf "[Downloading] %s\t%s\t%s\t%s\n" $facc $furl $msum $name
			#curl  -o $fbed -L $furl -C - -s
         curl  -o $fbed -L $furl 
			if [  -e $fbed ]; then
				md5sum $fbed | cut -d ' ' -f 1 > $fmd5
			fi
		fi
	done
done
    echo DNase-seq  .bam download finished
