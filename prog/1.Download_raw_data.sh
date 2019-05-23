#!/bin/bash
export LC_ALL=C
PPATH=`readlink -f $BASH_SOURCE | sed -r 's/^\.\///g'`
PPATH=`dirname $PPATH`
DPATH=`echo $PPATH | sed -r 's/prog/data/g'`
declare -A GNMAP
GNMAP=([hg19]=hg19 [hg38]=GRCh38)
gname=${GNMAP['hg19']}
#meta data from encode
mpath="$DPATH/20180113.txt.gz"

function printHelp()
{
    echo "Usage: bash $(basename "$0") [-c CELLID] [-h] [-r] [-p] [-b] "
    echo "-- a program to download RNA-seq data(.tsv), DNase-seq data(.narrowPeak and .bam) from the Encode project (https://www.encodeproject.org)"
    echo "OPTION:"
    echo "    -c  CELLID: pre-defined cell ID (from 1 to 60)"
    echo "    -h  show this help text"
    echo "    -r  download RNA-seq data (.tsv)"
    echo "    -p  download chromatin accessible peaks from DNase-seq data (.narrowPeak)"
    echo "    -b  download chromatin accessible readscount from DNase-seq data (.bam)"

}



#download RNA-seq data (.tsv file)
function downloadRseq()
{

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
}


#download DNase-seq data (.narrowPeak file)
function downloadDseq()
{
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
}


#download DNase-seq data (.bam file)
function downloadDseqbam()
{
zcat -f $mpath		|
tail -n +2			|
grep -P "$gname"	|
grep DNase-seq |
grep bam |
awk -F '\t' '$3=="alignments"' |
#awk -F '\t' '$30=="1" || $30=="2"' |
#awk -F '\t' '$30=="1"' |
awk -F '\t' '$46=="released"' |
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
	fbed=$Dseq/$exp/$facc.bam
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
}

while getopts "c:hrpb" Option;do
    case $Option in 
        h) printHelp
        exit 1
        ;;
        c) 
            CELLID=$OPTARG
            CELL=$(sed -n ${CELLID}p $DPATH/cells.txt|cut -f 2)
            mkdir -p "$DPATH/all/$CELLID"
            Rseq="$DPATH/all/$CELLID/rseq"
            Dseq="$DPATH/all/$CELLID/dseq"
        ;;
        r) downloadRseq
        ;;
        p) downloadDseq
        ;;
        b) downloadDseqbam
        ;;
        \?) printHelp 
        exit 1
        ;;
    esac
done


