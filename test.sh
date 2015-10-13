#!/bin/bash
set -e

check(){
	obs=$1
	exp=$2

if [[ "$exp" != "$obs" ]]; then
	echo "FAIL"
else
	echo "    OK"
fi
}

usage()
{
    cat << EOF

usage: $0 OPTIONS

OPTIONS can be:
    -d      Data directory

EOF
}

# Show usage when there are no arguments.
if test -z "$1"
then
    usage
    exit
fi

DATA=

# Check options passed in.
while getopts "h d:" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        d)
            DATA=$OPTARG
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

if [ -z "$DATA" ]
then
    echo "Data directory not set"
    usage
    exit
fi

if [ ! -d "$DATA" ]
then
    echo "Data directory does not exist"
    usage
    exit
fi

##############################################################################
echo "pmodel.t0... test that the number of transcripts in the exons file matches the number of transcripts after utr removal" 

if [ ! -f "$DATA/exons.bed" ]
then
    echo "COULD NOT FIND $DATA/exons.bed"
    exit
fi

if [ ! -f "$DATA/exons_sans_utrs.bed" ]
then
    echo "COULD NOT FIND $DATA/exons_sans_utrs.bed"
    exit
fi

exp=`cat $DATA/exons.bed \
    | grep -v "^#" \
    | grep transcript_id \
    | cut -f 7 \
    | cut -d ";" -f2 \
    | uniq \
    | sort -u \
    | wc -l`

obs=`cat $DATA/exons_sans_utrs.bed \
    | grep -v "^#" \
    | grep transcript_id \
    | cut -f 7 \
    | cut -d ";" -f2 \
    | uniq \
    | sort -u \
    | wc -l`
check $obs $exp


##############################################################################
echo "pmodel.t1... test that number of uniq genes before filtering on maximum length is equal to the number in transcripts.txt"

if [ ! -f "$DATA/appris_data.principal.txt" ]
then
    echo "COULD NOT FIND $DATA/appris_data.principal.txt"
    exit
fi

if [ ! -f "$DATA/transcripts.txt" ]
then
    echo "COULD NOT FIND $DATA/transcripts.txt"
    exit
fi

exp=$(cut -f 2 $DATA/appris_data.principal.txt | sort -u | wc -l | awk '{print $1;}')
obs=$(wc -l < $DATA/transcripts.txt)
check $obs $exp

##############################################################################
echo "pmodel.t2... test that the nodoms don't intersect the doms"

if [ ! -f "$DATA/doms-any-coverage.bed" ]
then
    echo "COULD NOT FIND $DATA/doms-any-coverage.bed"
    exit
fi

if [ ! -f "$DATA/nodoms-any-coverage.bed" ]
then
    echo "COULD NOT FIND $DATA/nodoms-any-coverage.bed"
    exit
fi

obs=$(bedtools intersect -a $DATA/doms-any-coverage.bed -b $DATA/nodoms-any-coverage.bed -wo \
	| python test-different-transcripts.py)
exp=""
check $obs $exp


##############################################################################
echo "pmodel.t3... test that exons-(doms+nodoms) gives empty results"

if [ ! -f "$DATA/doms-any-coverage.bed" ]
then
    echo "COULD NOT FIND $DATA/doms-any-coverage.bed"
    exit
fi

if [ ! -f "$DATA/nodoms-any-coverage.bed" ]
then
    echo "COULD NOT FIND $DATA/nodoms-any-coverage.bed"
    exit
fi

if [ ! -f "$DATA/exons_sans_utrs.bed" ]
then
    echo "COULD NOT FIND $DATA/exons_sans_utrs.bed"
    exit
fi

if [ ! -f "$DATA/transcripts.txt" ]
then
    echo "COULD NOT FIND $DATA/transcripts.txt"
    exit
fi

obs=$(bedtools subtract -a <(cat $DATA/doms-any-coverage.bed $DATA/nodoms-any-coverage.bed|cut -f 1-22) -b $DATA/exons_sans_utrs.bed)
exp=""
check $obs $exp

##############################################################################

echo "pmodel.t4... test that (doms+nodoms)-exons gives empty results"

if [ ! -f "$DATA/doms-any-coverage.bed" ]
then
    echo "COULD NOT FIND $DATA/doms-any-coverage.bed"
    exit
fi

if [ ! -f "$DATA/nodoms-any-coverage.bed" ]
then
    echo "COULD NOT FIND $DATA/doms-any-coverage.bed"
    exit
fi

if [ ! -f "$DATA/exons_sans_utrs.bed" ]
then
    echo "COULD NOT FIND $DATA/exons_sans_utrs.bed"
    exit
fi
cut -f 2 $DATA/transcripts.txt > $DATA/foo8738
obs=$(bedtools subtract -b <(cat $DATA/doms-any-coverage.bed $DATA/nodoms-any-coverage.bed | cut -f 1-22) -a <(grep -wF -f $DATA/foo8738 $DATA/exons_sans_utrs.bed))
exp=""
check $obs $exp

###############################################################################

echo "pmodel.t5... test that doms and nodoms after coverage filtering are still exclusive"

if [ ! -f $DATA/doms.bed ]
then
	echo "COULD NOT FIND $DATA/doms.bed"
	exit
fi

if [ ! -f $DATA/nodoms.bed ]
then
	echo "COULD NOT FIND $DATA/nodoms.bed"
	exit
fi

obs=$(bedtools intersect -a $DATA/doms.bed -b $DATA/nodoms.bed -wo | python test-different-transcripts.py)
exp=""
check $obs $exp

################################################################################

echo "pmodel.t6... test that dom variants and nodom variants are mutually exclusive"

if [ ! -f "$DATA/VEPEXAC3filter.vcf.gz" ]
then
    echo "COULD NOT FIND $DATA/VEPEXAC3filter.vcf.gz"
    exit
fi

if [ ! -f "$DATA/nodoms.bed" ]
then
    echo "COULD NOT FIND $DATA/nodoms.bed"
    exit
fi

if [ ! -f "$DATA/doms.bed" ]
then
    echo "COULD NOT FIND $DATA/doms.bed"
    exit
fi

obs=$(bedtools intersect -a $DATA/nodomintfilter.bed -b $DATA/domintfilter.bed)
exp=""
check $obs $exp
