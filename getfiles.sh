#!/bin/bash
set -e
#getfiles.sh - gets files and formats the initial files used to generate all
#other files

# bill's domain count file:

## BEGIN SCRIPT
usage()
{
    cat << EOF

usage: $0 OPTIONS

OPTIONS can be:
    -d      Data directory
    -s      Software directory

EOF
}

# Show usage when there are no arguments.
if test -z "$1"
then
    usage
    exit
fi

DATA=
SOFTWARE=

# Check options passed in.
while getopts "h d:s:" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        d)
            DATA=$OPTARG
            ;;
        s)
            SOFTWARE=$OPTARG
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

if [ -z "$SOFTWARE" ]
then
    echo "Software directory not set"
    usage
    exit
fi
if [ ! -d "$SOFTWARE" ]
then
    echo "Software directory does not exist"
    usage
    exit
fi

if [ ! -f "$DATA/human_pfam.counts" ]
then
    echo "$DATA/human_pfam.counts not found and cannot be recreated."
    exit
#    mysql \
#        -N \
#        --raw \
#        -h wrpxdb.its.virginia.edu \
#        -u web_user -pfasta_www pfam27 
#        < count_human_pfam.sql \
#            > $DATA/human_pfam.counts
fi


# for vcf:

if [ ! -f "$DATA/ExAC.r0.3.sites.vep.vcf.gz" ]
then
    wget -P $DATA ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz
fi

if [ ! -f "$DATA/VEPEXAC3.vcf.gz" ]
then
    perl $SOFTWARE/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/ExAC.r0.3.sites.vep.vcf.gz --cache --merged --sift b --polyphen b --symbol --numbers --biotype --total_length --allele_number -o $DATA/VEPEXAC3.vcf --vcf --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM --force_overwrite --port 3337 --dir ~/.vep/ --compress "gunzip -c" --offline
fi

# for exons: 

if [ ! -f "$DATA/Homo_sapiens.GRCh37.75.gtf.gz" ]
then
    echo "$DATA/Homo_sapiens.GRCh37.75.gtf.gz not found"
    exit
fi

if [ ! -f "$DATA/exons.bed" ]
then
    zcat < $DATA/Homo_sapiens.GRCh37.75.gtf.gz \
        | grep protein_coding$'\t'exon  \
        | perl -pe 's/protein_coding\texon\t//g' \
        | grep -e '^\d*[0-9X-Y]\t' \
        | perl gtf2bed.pl \
        | sort -k10,10 \
        > $DATA/exons.bed
fi

if [ ! -f "$DATA/utrs.bed" ]
then
    zcat < $DATA/Homo_sapiens.GRCh37.75.gtf.gz \
        | grep protein_coding$'\t'UTR \
        | perl -pe 's/protein_coding\tUTR\t//g' \
        | grep -e '^\d*[0-9X-Y]\t' \
        | perl gtf2bed.pl \
        | sort -k10,10 \
        > $DATA/utrs.bed
fi

python filterutrs.py $DATA/exons.bed $DATA/utrs.bed > $DATA/GRCh37.bed

# for correct transcripts:

wget -P $DATA http://apprisws.bioinfo.cnio.es/pub/current_release/ensembl_datafiles/species/homo_sapiens/GRCh37/appris_data.principal.txt

# for vcf coverage:

for chrom in {1..22} X Y ; do wget -P $DATA ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/coverage/Panel.chr$chrom.coverage.txt.gz; done
gzcat $DATA/Panel.chr*.coverage.txt.gz | awk '/^#/ {sub(/#.*/,"");getline;}1 {print $1,$2-1,$2,$3}' | tr -s " " "\t" > $DATA/coverage.bed

# bill's gtfs:

for chrom in {1..22} X Y ; do wget -P $DATA http://fastademo.bioch.virginia.edu/pfam_dna/Homo_sapiens.chr$chrom.pfam.gtf; done
for chrom in {1..22} X Y
do
sort -k2,2n $DATA/Homo_sapiens.chr$chrom.pfam.gtf | python rearrange.py | sort -k1,1 -k2,2n  > $DATA/chr$chrom.bed
done
cat $DATA/chr*.bed | sort -k1,1 -k2,2n > $DATA/all.bed
rm $DATA/chr*
