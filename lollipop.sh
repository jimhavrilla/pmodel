#lollipop.sh
DATA=$1; OUT=$2; GENE=$3; MOD=$4; MAF1=$5; MAF2=$6
#sort -k6,6 -k8,12 $DATA/allintfilter.bed | sed '1d' | sed -e '/[[:space:]]na[[:space:]]/ {d;}' > $DATA/lollipop.bed
grep $GENE $DATA/lollipop.bed > $OUT/$GENE'raw'
if [ -z "$MAF2"]
then
	MAF2=1
fi
OF=$OUT/$GENE.$MOD.$MAF1.$MAF2
awk -v maf1="$MAF1" -v maf2="$MAF2" '{if ($12 >= maf1 && $12 <=maf2) print}' $OUT/$GENE'raw' > $OUT/$GENE'filter'
expr=$(perl -pe 's/(\S*\s){11}(\S*)\s(\S*\s){3}(\S*)\s(\S*\s){6}(\S*).*\n/$4 $6 $2 /g' $OUT/$GENE'filter' | perl -pe 's/(\w*|\*|\.)\/?(\S*) (\?|\d*-?\d*)\/\d* (\d*\.\d*e?-?\d*)/$1$3$2 $4/g' | paste <(printf $GENE) - )
echo $expr > foo
lollipops -m=1 -o=$OF.svg $expr
convert $OF.svg $OF.png # uses ImageMagick convert
