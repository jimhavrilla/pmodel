#lollipop.sh
DATA=$1; OUT=$2; GENE=$3; MOD=$4; MAF1=$5; MAF2=$6
#sort -k6,6 -k8,12 $DATA/allintfilter.bed | sed '1d' | sed -e '/[[:space:]]na[[:space:]]/ {d;}' > $DATA/lollipop.bed
grep $GENE $DATA/lollipop.bed > $OUT/$GENE'fooint'
if [[ $MOD -eq g && $MAF1 -eq 0 ]]
then
	OF=$OUT/$GENE
else 
	OF=$OUT/$GENE.$MOD.$MAF1.$MAF2
fi
expr=$(perl -pe 's/(\S*\s){11}(\S*)\s(\S*\s){3}(\S*)\s(\S*\s){6}(\S*).*\n/$4 $6 $2 /g' $OUT/$GENE'fooint' | perl -pe 's/(\w*|\*|\.)\/?(\S*) (\?|\d*-?\d*)\/\d* (\d*\.\d*e?-?\d*)/$1$3$2 $4/g' | paste <(printf $GENE) - )
echo $expr > foo
lollipops -m=1 -o=$OF.svg $expr
convert $OF.svg $OF.png # uses ImageMagick convert
