#lollipop.sh
DATA=$1; OUT=$2; GENE=$3; MOD=$4; MAF1=$5; MAF2=$6
bash query.sh $DATA/var.db $GENE $MOD $MAF1 $MAF2 | python formatvar.py | sort -k6,6 -k8,8 | grep -v -w "na" > $DATA/$GENE'fooint'
lollipops -m=1 -o=$OUT/$GENE.svg $(perl -pe 's/(\S*\s){8}(\S*)\s(\S*\s){3}(\S*)\s(\S*\s){6}(\S*).*\n/$4,$6,$2 /g' $DATA/$GENE'fooint' | perl -pe 's/(\w*|\*|\.)\/?(\S*),(\d*-?\d*)\/\d*,(\d*\.\d*e?-?\d*)/$1$3$2 $4/g' | paste <(printf $GENE) - )
convert $OUT/$GENE.svg $OUT/$GENE.png # uses ImageMagick convert
#rm $DATA/$GENE'fooint'