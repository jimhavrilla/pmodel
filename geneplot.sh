#geneplot.sh
GENE=$1; MAF=$2; DATA=$3
bash query.sh $DATA/lolli.db $GENE $MAF | python formatvar.py | sort -k6,6 -k8,8 | grep -v -w "na" > $DATA/fooint
lollipops $(perl -pe 's/(\S*\s){12}(\S*)\s(\S*\s){6}(\S*).*\n/$2,$4 /g' $DATA/fooint | perl -pe 's/(\w)\/?(\S)?,(\d*?)\/\d*/$1$3$2/g' | paste <(printf $GENE) - )
convert $GENE.svg $GENE.png # uses ImageMagick convert