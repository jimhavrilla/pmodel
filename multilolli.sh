DATA=$1; INFILE=$2; OUT=$3; MOD=$4; MAF1=$5; MAF2=$6;
mkdir $OUT/$INFILE
parallel bash lollipop.sh $DATA $OUT/$INFILE {} $MOD $MAF1 $MAF2 ::: $(cat $INFILE)