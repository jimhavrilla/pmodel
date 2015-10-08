#!/bin/bash

DATA=$1
SOFTWARE=$2

function usage () {
	cat << EOF
	
usage: $0 $DATA, $SOFTWARE opts
	where opts are:
		-h	Display this help message
		-g	get files - pfam defs, exons, vcf - and VEP annotate:
				0 means do all
				1 means don't redownload or vep annotate
				2 means redownload, but don't vep annotate
		-c	generate domain and nodom files
				value input is the coverage mean chosen,
				default	is 5, for both domain and nodom
		-i	run intersection:
		-t	make tables:
				value input is the maf cutoff chosen,
				this only modifies table values
		-l	make lollipop:
				value input is out folder, gene,
				relational operator (g,l,r), maf1,
				maf2 (optional)
EOF
}

if test -z "$1"
then
	usage
	exit
fi

while getopts "hg:cc:it:l:" OPTION
do
	case $OPTION in
	h)
		usage
		exit
		;;
	g)
		bash getfiles.sh $OPTARG
		;;
	c)
		if [ !$OPTARG ];
		then
			OPTARG=5
		fi
	#	bash create.sh $OPTARG
		;;
	i)
		bash intersect.sh
		;;
	t)
		bash tables.sh $OPTARG
		;;
	l)
		array=(${OPTARG//,/ })
		if [ ${#array[@]} == 5 ];
		then
			bash lollipop.sh $DATA ${array[0]} ${array[1]} ${array[2]} ${array[3]} ${array[4]}
		elif [ ${#array[@]} == 4 ];
		then
			bash lollipop.sh $DATA ${array[0]} ${array[1]} ${array[2]} ${array[3]}
		else
			usage
			exit 1
		fi
		;;
	esac
done
