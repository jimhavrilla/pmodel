if [[ $2 == '%' ]]
	then
	if [[ $3 == g ]] # greater than
		then sqlite3 $1 "select * from variants where maf>'$4';"
	fi
	if [[ $3 == l ]] #less than
		then sqlite3 $1 "select * from variants where maf<'$4';"
	fi
	if [[ $3 == r ]] # range
		then sqlite3 $1 "select * from variants where maf>'$4' and maf <'$5';"
	fi
fi
if [[ $2 != '%' ]]
then
	if [[ $3 == g ]] # greater than
		then sqlite3 $1 "select * from variants where gene_symbol='$2' and maf>'$4';"
	fi
	if [[ $3 == l ]] #less than
		then sqlite3 $1 "select * from variants where gene_symbol='$2' and maf<'$4';"
	fi
	if [[ $3 == r ]] # range
		then sqlite3 $1 "select * from variants where gene_symbol='$2' and maf>'$4' and maf <'$5';"
	fi
fi