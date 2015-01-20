if [[ $3 == g ]] # greater than
	then
	if [[ $2 == *'%'* ]]
		then sqlite3 $1 "select * from variants where gene_symbol LIKE '$2' and maf>'$4';"
	else
		sqlite3 $1 "select * from variants where gene_symbol='$2' and maf>'$4';"
	fi
fi
if [[ $3 == l ]] #less than
	then
	if [[ $2 == *'%'* ]]
		then sqlite3 $1 "select * from variants where gene_symbol LIKE '$2' and maf<'$4';"
	else
		sqlite3 $1 "select * from variants where gene_symbol='$2' and maf<'$4';"
	fi
fi
if [[ $3 == r ]] # range
	then
	if [[ $2 == *'%'* ]]
		then sqlite3 $1 "select * from variants where gene_symbol LIKE '$2' and maf>'$4' and maf <'$5';"
	else
		sqlite3 $1 "select * from variants where gene_symbol='$2' and maf>'$4' and maf <'$5';"
	fi
fi