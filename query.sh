if [[ $3 == g ]] # greater than
	then sqlite3 $1 "select * from variants where gene_symbol LIKE '$2' and maf>'$4';"
fi
if [[ $3 == l ]] #less than
	then sqlite3 $1 "select * from variants where gene_symbol LIKE '$2' and maf<'$4';"
fi
if [[ $3 == r ]] # range
	then sqlite3 $1 "select * from variants where gene_symbol LIKE '$2' and maf>'$4' and maf <'$5';"
fi