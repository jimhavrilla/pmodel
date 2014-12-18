if [[ $2 == *'%'* ]]
	then sqlite3 $1 "select * from variants where gene_symbol LIKE '$2' and maf>'$3';"
else
	sqlite3 $1 "select * from variants where gene_symbol='$2' and maf>'$3';"
fi