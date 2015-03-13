sqlite3 $1 "create table variants(chr text,start integer,end integer,ref text,alt text,pfamA_id text,uniqid text,covct text,length_of_region text,covratio text,gene_symbol text,maf float,impact text,codons text,amino_acids text,gene_id_csq text,gene_symbol_csq text,transcript_id_csq text,exon_number_csq text,polyphen text,sift text,protein_position text,biotype text);"
sqlite3 -separator $'\t' $1 ".import $2 variants"
sqlite3 $1 <<"EOF"
create index mafidx on variants(maf);
create index geneidx on variants(gene_symbol);
create index startidx on variants(start);
create index endidx on variants(end);
EOF
