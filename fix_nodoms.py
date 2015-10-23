import sys

cols={ 0: 'chr',
       1: 'start',
       2: 'end',
       3: 'len',
       5: 'gene_id',
       7: 'transcript_id',
       9: 'exon_number',
       11: 'gene_name',
       13: 'gene_source',
       15: 'gene_biotype',
       17: 'transcript_name',
       19: 'transcript_source',
       21: 'exon_id',
       23: 'uniq_id',
       24: 'num_bases_at_5x' ,
       25: 'auto_reg_length' ,
       26: 'ratio_of_cover_to_total' }

O = []
for col in cols:
    O.append(cols[col])

print '#' + '\t'.join(O)

for l in sys.stdin:
    A = l.rstrip().split('\t')
    O = []
    for col in cols:
        O.append(A[col])
    print '\t'.join(O)
