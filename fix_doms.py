import sys

cols={ 0: 'chr',
       1: 'start',
       2: 'end',
       3: 'len',
       7: 'strand',
       10: 'pfamA_id',
       12: 'gene_name',
       14: 'transcript_id',
       16: 'protein_id',
       22: 'pfamA_acc',
       24: 'pfamA_auto_reg',
       26: 'gene_id',
       30: 'exon_number',
       34: 'gene_source',
       36: 'gene_biotype',
       38: 'transcript_name',
       40: 'transcript_source',
       42: 'exon_id',
       44: 'uniq_id' ,
       45: 'num_bases_at_5x' ,
       46: 'auto_reg_length' ,
       47: 'ratio_of_cover_to_total' }

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
