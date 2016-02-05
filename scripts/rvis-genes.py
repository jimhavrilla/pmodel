"""
run as: python rvis-genes.py Homo_sapiens.GRCh38.76.gtf.gz
/scratch/ucgd/lustre/u1021864/serial/rvis.txt | sort -k1,1V -k2,2n > rvis.bed

to get a bed file with the positions. and rvis scores.
"""

import gzip
import sys
from collections import defaultdict
gtf = sys.argv[1]

exons = defaultdict(list)
for line in gzip.open(gtf):
    if line[0] == "#": continue
    toks = line.rstrip().split("\t")
    if toks[2] != "gene": continue

    info = toks[-1]
    gene = info.split('gene_name "')[1].split('"')[0]
    chrom, start, end = toks[0], toks[3], toks[4]

    exons[gene].append([chrom, int(start), int(end), gene])

rvis = sys.argv[2]

for i, line in enumerate(open(rvis)):
    if i == 0: continue
    if not line.strip(): continue
    gene, score, pct = line.rstrip().split("\t")

    if not gene in exons:
        continue
    for li in exons[gene]:
        li.extend([score, pct])


print "#chrom\tstart\tend\tgene\tscore\tpct"
out = []
for gene, items in exons.iteritems():
    for li in items:
        if len(li) < 5: continue
        out.append(li)
out.sort()

for li in out:
    print "\t".join(map(str, li))
