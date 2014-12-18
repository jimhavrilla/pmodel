import string
import sys
from numpy import std

class Record1(object):
	def __init__(self, fields):
		self.domain = fields[0]
		self.gene = fields[1]
		self.nct = fields[2]
		self.sct = fields[3]
		self.ct = fields[4]
		self.domcount = fields[5]
		self.bp = fields[6]
		self.mmaf = fields[7]
		self.dnds = fields[8]
		self.vbpratio = fields[9]

gene={}

#domain list, count of domains, bps, bp sum, domcounts, domcount sum, dnds, min/max/stdev dnds, rank metric
for line in sys.stdin:
	fields=line.rstrip().split("\t")
	r_=Record1(fields)
	for x in r_.gene.split(","):
		try:
			foo=gene[x]
			foo[0]=x
			foo[1].add(r_.domain)
			foo[2]=foo[2]+1
			foo[3].append(int(r_.bp))
			foo[4]=sum(foo[3])
			foo[5].append(int(r_.domcount))
			foo[6]=sum(foo[5])
			foo[7].append(float(r_.dnds))
			foo[8]=min(foo[7])
			foo[9]=max(foo[7])
			foo[10]=std(foo[7])
			foo[11]=(foo[9]-foo[8])*foo[10]*foo[4]*foo[6]/100 #(max-min)*stdev*bp*domcount/100
			gene[x]=foo
		except KeyError:
			gene[x]=[x,set([r_.domain]),1,[int(r_.bp)],0,[int(r_.domcount)],0,[float(r_.dnds)],0,0,0,0]

for x in gene:
	sys.stdout.write('\t'.join([gene[x][0],','.join(gene[x][1]),str(gene[x][2]),str(gene[x][4]),str(gene[x][6]),','.join([str(y) for y in gene[x][7]]),str(gene[x][8]),str(gene[x][9]),str(gene[x][10]),str(gene[x][11])])+'\n')
