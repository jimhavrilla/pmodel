import string
import sys

class Record(object):
	def __init__(self, fields):
		self.gene = fields[5]
		self.phenotype = fields[11]

f1=open(sys.argv[1],"r");f2=open(sys.argv[2],"w");
for line in f1:
	fields=line.rstrip().split("|")
	r_=Record(fields)
	for x in r_.gene.split(","):
		x=x.strip(" ")
		count=str(len(r_.phenotype.split(";")))
		if r_.phenotype == '':
			r_.phenotype="none"
			count=str(0)
		f2.write("\t".join([x,r_.phenotype,count])+"\n")