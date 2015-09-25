#genetable.py

import sys
import string
import numpy as np
import collections as col

class Vars(object):
	def __init__(self,fields):
		self.chr=fields[0]
		self.start=fields[1]
		self.end=fields[2]
		self.autoreg=fields[6]
		self.length=float(fields[10])
		self.gene=fields[11]
		self.maf=float(fields[12])
		self.type=fields[14]

f1=open(sys.argv[1],"r")

genes=set()
auto=set()
sct={}
nct={}
ct={}
maf=col.defaultdict(list)
length={}
glen={}

for line in f1:
	v=Vars(line.rstrip().split("\t"))
	genes.add(v.gene)
	auto.add(v.autoreg)
	if v.type=="ds":
		sct[v.gene]=sct.get(v.gene,float(0))+float(1)
		ct[v.gene]=ct.get(v.gene,float(0))+float(1)
	if v.type=="dn":
		nct[v.gene]=nct.get(v.gene,float(0))+float(1)
		ct[v.gene]=ct.get(v.gene,float(0))+float(1)
	maf[v.gene].append({v.chr+v.start+v.end+v.type:v.maf})
	length[(v.gene,v.autoreg)]=v.length

f1.close()

for i in genes:
	for j in auto:
		try:
			glen[i]+=float(length[(i,j)])
		except KeyError:
		 	glen[i]=glen.get(i,float(0))
	try:
		if sct[i]==0:
			s=1
		else:
			s=sct[i]
	except KeyError:
		sct[i]=0; s=float(1)
	try:
		assert nct[i]
	except KeyError:
		nct[i]=0
	if glen[i]==0:
		g=float(1)
	else:
		g=glen[i]
	print i,int(glen[i]),int(nct[i]),int(sct[i]),int(nct[i]+sct[i]),nct[i]/s,(nct[i]+sct[i])/g,np.median([k.values()[0] for k in maf[i]])