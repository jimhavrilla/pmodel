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
		self.ccdsid=fields[3]
		self.autoreg=fields[6]
		self.length=fields[10]
		self.gene=fields[11]
		self.maf=fields[12]
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

#dict.update if tuple or list, dict.add if single value

for line in f1:
	v=Vars(line.rstrip().split("\t"))
	genes.add(v.gene)
	auto.add(v.autoreg)
	if v.type="ds":
		sct[v.gene]=foo.get(sct,0)+1
		ct[v.gene]=foo.get(ct,0)+1
	if v.type="dn":
		nct[v.gene]=foo.get(nct,0)+1
		ct[v.gene]=foo.get(ct,0)+1
	maf[v.gene].append({v.chr+v.start+v.end+v.type:v.maf})
	length[v.gene,v.autoreg]=v.length

f1.close()

for i in genes:
	for j in auto:
		glen[i]+=length[i,j]
	if sct[i]==0:
		s=1
	else:
		s=sct[i]
	if glen[i]==0:
		g=1
	else:
		g=glen[i]
	print i,glen[i],nct[i],sct[i],nct[i]+sct[i],nct[i]/s,(nct[i]+sct[i])/glen[i],np.median([k.values()[0] for k in maf[i]])