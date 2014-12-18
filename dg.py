# dg.py - 
from __future__ import print_function
import string
import sys
import numpy
import warnings

class Record3(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.ref = fields[3]
		self.alt = fields[4]
		self.domain = fields[5].rstrip(";").strip("\"")
		self.uniqid = fields[6]
		self.gene = fields[7]
		self.maf = float(fields[8])
		self.impact = fields[9]
		self.type = fields[10]
		self.info = fields[11:21]

def ratiocalc(nct,sct):
	try:
		ns=round(float(nct)/float(sct),2)
	except ZeroDivisionError:
		ns=round(float(nct)/float(sct+1),2)

	return ns

def count(old_r,maf,ct,nct,sct):
	with warnings.catch_warnings(): #catches NaN warnings and empty slices
		warnings.simplefilter("ignore")
		mmaf=numpy.median(maf)
	table.append([old_r.domain,old_r.gene,str(ct),str(nct),str(sct),ratiocalc(nct,sct),str(mmaf)])
	return table

ct=0
nct=0
sct=0
maf=[]
table=[]
old_r=None
f1=open(sys.argv[1],"r")
for line in f1:
	fields=line.rstrip().split("\t")
	r_=Record3(fields)
	if old_r!=None and old_r.gene!=r_.gene: # assumes data is sorted by domain then gene, which it should be
		count(old_r,maf,ct,nct,sct)
		nct=0;sct=0;ct=0;
		maf=[]
	if r_.type=="ds":
		sct=sct+1
		ct=ct+1
		maf.append(r_.maf)
	if r_.type=="dn":
		nct=nct+1
		ct=ct+1
		maf.append(r_.maf)
	
	old_r=r_

count(r_,maf,ct,nct,sct)

f2=open(sys.argv[2],"w")
xlist=[]
old_x=None #x[0] domain, x[1] gene, x[2] ct, x[3] nct, x[4] sct, x[5] dn/ds, x[6] mmaf
for x in table:
	if old_x!=None and old_x[0] != x[0]:
		l=[y[5] for y in xlist[0:len(xlist)]]
		m=numpy.mean(l)
		s=numpy.std(l)
		for y in xlist:
			if s==0.0:
				s=1.0
			print("\t".join([y[0],y[1],y[2],y[3],y[4],str(y[5]),str((y[5]-m)/s),y[6]])+"\n",file=f2)
		xlist=[]
	xlist.append(x)
	old_x=x

l=[y[5] for y in xlist[0:len(xlist)]]
m=numpy.mean(l)
s=numpy.std(l)
for y in xlist:
	if s==0.0:
		s=l[0]
		print("\t".join([y[0],y[1],y[2],y[3],y[4],str(y[5]),str((y[5]-m)/s),y[6]])+"\n",file=f2)