# maketable.py - to make a table for counts vs. dn/ds plot and bp vs. var/bp plot
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

def count(old_r,maf,ct,tnct,tsct,geneset):
	with warnings.catch_warnings(): #catches NaN warnings and empty slices
		warnings.simplefilter("ignore")
		mmaf=numpy.median(maf)
	foo=''
	for x in geneset:
		foo=foo+x+","
	foo=foo.rstrip(",")
	sys.stdout.write("\t".join([old_r.domain,foo,str(tnct),str(tsct),str(ct),str(ratiocalc(tnct,tsct)),str(mmaf)])+"\n")
	ct=0
	tnct=0
	tsct=0
	maf=[]
	return [maf,mmaf,ct,nct,sct]

f1=open(sys.argv[1])
geneset = set([])
ct=0
nct=0
sct=0
tnct=0
tsct=0
maf=[]
dnds=[]
old_r=None
for line in f1:
	fields=line.rstrip().split("\t")
	r_=Record3(fields)
	if r_.type=="ds":
		sct=sct+1
		tsct=tsct+1
		ct=ct+1
		maf.append(r_.maf)
	if r_.type=="dn":
		nct=nct+1
		tnct=tnct+1
		ct=ct+1
		maf.append(r_.maf)
	if r_.gene not in geneset:
		dnds.append(ratiocalc(nct,sct))
		nct=0;sct=0
		geneset.add(r_.gene)
	if old_r!=None and old_r.domain != r_.domain:
		[maf,mmaf,ct,tnct,tsct]=count(old_r,maf,ct,tnct,tsct,geneset)
		geneset.clear()
	
	old_r=r_

count(r_,maf,ct,tnct,tsct,geneset)
