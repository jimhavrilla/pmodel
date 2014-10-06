# maketable.py - to make a table for counts vs. dn/ds plot and bp vs. var/bp plot
import string
import sys
import fileinput
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
		self.eamaf = float(fields[8])
		self.aamaf = float(fields[9])
		self.maf = float(fields[10])
		self.impact = fields[11]
		self.type = fields[12]
		self.info = fields[13:23]

def count(r_):



ct=0
nct=0
sct=0
mmaf=[]
mmafct=0
for line in fileinput.input():
	fields=line.rstrip().split("\t")
	try:
		old_r=r_
	except NameError:
		pass
	r_=Record3(fields)
	try:
		if old_r.domain != r_.domain:
			with warnings.catch_warnings(): #catches NaN warnings and empty slices
				warnings.simplefilter("ignore")
				mmaf=numpy.median(mmaf)
			if ct==0:
				mmaf=[]
				continue
			try:
				sys.stdout.write("\t".join([old_r.domain,old_r.gene,str(nct),str(sct),str(ct),str(round(float(nct)/float(sct),2)),str(mmaf)])+"\n")
			except ZeroDivisionError:
				sys.stdout.write("\t".join([old_r.domain,old_r.gene,str(nct),str(sct),str(ct),str(round(float(nct)/float(sct+1),2)),str(mmaf)])+"\n")
			ct=0
			nct=0
			sct=0
			mmaf=[]
			continue
	except NameError:
		pass
	if r_.type=="ds":
		sct=sct+1
		ct=ct+1
		mmaf.append(r_.maf)
	if r_.type=="dn":
		nct=nct+1
		ct=ct+1
		mmaf.append(r_.maf)
