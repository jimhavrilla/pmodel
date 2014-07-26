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
		self.strand = fields[3]
		self.domain = fields[4].strip("\"")
		self.gene = fields[5].strip("\"")
		self.maf = float(fields[6])
		self.type = fields[7]

ct=0
nct=0
sct=0
mmaf=[]
mmafct=0
for line in fileinput.input():
	fields=line.rstrip().split("\t")
	try:
		old_r=r
	except NameError:
		pass
	r=Record3(fields)
	try:
		if old_r.domain != r.domain:
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
	if r.type=="coding" or r.type=="coding-notMod3" or r.type=="coding-synonymous" or r.type=="coding-synonymous-near-splice" or r.type=="codingComplex":
		sct=sct+1
		ct=ct+1
		mmaf.append(r.maf)
	if r.type=="missense" or r.type=="missense-near-splice" or r.type=="stop-gained" or r.type=="stop-gained-near-splice" or r.type=="stop-lost":
		nct=nct+1
		ct=ct+1
		mmaf.append(r.maf)