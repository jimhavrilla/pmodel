import string
import fileinput
import sys
import re

class Record1(object):
	def __init__(self,fields):
		self.domain=fields[0]
		self.gene=fields[1]
		self.nct=fields[2]
		self.sct=fields[3]
		self.ct=fields[4]
		self.dct=fields[5]
		self.bp=fields[6]
		self.maf=fields[7]
		self.dnds=fields[8]
		self.vbp=fields[9]

class Record2(object):
	def __init__(self,fields):
		self.pfamacc=fields[0]
		self.domain=fields[1]
		self.clanacc=fields[2]
		self.clan=fields[3]
		self.dct=fields[4]

DATA="/Users/jmh2tt/work/data/pmodeldata"
f1=open(DATA+"/table.txt","r")#sys.argv[1]
f2=open(DATA+"/count_human_pfam_clan.tab","r")#sys.argv[2]
for line1 in f1:
	r1=Record1(line1.rstrip().split("\t"))
	for line2 in f2:
		r2=Record2(line2.rstrip().split("\t"))
		if r1.domain==r2.domain:
			print "\t".join([r1.domain,r1.gene,r1.nct,r1.sct,r1.ct,r1.dct,r1.bp,r1.maf,r1.dnds,r1.vbp,r2.pfamacc,r2.domain,r2.clanacc,r2.clan,r2.dct])
	f2.seek(0)

f1.close()
f2.close()
