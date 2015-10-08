import string
import sys

class RecordUniq(object):
	def __init__(self, fields):
		self.domain = fields[0]
		self.gene = fields[1]
		self.autoreg = fields[2]
		self.uniqid = fields[3]
		self.covratio = fields[4]
		self.length = fields[5]
		self.dn = fields[6]
		self.ds = fields[7]
		self.na = fields[8]
		self.dnds = fields[9]
		self.density = fields[10]

class RecordNone(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.domain = fields[3]
		self.gene = fields[4]
		self.autoreg = fields[5]
		self.uniqid = fields[6]
		self.covratio = fields[7]
		self.length = fields[8]

uset=set()

with open(sys.argv[1],"r") as f1:
	for line in f1:
		uniqs=line.rstrip()
		u=RecordUniq(uniqs.split())
		uset.add(u.autoreg)

with open(sys.argv[2],"r") as f2:
	for line in f2:
		nones=line.rstrip()
		n=RecordNone(nones.split())
		if n.autoreg not in uset:
			with open(sys.argv[1],"a") as f3:
				f3.write("\t".join([n.domain,n.gene,n.autoreg,n.uniqid,n.covratio,n.length,str(0),str(0),str(0),str(0),str(0)])+"\n")