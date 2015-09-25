import string
import sys

class Record(object):
	def __init__(self,fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.transid = fields[3]
		self.domain = fields[4]
		self.gene = fields[5]
		self.autoreg = fields[6]
		self.length = fields[7]
		self.dn = fields[8]
		self.ds = fields[9]
		self.dnds = fields[10]
		self.covratio = fields[11]
		self.mafs = fields[12]
		self.impacts = fields[13]
		self.types = fields[14]
		self.starts = fields[15]
		self.ends = fields[16]

rare_mafs=[]
for line in sys.stdin:
	r_=Record(line.rstrip().split("\t"))
	mafs=r_.mafs.split(",")
	for x in mafs:
		if x==".":
			mafs=[]
			continue
		if (float(x) <=0.001):
			rare_mafs.append(x)
	density=str(float(len(mafs))/float(r_.length))
	try:
		fvrv=str(float(len(rare_mafs))/float(len(mafs)))
	except ZeroDivisionError:
		fvrv="0.0"
	sys.stdout.write("\t".join([r_.chr,r_.start,r_.end,r_.transid,r_.domain,r_.gene,r_.autoreg,r_.length,r_.dn,r_.ds,r_.dnds,r_.covratio,density,fvrv,r_.mafs,r_.impacts,r_.types,r_.starts,r_.ends,'0'])+"\n")
	rare_mafs=[]