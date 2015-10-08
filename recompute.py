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
		self.uniqid = fields[7]
		self.covratio = fields[8]
		self.length = fields[9]
		self.dn = fields[10]
		self.ds = fields[11]
		self.na = fields[12]
		self.dnds = fields[13]
		self.density = fields[14] #ignored for recomputing
		self.prevalence = fields[15]
		self.mafs = fields[16]
		self.impacts = fields[17]
		self.types = fields[18]
		self.starts = fields[19]
		self.ends = fields[20]

true_mafs=[] #includes alternate alleles for density count
rare_mafs=[]
for line in sys.stdin:
	r_=Record(line.rstrip().split("\t"))
	mafs=r_.mafs.split(",")
	for x in mafs:
		if x==".": #takes care of non-intersecting regions
			mafs=[]
			continue
		for y in x.split('|'):
			true_mafs.append(y)
			if (float(y) <=0.001):
				rare_mafs.append(y)
	density=str(float(len(true_mafs))/float(r_.length))
	try:
		fvrv=str(float(len(rare_mafs))/float(len(true_mafs)))
	except ZeroDivisionError:
		fvrv="0.0"
	sys.stdout.write("\t".join([r_.chr,r_.start,r_.end,r_.transid,r_.domain,r_.gene,r_.autoreg,r_.covratio,r_.length,r_.dn,r_.ds,r_.na,r_.dnds,density,fvrv,r_.prevalence,r_.mafs,r_.impacts,r_.types,r_.starts,r_.ends,'0'])+"\n")
	rare_mafs=[]
	true_mafs=[]