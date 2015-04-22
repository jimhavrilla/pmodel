import string
import sys
import base64

class Record1(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.domain = fields[3]
		self.gene = fields[4]
		self.autoreg = fields[5]
		self.geneid = fields[6]
		self.exonid = fields[7]
		self.uniqid = fields[8] 
		self.covct = fields[9]
		self.len = fields[10]
		self.covratio = fields[11]

row=[]
length={}
for line in sys.stdin:
	fields=line.rstrip().split("\t")
	r_=Record1(fields)
	row.append([r_.chr,r_.start,r_.end,r_.domain,r_.gene,r_.autoreg,r_.uniqid,r_.covct,r_.covratio])
	key=base64.encodestring(r_.gene+r_.autoreg)
	try:
		length[key]=int(length[key])+int(r_.len)
	except KeyError:
		length[key]=int(r_.len)

for y in row:
	key=base64.encodestring(y[4]+y[5])
	print "\t".join(y)+"\t"+str(length[key])