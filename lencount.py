import string
import sys
import base64

class Record1(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.domain = fields[3]
		self.autoreg = fields[4]
		self.geneid = fields[5]
		self.gene = fields[6]
		self.exonid = fields[7]
		self.uniqid = fields[8]
		self.len = fields[9]
		self.covratio = fields[10]

row=[]
length={}
for line in sys.stdin:
	fields=line.rstrip().split("\t")
	r_=Record1(fields)
	row.append([r_.chr,r_.start,r_.end,r_.domain,r_.gene,r_.autoreg,r_.uniqid,r_.covratio])
	key=base64.encodestring(r_.gene+r_.autoreg)
	try:
		length[key]=int(length[key])+int(r_.len)
	except KeyError:
		length[key]=int(r_.len)

for y in row:
	key=base64.encodestring(y[4]+y[5])
	print "\t".join(y)+"\t"+str(length[key])