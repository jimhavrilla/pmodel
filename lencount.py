import string
import sys
import base64

class Record1(object):
	def __init__(self, fields):
		self.info = fields[0:46]
		self.gene = fields[12]
		self.autoreg = fields[24]
		self.len = fields[46]
		self.covratio = fields[47]

row=[]
row2=[]
keys=[]
length={}
for line in sys.stdin:
	fields=line.rstrip().split("\t")
	r_=Record1(fields)
	row.append(r_.info)
	keys.append(r_.gene+r_.autoreg)
	row2.append(r_.covratio)
	key=base64.encodestring(r_.gene+r_.autoreg)
	try:
		length[key]=int(length[key])+int(r_.len)
	except KeyError:
		length[key]=int(r_.len)

for y in range(1,len(row)):
	key=base64.encodestring(keys[y])
	print "\t".join(row[y])+"\t"+str(length[key])+"\t"+str(row2[y])