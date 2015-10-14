import string
import sys

class DomInts(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = int(fields[1])
		self.end = int(fields[2])
		self.ref = fields[3]
		self.alt = fields[4]
		self.info1 = fields[5:25]

class NodomInts(object):
	def __init__(self,fields):
		self.chr = fields[25]
		self.start = int(fields[26])
		self.end = int(fields[27])
		self.ref = fields[28]
		self.alt = fields[29]
		self.info2 = fields[30:50]

for line in sys.stdin:
	fields=line.rstrip().split("\t")
	wa=DomInts(fields)
	wb=NodomInts(fields)
	if not ((wb.start-len(wa.ref)+1)<wa.end or (wa.start-len(wb.ref)+1)<wb.end):
		print line
