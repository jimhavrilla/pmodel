import re
import string
import sys
import fileinput

class Record2(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.info = fields[9]
		self.interval = fields[3]
		self.database = fields[4]
		self.seqtype = fields[5]
		self.field7 = fields[6]
		self.field9 = fields[8]
		self.strand = fields[7]
		self.transid = re.search("transcript_id (.*?);",fields[10]).group(1).strip("\"")

for line in fileinput.input():
	fields=line.rstrip().split("\t")
	try:
		old_r=r
	except NameError:
		pass
	r=Record2(fields)
	try:
		if re.search("transcript_id (.*?);",r.info).group(1).strip("\"") != r.transid:
			continue
	except NameError:
		pass
	try:
		if re.search("pfamA_auto_reg.*?;",r.info).group(0) == re.search("pfamA_auto_reg.*?;",old_r.info).group(0) and re.search("gene_name.*?;",r.info).group(0) == re.search("gene_name.*?;",old_r.info).group(0):
			continue
	except NameError:
		pass
	sys.stdout.write("\t".join([r.chr,r.start,r.end,str(int(r.end)-int(r.start)+1),r.database,r.seqtype,r.field7,r.strand,r.field9,r.info])+"\n")