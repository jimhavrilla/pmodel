import re
import string
import sys
import fileinput

class Record(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.info = fields[8]+";"
		self.start = fields[3]
		self.end = fields[4]
		self.database = fields[1]
		self.seqtype = fields[2]
		self.field6 = fields[5]
		self.field8 = fields[7]
		self.strand = fields[6]

for line in fileinput.input():
	fields=line.rstrip().split("\t")
	try:
		old_r=r
	except NameError:
		pass
	r=Record(fields)
	if int(r.start)>int(r.end):
		tmp=r.start
		r.start=r.end
		r.end=tmp
	m=re.search("pfamA_id.*?; ",r.info)
	n=re.sub(" pfamA_id.*?;","",r.info)
	r.info=m.group(0)+n
	if("ccds_id" in r.info):
		m=re.search(" ccds_id.*?;",r.info)
		n=re.sub("ccds_id.*?; ","",r.info)
		r.info=n+m.group(0)
	# try:
	# 	if re.search("pfamA_id.*?;",r.info).group(0) == re.search("pfamA_id.*?;",old_r.info).group(0):
	# 		continue
	# except NameError:
	# 	pass
	sys.stdout.write("\t".join([r.chr,str(int(r.start)-1),r.end,str(int(r.end)-int(r.start)+1),r.database,r.seqtype,r.field6,r.strand,r.field8,r.info])+"\n")