import re
import string
import sys
import fileinput

class Record1(object):
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
		self.uniqid = fields[10].split(";")[9].split(" ")[2]

# class Record2(object):
# 	def __init__(self,fields):
# 		self.chr = fields[0]
# 		self.start = fields[1]
# 		self.end = fields[2]
# 		self.field4 = fields[3]
# 		self.strand = fields[4]
# 		self.field6 = fields[5]
# 		self.info = fields[6]
# 		if("ccds_id" in self.info):
# 			m=re.search(" ccds_id.*?;",self.info)
# 			n=re.sub("ccds_id.*?; ","",self.info)
# 			self.info=n+m.group(0)
# 		if("tag" in self.info):
# 			m=re.search(" tag.*?;",self.info)
# 			n=re.sub("tag.*?; ","",self.info)
# 			self.info=n+m.group(0)
# 		foo=self.info.rstrip().split(";")
# 		bar=(foo[2]+foo[0]+foo[8]).split("\"")
# 		self.uniqid=bar[3]+"_"+bar[5]

f1=open(sys.argv[2],"w")
# f2=open(sys.argv[4],"w+")

for line in fileinput.input(sys.argv[1]):
	fields=line.rstrip().split("\t")
	try:
		old_r=r
	except NameError:
		pass
	r=Record1(fields)
	try:
		if re.search("transcript_id (.*?);",r.info).group(1).strip("\"") != r.transid:
			continue
	except NameError:
		pass
	try:
		if r.uniqid == old_r.uniqid or re.search("pfamA_auto_reg.*?;",r.info).group(0) == re.search("pfamA_auto_reg.*?;",old_r.info).group(0):
			continue
	except NameError:
		pass
	f1.write("\t".join([r.chr,r.start,r.end,str(int(r.end)-int(r.start)+1),r.database,r.seqtype,r.field7,r.strand,r.field9,"uniq_id "+r.uniqid+"; "+r.info])+"\n")

# del old_r,r

# ct=1;

# for line in fileinput.input(sys.argv[2]):
# 	fields=line.rstrip().split("\t")
# 	try:
# 		old_r=r
# 	except NameError:
# 		pass
# 	r=Record2(fields)
# 	try:
# 		if re.search("gene_id (.*?);",r.info).group(1).strip("\"") == re.search("gene_id (.*?);",old_r.info).group(1).strip("\""):
# 			ct=ct+1;
# 		else:
# 			ct=1;
# 	except NameError:
# 		pass
# 	r.uniqid=r.uniqid+"_ND"+str(ct)
# 	try:
# 		if re.search("exon_id.*?;",r.info).group(0) == re.search("exon_id.*?;",old_r.info).group(0) and re.search("gene_id.*?;",r.info).group(0) == re.search("gene_id.*?;",old_r.info).group(0):
# 			continue
# 	except NameError:
# 		pass
# 	f2.write("\t".join([r.chr,r.start,r.end,str(int(r.end)-int(r.start)+1),r.field4,r.strand,r.field6,"uniq_id "+r.uniqid+"; "+r.info])+"\n")

f1.close();#f2.close()
