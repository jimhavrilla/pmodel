#nodom.py
import string
import sys
import pybedtools

class Record1(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.info = " ".join(fields[9:len(fields)-11])+" ".join(fields[len(fields)-1:len(fields)])
		self.chr2 = fields[len(fields)-10]
		self.start2 = fields[len(fields)-9]
		self.end2 = fields[len(fields)-8]

f2=open(sys.argv[1],"w")
bed=''
old_r=None
for line in sys.stdin:
	fields=line.rstrip().split("\t")
	r_=Record1(fields)
	if old_r != None and r_.domain != old_r.domain:
		bed=str(pybedtools.BedTool(bed.rstrip('\n'),from_string=True).merge())
		for i in bed.splitlines():
			i=i.split("\t")
			old_r.chr=i[0];old_r.start=i[1];old_r.end=i[2]
			f2.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)+1),old_r.database,old_r.seqtype,old_r.field7,old_r.strand,old_r.field9,old_r.info])+"\n")
		bed=''
	bed=bed+r_.chr+"\t"+r_.start+"\t"+r_.end+"\n"
	old_r=r_

bed=str(pybedtools.BedTool(bed.rstrip('\n'),from_string=True).merge())
for i in bed.splitlines():
	i=i.split("\t")
	old_r.chr=i[0];old_r.start=i[1];old_r.end=i[2]
	f2.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)+1),old_r.database,old_r.seqtype,old_r.field7,old_r.strand,old_r.field9,old_r.info])+"\n")

f2.close();