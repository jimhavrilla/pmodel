#nodom.py - gets nodom regions
import string
import sys
import pybedtools
import re

class Record1(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.info = " ".join(fields[6:24])
		self.info=re.sub("\"|;","",self.info)
		if len(fields)>24:
			self.chr2 = fields[len(fields)-11]
			self.start2 = fields[len(fields)-10]
			self.end2 = fields[len(fields)-9]
			#self.info=self.info+" near_"+" ".join([fields[len(fields)-2],fields[len(fields)-1]])
		else:
			self.chr2 = None
			self.start2 = None
			self.end2 = None
		self.uniqid = fields[7].strip(";").strip("\"")+"_"+fields[23].strip(";").strip("\"")+"_ND"

f2=open(sys.argv[1],"w")
old_r=None
ct=0
ct2=0
for line in sys.stdin:
	fields=line.rstrip().split("\t")
	r_=Record1(fields)
	if old_r!=None and old_r.start+old_r.end != r_.start+r_.end:
		if old_r.uniqid==r_.uniqid:
			ct=ct+1
		else:
			ct=0
		if old_r.chr2!=None and old_r.start2!=None and old_r.end2!=None:
			a=pybedtools.BedTool(old_r.chr+"\t"+old_r.start+"\t"+old_r.end,from_string=True)
			b=pybedtools.BedTool(old_r.chr2+"\t"+old_r.start2+"\t"+old_r.end2,from_string=True)
			bed=str(a.subtract(b))
		else:
			bed=old_r.chr+"\t"+old_r.start+"\t"+old_r.end+"\n"
		for i in bed.splitlines():
			ct2=ct2+1
			i=i.split("\t")
			old_r.chr=i[0];old_r.start=i[1];old_r.end=i[2]
			f2.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)+1),old_r.info+" uniqid "+old_r.uniqid+str(ct+ct2)])+"\n")
	ct2=0
	old_r=r_

if old_r.chr2!=None and old_r.start2!=None and old_r.end2!=None:
	a=pybedtools.BedTool(old_r.chr+"\t"+old_r.start+"\t"+old_r.end,from_string=True)
	b=pybedtools.BedTool(old_r.chr2+"\t"+old_r.start2+"\t"+old_r.end2,from_string=True)
	bed=str(a.subtract(b))
else:
	bed=old_r.chr+"\t"+old_r.start+"\t"+old_r.end+"\n"
for i in bed.splitlines():
	ct2=ct2+1
	i=i.split("\t")
	old_r.chr=i[0];old_r.start=i[1];old_r.end=i[2]
	f2.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)+1),old_r.info+" uniqid "+old_r.uniqid+str(ct+ct2)])+"\n")

f2.close();