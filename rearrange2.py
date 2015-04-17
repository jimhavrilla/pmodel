import string
import sys
import fileinput
import pybedtools
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d","--domain",dest="bool",help="domains",action="store_false")
parser.add_option("-u","--uniqid",dest="bool",help="default,uniqids",action="store_true")
(options,args)=parser.parse_args()

if options.bool==True:
	class Record1(object):
		def __init__(self, fields):
			self.chr = fields[0]
			self.start = fields[1]
			self.end = fields[2]
			self.info = " ".join(fields[9:len(fields)])
			self.interval = fields[3]
			self.database = fields[4]
			self.seqtype = fields[5]
			self.field7 = fields[6]
			self.field9 = fields[8]
			self.strand = fields[7]
			self.transid = fields[14]
			self.transid2 = fields[28]
			self.uniqid = fields[44].rstrip(';')
			self.autoreg = fields[24]
			self.gene = fields[12]

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
	f1=open(sys.argv[2],"r")
	f2=open(sys.argv[3],"w")
	bed=''
	old_r=None
	for line in f1:
		fields=line.rstrip().split("\t")
		r_=Record1(fields)
		if old_r != None and r_.autoreg != old_r.autoreg:
			bed=str(pybedtools.BedTool(bed.rstrip('\n'),from_string=True).merge())
			for i in bed.splitlines():
				i=i.split("\t")
				old_r.chr=i[0];old_r.start=i[1];old_r.end=i[2]
				f2.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)),old_r.database,old_r.seqtype,old_r.field7,old_r.strand,old_r.field9,old_r.info])+"\n")
			bed=''
		bed=bed+r_.chr+"\t"+r_.start+"\t"+r_.end+"\n"
		old_r=r_

	bed=str(pybedtools.BedTool(bed.rstrip('\n'),from_string=True).merge())
	for i in bed.splitlines():
		i=i.split("\t")
		old_r.chr=i[0];old_r.start=i[1];old_r.end=i[2]
		f2.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)),old_r.database,old_r.seqtype,old_r.field7,old_r.strand,old_r.field9,old_r.info])+"\n")

	f1.close();f2.close();

if options.bool==False:
	class Record1(object):
		def __init__(self, fields):
			self.chr = fields[0]
			self.start = fields[1]
			self.end = fields[2]
			self.info = " ".join(fields[9:len(fields)])
			self.interval = fields[3]
			self.database = fields[4]
			self.seqtype = fields[5]
			self.field7 = fields[6]
			self.field9 = fields[8]
			self.domain = fields[10]
			self.strand = fields[7]
			self.transid = fields[14]
			self.transid2 = fields[28]
			self.uniqid = fields[44].rstrip(';')
			self.autoreg = fields[24]

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
				f2.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)),old_r.database,old_r.seqtype,old_r.field7,old_r.strand,old_r.field9,old_r.info])+"\n")
			bed=''
		bed=bed+r_.chr+"\t"+r_.start+"\t"+r_.end+"\n"
		old_r=r_

	bed=str(pybedtools.BedTool(bed.rstrip('\n'),from_string=True).merge())
	for i in bed.splitlines():
		i=i.split("\t")
		old_r.chr=i[0];old_r.start=i[1];old_r.end=i[2]
		f2.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)),old_r.database,old_r.seqtype,old_r.field7,old_r.strand,old_r.field9,old_r.info])+"\n")

	f2.close();