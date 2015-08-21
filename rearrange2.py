import string
import sys
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
			self.info1 = " ".join(fields[9:14])
			self.info2 = " ".join(fields[15:24])
			self.info3 = " ".join(fields[25:42])
			self.info4 = " ".join(fields[43:44])
			self.info5 = " ".join(fields[45:len(fields)])
			self.interval = fields[3]
			self.database = fields[4]
			self.seqtype = fields[5]
			self.field7 = fields[6]
			self.field9 = fields[8]
			self.strand = fields[7]
			self.domain = fields[10]
			self.transid = fields[14]
			self.transid2 = fields[28]
			self.uniqid = fields[44].rstrip(';')
			self.autoreg = fields[24]
			self.geneid = fields[26]
			self.gene = fields[32]
			self.exonid = fields[42]
			
	f1=open(sys.argv[2],"r")
	bed=''
	old_r=None
	olduniq=None
	ct=1
	for line in f1:
		fields=line.rstrip().split("\t")
		r_=Record1(fields)
		if old_r != None and (r_.domain != old_r.domain or r_.geneid!=old_r.geneid): # check why ensl ids are different and why transcript pfam != transcript ensl
			bed=str(pybedtools.BedTool(bed.rstrip('\n'),from_string=True).merge(o=['distinct','distinct','distinct'],c=[4,5,6]))
			for i in bed.splitlines():
				i=i.split("\t")
				old_r.chr=i[0];old_r.start=i[1];old_r.end=i[2];old_r.transid=i[3];old_r.exonid=i[4];old_r.autoreg=i[5]
				if olduniq != None and olduniq == old_r.geneid+"_"+old_r.exonid+"_"+old_r.autoreg+"_"+old_r.domain:
					ct=ct+1 # tracks potential duplicates
				else:
					olduniq=old_r.geneid+"_"+old_r.exonid+"_"+old_r.autoreg+"_"+old_r.domain
					ct=1
				sys.stdout.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)),old_r.database,old_r.seqtype,old_r.field7,old_r.strand,old_r.field9," ".join([old_r.info1,old_r.transid,old_r.info2,old_r.autoreg,old_r.info3,old_r.exonid,old_r.info4,old_r.geneid+"_"+old_r.exonid+"_"+old_r.autoreg+"_"+old_r.domain+"_"+str(ct),old_r.info5])])+"\n")
			bed=''
		bed=bed+r_.chr+"\t"+r_.start+"\t"+r_.end+"\t"+r_.transid2+"\t"+r_.exonid+"\t"+r_.autoreg+"\n" # transid2, because ENST is always correct from ensembl file, not pfam
		old_r=r_

	bed=str(pybedtools.BedTool(bed.rstrip('\n'),from_string=True).merge(o=['distinct','distinct','distinct'],c=[4,5,6]))
	for i in bed.splitlines():
		i=i.split("\t")
		old_r.chr=i[0];old_r.start=i[1];old_r.end=i[2];old_r.transid=i[3];old_r.exonid=i[4];old_r.autoreg=i[5]
		if olduniq != None and olduniq == old_r.geneid+"_"+old_r.exonid+"_"+old_r.autoreg+"_"+old_r.domain:
			ct=ct+1
		else:
			olduniq=old_r.geneid+"_"+old_r.exonid+"_"+old_r.autoreg+"_"+old_r.domain
			ct=1
		sys.stdout.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)),old_r.database,old_r.seqtype,old_r.field7,old_r.strand,old_r.field9," ".join([old_r.info1,old_r.transid,old_r.info2,old_r.autoreg,old_r.info3,old_r.exonid,old_r.info4,old_r.geneid+"_"+old_r.exonid+"_"+old_r.autoreg+"_"+old_r.domain+"_"+str(ct),old_r.info5])])+"\n")

	f1.close();

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
				sys.stdout.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)),old_r.database,old_r.seqtype,old_r.field7,old_r.strand,old_r.field9,old_r.info])+"\n")
			bed=''
		bed=bed+r_.chr+"\t"+r_.start+"\t"+r_.end+"\n"
		old_r=r_

	bed=str(pybedtools.BedTool(bed.rstrip('\n'),from_string=True).merge())
	for i in bed.splitlines():
		i=i.split("\t")
		old_r.chr=i[0];old_r.start=i[1];old_r.end=i[2]
		sys.stdout.write("\t".join([old_r.chr,old_r.start,old_r.end,str(int(old_r.end)-int(old_r.start)),old_r.database,old_r.seqtype,old_r.field7,old_r.strand,old_r.field9,old_r.info])+"\n")

	f2.close();