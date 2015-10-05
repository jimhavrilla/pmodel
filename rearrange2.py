import string
import sys
import pybedtools

class Record1(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.info1 = "\t".join(fields[9:14])
		self.info2 = "\t".join(fields[15:24])
		self.info3 = "\t".join(fields[25:42])
		self.info4 = "\t".join(fields[43:44])
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
		
f1=open(sys.argv[1],"r")
bed=''
old_r=None
olduniq=None
ct=1
rangerep=[['distinct']*len(range(4,44)),range(4,44)]
for line in f1:
	fields=line.rstrip().split()
	r_=Record1(fields)
	#if r_.transid == r_.transid2:
	if old_r != None and (r_.autoreg != old_r.autoreg or r_.transid!=old_r.transid): # check why ensl ids are different and why transcript pfam != transcript ensl
		bed=str(pybedtools.BedTool(bed.rstrip('\n'),from_string=True).merge(o=rangerep[0],c=rangerep[1]))  # in the future, just in case add merging on all fields
		for i in bed.splitlines():
			i=i.split()
			domain=i[9];transid=i[13];autoreg=i[23];geneid=i[25];exonid=i[41]
			uniqid=geneid+"_"+transid+"_"+exonid+"_"+autoreg+"_"+domain
			sys.stdout.write("\t".join(i[0:3])+"\t"+str(int(i[2])-int(i[1]))+"\t"+"\t".join(i[3:43])+"\t"+uniqid+"\n")
		bed=''
	bed=bed+r_.chr+"\t"+r_.start+"\t"+r_.end+"\t"+r_.database+"\t"+r_.seqtype+"\t"+r_.field7+"\t"+r_.strand+"\t"+r_.field9+"\t"+r_.info1+"\t"+r_.transid+"\t"+r_.info2+"\t"+r_.autoreg+"\t"+r_.info3+"\t"+r_.exonid+"\t"+r_.info4+"\n" # in the future, just in case add merging on all fields
	old_r=r_


bed=str(pybedtools.BedTool(bed.rstrip('\n'),from_string=True).merge(o=rangerep[0],c=rangerep[1])) # in the future, just in case add merging on all fields
for i in bed.splitlines():
	i=i.split()
	domain=i[9];transid=i[13];autoreg=i[23];geneid=i[25];exonid=i[41]
	uniqid=geneid+"_"+transid+"_"+exonid+"_"+autoreg+"_"+domain
	sys.stdout.write("\t".join(i[0:3])+"\t"+str(int(i[2])-int(i[1]))+"\t"+"\t".join(i[3:43])+"\t"+uniqid+"\n")


f1.close();