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

def grouper(o):
    return o.autoreg, o.transid

grouper = lambda o: (o.autoreg, o.transid)

from itertools import groupby

for grp, records in groupby((Record1(line.rstrip().split()) for line in f1), grouper):
    records = list(records)

    bed = "\n".join(r.chr+"\t"+r.start+"\t"+r.end+"\t"+r.database+"\t"+r.seqtype+"\t"+r.field7+"\t"+r.strand+"\t"+r.field9+"\t"+r.info1+"\t"+r.transid+"\t"+r.info2+"\t"+r.autoreg+"\t"+r.info3+"\t"+r.exonid+"\t"+r.info4 for r in records)

    bed = str(pybedtools.BedTool(bed, from_string=True))

    for i in (x.split() for x in bed.splitlines()):

         domain=i[9];transid=i[13];autoreg=i[23];geneid=i[25];exonid=i[41]
         uniqid=geneid+"_"+transid+"_"+exonid+"_"+autoreg+"_"+domain
         sys.stdout.write("\t".join(i[0:3])+"\t"+str(int(i[2])-int(i[1]))+"\t"+"\t".join(i[3:43])+"\t"+uniqid+"\n")

f1.close();
