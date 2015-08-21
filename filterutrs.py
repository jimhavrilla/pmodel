import string
import sys
import pybedtools

class RecordExon(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.crapandstrand = "\t".join(fields[3:6])
		self.info=fields[6]
		self.transid = self.info.replace("\"","").replace(";","").split(" ")[3]

class RecordUtrs(object):
	def __init__(self,fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.crapandstrand = "\t".join(fields[3:6])
		self.info = fields[6]
		self.transid = self.info.replace("\"","").replace(";","").split(" ")[3]

f1=open(sys.argv[1],"r")
f2=open(sys.argv[2],"r")
exons=f1.readline()
utrs=f2.readline()
e_bed=''
u_bed=''

while exons and utrs:
	e=RecordExon(exons.split("\t"))
	u=RecordUtrs(utrs.split("\t"))
	old_e_trans=e.transid
	old_u_trans=u.transid
	if e.transid == u.transid:
		while e.transid == old_e_trans:
			e_bed=e_bed+e.chr+"\t"+e.start+"\t"+e.end+"\t"+e.crapandstrand+"\t"+e.info+"\n"
			exons=f1.readline()
			try:
				e=RecordExon(exons.split("\t"))
			except IndexError:
				old_e_trans=None
		while u.transid == old_u_trans:
			u_bed=u_bed+u.chr+"\t"+u.start+"\t"+u.end+"\n"
			utrs=f2.readline()
			try:
				u=RecordUtrs(utrs.split("\t"))
			except IndexError:
				old_u_trans=None
		a=pybedtools.BedTool(e_bed,from_string=True)
		b=pybedtools.BedTool(u_bed,from_string=True)
		bed=str(a.subtract(b))
		old_chrom=None
		old_start=None
		old_end=None
		for i in bed.splitlines():
			i=i.split()
			chrom=i[0];start=i[1];end=i[2];crapandstrand="\t".join(i[3:6]);info=" ".join(i[6:len(i)])
			old_chrom=chrom
			old_start=start
			old_end=end
			sys.stdout.write("\t".join([chrom,start,end,crapandstrand,info])+"\n")
		e_bed=''
		u_bed=''
	elif e.transid > u.transid:
		utrs=f2.readline()
	elif e.transid < u.transid:
		exons=f1.readline()
	if not utrs or not exons:
		break

f1.close()
f2.close()