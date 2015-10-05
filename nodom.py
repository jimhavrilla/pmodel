import string
import sys
import pybedtools
import re

class RecordExon(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.info = "\t".join(fields[6:24]).rstrip()
		self.uniqid = fields[7].strip(";").strip("\"")+"_"+fields[23].rstrip().strip(";").strip("\"")+"_NoDom"
		self.transid = fields[9].strip(";").strip("\"")

class RecordUniq(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.gene = fields[3]
		self.autoreg = fields[4]
		self.transid = fields[5]
		self.exonid = fields[6]

f1=open(sys.argv[1],"r")
f2=open(sys.argv[2],"r")
exons=f1.readline()
uniqs=f2.readline()
e_bed=''
u_bed=''
utrans=set()

while exons and uniqs:
	e=RecordExon(exons.split("\t"))
	u=RecordUniq(uniqs.split("\t"))
	utrans.add(u.transid)
	old_e_trans=e.transid
	old_u_trans=u.transid
	if e.transid == u.transid:
		while e.transid == old_e_trans:
			e_bed=e_bed+e.chr+"\t"+e.start+"\t"+e.end+"\t"+e.info+"\t"+e.uniqid+"\n"
			exons=f1.readline()
			try:
				e=RecordExon(exons.split("\t"))
			except IndexError:
				old_e_trans=None
		while u.transid == old_u_trans:
			u_bed=u_bed+u.chr+"\t"+u.start+"\t"+u.end+"\n"
			uniqs=f2.readline()
			try:
				u=RecordUniq(uniqs.split("\t"))
			except IndexError:
				old_u_trans=None
		a=pybedtools.BedTool(e_bed,from_string=True)
		b=pybedtools.BedTool(u_bed,from_string=True)
		bed=str(a.subtract(b))
		bed2=str(a.intersect(b,v=True))
		bed=bed+"\n"+bed2
		bed=str(pybedtools.BedTool(bed,from_string=True).sort())
		old_chrom=None
		old_start=None
		old_end=None
		old_uniq=None
		ct=1
		for i in bed.splitlines():
			i=i.split()
			chrom=i[0];start=i[1];end=i[2];info="\t".join(i[3:len(i)-1]);uniqid=i[len(i)-1]
			if old_chrom!=None and old_uniq==uniqid and not (chrom==old_chrom and start==old_start and end==old_end):
				ct=ct+1
			else:
				ct=1
			old_chrom=chrom
			old_start=start
			old_end=end
			old_uniq=uniqid
			sys.stdout.write("\t".join([chrom,start,end,str(int(end)-int(start)),info+" uniqid "+uniqid+str(ct)])+"\n")
		e_bed=''
		u_bed=''
	elif e.transid > u.transid:
		uniqs=f2.readline()
	elif e.transid < u.transid:
		exons=f1.readline()
	if not uniqs or not exons:
		break

f2.close()

f1.seek(0)
exons=f1.readline()
e_bed=''

while exons:
	e=RecordExon(exons.split("\t"))
	old_e_trans=e.transid
	if e.transid not in utrans:
		while e.transid == old_e_trans:
			e_bed=e_bed+e.chr+"\t"+e.start+"\t"+e.end+"\t"+e.info+"\t"+e.uniqid+"\n"
			exons=f1.readline()
			try:
				e=RecordExon(exons.split("\t"))
			except IndexError:
				old_e_trans=None
		bed=str(pybedtools.BedTool(e_bed,from_string=True).sort())
		old_chrom=None
		old_start=None
		old_end=None
		old_uniq=None
		ct=1
		for i in bed.splitlines():
			i=i.split()
			chrom=i[0];start=i[1];end=i[2];info="\t".join(i[3:len(i)-1]);uniqid=i[len(i)-1]
			if old_chrom!=None and old_uniq==uniqid and not (chrom==old_chrom and start==old_start and end==old_end):
				ct=ct+1
			else:
				ct=1
			old_chrom=chrom
			old_start=start
			old_end=end
			old_uniq=uniqid
			sys.stdout.write("\t".join([chrom,start,end,str(int(end)-int(start)),info+" uniqid "+uniqid+str(ct)])+"\n")
		e_bed=''
	else:
		exons=f1.readline()
	if not exons:
		break

f1.close()
