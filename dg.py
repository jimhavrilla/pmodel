# dg.py - 
import string
import sys
import numpy
import warnings
from argparse import ArgumentParser

warnings.simplefilter("error") #catches NaN warnings and empty slices

parser = ArgumentParser()
parser.add_argument("--domain","-d",help="run the domain protocol",action="store_true")
parser.add_argument("--uniqid","-u",help="run the uniqid protocol",action="store_true")
parser.add_argument("--domcount","-c",help="default=>0 domcount",default=0,type=int)
parser.add_argument("--bp","-b",help="default=>0 bp",default=0,type=int)
parser.add_argument("--files","-f",help="three input files, the maketable results, the pfam domain counts, and the sum of bp for each domain, and the output filename",nargs="*",default=['allint2.bed','human_pfam.counts','sumlist.bed','dngpair.txt'],type=str)
args=parser.parse_args()

dpar=args.domcount
bpar=args.bp

if args.uniqid:
	class Record(object):
		def __init__(self, fields):
			self.domain = fields[0]
			self.gene = fields[1]
			self.autoreg = fields[2]
			self.uniqid = fields[3]
			self.coverage = fields[4]
			self.truelength = fields[5]
			self.covratio = fields[6]
			self.nct = fields[7]
			self.sct = fields[8]
			self.ct = fields[9]
			self.dnds = float(fields[10])
			self.density = fields[11]
			self.mmaf = fields[12]
			self.domcount = fields[13]

	xlist=[]
	old_r=None
	for line in sys.stdin:
		fields=line.rstrip().split(" ")
		r_=Record(fields)
		if old_r!=None and old_r.gene != r_.gene:
			l=[y.dnds for y in xlist[0:len(xlist)]]
			m=numpy.mean(l)
			s=numpy.std(l)
			for y in xlist:
				if s==0.0:
					s=1.0
				if bpar>r_.truelength or dpar>r_.domcount:
					continue
				else:
					sys.stdout.write("\t".join([y.domain,y.gene,y.autoreg,y.uniqid,y.coverage,y.covratio,y.truelength,y.nct,y.sct,y.ct,str(y.dnds),str((y.dnds-m)/s),y.density,y.mmaf,y.domcount,"\n"]))
			xlist=[]
		xlist.append(r_)
		old_r=r_

	l=[y.dnds for y in xlist[0:len(xlist)]]
	m=numpy.mean(l)
	s=numpy.std(l)
	for y in xlist:
		if s==0.0:
			s=1.0
		if bpar>r_.truelength or dpar>r_.domcount:
			continue
		else:
			sys.stdout.write("\t".join([y.domain,y.gene,y.autoreg,y.uniqid,y.coverage,y.covratio,y.truelength,y.nct,y.sct,y.ct,str(y.dnds),str((y.dnds-m)/s),y.density,y.mmaf,y.domcount,"\n"]))

if args.domain:
	class RecordA(object):
		def __init__(self, fields):
			self.chr = fields[0]
			self.start = fields[1]
			self.end = fields[2]
			self.ref = fields[3]
			self.alt = fields[4]
			self.domain = fields[5].rstrip(";").strip("\"")
			self.uniqid = fields[6]
			self.gene = fields[7]
			self.maf = float(fields[8])
			self.impact = fields[9]
			self.type = fields[10]
			self.info = fields[11:21]
	class RecordB(object):
		def __init__(self,fields):
			self.pfamacc = fields[0]
			self.domain = fields[1]
			self.count = fields[2].strip()
	class RecordC(object):
		def __init__(self,fields):
			self.domain = fields[0].strip("\"")
			self.totalbp = fields[1].strip()

	def ratiocalc(nct,sct):
		try:
			ns=round(float(nct)/float(sct),2)
		except ZeroDivisionError:
			ns=round(float(nct)/float(sct+1),2)

		return ns

	def count(old_r,maf,ct,nct,sct):
		mmaf=numpy.median(maf)
		table.append([old_r.domain,old_r.gene,str(ct),str(nct),str(sct),ratiocalc(nct,sct),str(mmaf)])
		return table

	ct=0
	nct=0
	sct=0
	maf=[]
	table=[]
	old_r=None
	f1=open(args.files[0],"r");f2=open(args.files[1],"r");f3=open(args.files[2],"r")
	for line in f1:
		fields=line.rstrip().split("\t")
		r_=RecordA(fields)
		if old_r!=None and old_r.gene!=r_.gene: # assumes data is sorted by domain then gene, which it should be
			count(old_r,maf,ct,nct,sct)
			nct=0;sct=0;ct=0;
			maf=[]
		if r_.type=="ds":
			sct=sct+1
			ct=ct+1
			maf.append(r_.maf)
		if r_.type=="dn":
			nct=nct+1
			ct=ct+1
			maf.append(r_.maf)
		
		old_r=r_

	count(r_,maf,ct,nct,sct)

	i=0
	count=f2.readline()
	totsum=f3.readline()

	while count and totsum:
		a=table[i][0]
		b=RecordB(count.split("\t"))
		c=RecordC(totsum.split(" "))
		if a<b.domain:
			i=i+1
			a=table[i][0]
		if a>b.domain:
			count=f2.readline()
			b=RecordB(count.split("\t"))
		if a<c.domain:
			i=i+1
			a=table[i][0]
		if a>c.domain:
			totsum=f3.readline()
			c=RecordC(totsum.split(" "))
		if b.domain<c.domain:
			count=f2.readline()
			b=RecordB(count.split("\t"))
		if b.domain>c.domain:
			totsum=f3.readline()
			c=RecordC(totsum.split(" "))
		if a==b.domain==c.domain:
			if int(b.count) < dpar:
				count=f2.readline()
				continue
			if int(c.totalbp) < bpar:
				totsum=f3.readline()
				continue
			table[i].extend((b.count,c.totalbp))
			try:
				while table[i+1][0]==a:
					i=i+1
					table[i].append(b.count)
					table[i].append(c.totalbp)
				else:
					i=i+1
					count=f2.readline()
					totsum=f3.readline()
			except IndexError:
				break

	f4=open(args.files[3],"w")
	xlist=[]
	old_x=None #x[0] domain, x[1] gene, x[2] ct, x[3] nct, x[4] sct, x[5] dn/ds, x[6] mmaf
	for x in table:
		if old_x!=None and old_x[0] != x[0]:
			l=[y[5] for y in xlist[0:len(xlist)]]
			m=numpy.mean(l)
			s=numpy.std(l)
			for y in xlist:
				if s==0.0:
					s=1.0
				if len(y)<8 and bpar==0 and dpar==0:
					f4.write("\t".join([y[0],y[1],y[2],y[3],y[4],str(y[5]),str((y[5]-m)/s),y[6],"0","0","\n"])) #nodom, gene, ct, nct, sct, dn/ds, z-score, mmaf
				elif len(y)<8 and (bpar>0 or dpar>0):
					continue
				else:
					f4.write("\t".join([y[0],y[1],y[2],y[3],y[4],str(y[5]),str((y[5]-m)/s),y[6],y[7],y[8],"\n"])) #domain, gene, ct, nct, sct, dn/ds, z-score, mmaf, domcount, totalbp
			xlist=[]
		xlist.append(x)
		old_x=x

	l=[y[5] for y in xlist[0:len(xlist)]]
	m=numpy.mean(l)
	s=numpy.std(l)
	for y in xlist:
		if s==0.0:
			s=l[0]
			if len(y)<8 and bpar==0 and dpar==0:
				f4.write("\t".join([y[0],y[1],y[2],y[3],y[4],str(y[5]),str((y[5]-m)/s),y[6],"0","0","\n"])) #nodom, gene, ct, nct, sct, dn/ds, z-score, mmaf
			elif len(y)<8 and (bpar>0 or dpar >0):
				continue
			else:
				f4.write("\t".join([y[0],y[1],y[2],y[3],y[4],str(y[5]),str((y[5]-m)/s),y[6],y[7],y[8],"\n"])) #domain, gene, ct, nct, sct, dn/ds, z-score, mmaf, domcount, totalbp