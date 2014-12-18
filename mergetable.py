# mergetable.py - to make a table for counts vs. dn/ds plot and bp vs. var/bp plot
import string
import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--domcount","-d",help="default=>0 domcount",default=0,type=int)
parser.add_argument("--bp","-b",help="default=>0 bp",default=0,type=int)
parser.add_argument("--files","-f",help="three files, the maketable results, the pfam domain counts, and the sum of bp for each domain",nargs=3,default=['foo.bed','human_pfam.counts','sumlist.bed'],type=str)
args=parser.parse_args()

dpar=args.domcount
bpar=args.bp

class RecordA(object):
	def __init__(self, fields):
		self.domain = fields[0]
		self.gene = fields[1]
		self.nonsynct = fields[2]
		self.synct = fields[3]
		self.totalct = fields[4]
		self.dnds = fields[5]
		self.mmaf = fields[6].strip()
class RecordB(object):
	def __init__(self,fields):
		self.pfamacc = fields[0]
		self.domain = fields[1]
		self.count = fields[2].strip()
class RecordC(object):
	def __init__(self,fields):
		self.domain = fields[0].strip("\"")
		self.totalbp = fields[1].strip()

f1=open(args.files[0],"r")
f2=open(args.files[1],"r")
f3=open(args.files[2],"r")
table=f1.readline()
count=f2.readline()
totsum=f3.readline()
sys.stdout.write("\t".join(["domain","gene","nonsynct","synct","totalvarct","domcount","totalbp","mmaf","dn/ds","var/bp ratio"])+"\n")
while table and count and totsum:
	a=RecordA(table.split("\t"))
	b=RecordB(count.split("\t"))
	c=RecordC(totsum.split(" "))
	if a.domain<b.domain:
		table=f1.readline()
		a=RecordA(table.split("\t"))
	if a.domain>b.domain:
		count=f2.readline()
		b=RecordB(count.split("\t"))
	if a.domain<c.domain:
		table=f1.readline()
		a=RecordA(table.split("\t"))
	if a.domain>c.domain:
		totsum=f3.readline()
		c=RecordC(totsum.split(" "))
	if b.domain<c.domain:
		count=f2.readline()
		b=RecordB(count.split("\t"))
	if b.domain>c.domain:
		totsum=f3.readline()
		c=RecordC(totsum.split(" "))
	if a.domain==b.domain==c.domain:
		if int(b.count) < dpar:
			count=f2.readline()
			continue
		if int(c.totalbp) < bpar:
			totsum=f3.readline()
			continue
		print "\t".join([a.domain,a.gene,a.nonsynct,a.synct,a.totalct,b.count,c.totalbp,a.mmaf,a.dnds,str(float(a.totalct)/float(c.totalbp))])
		table=f1.readline()
		count=f2.readline()
		totsum=f3.readline()