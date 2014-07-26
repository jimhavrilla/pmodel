# mergetable.py - to make a table for counts vs. dn/ds plot and bp vs. var/bp plot
import string
import sys

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

f1=open(sys.argv[1],"r")#f1=open("foo.txt","r")
f2=open(sys.argv[2],"r")#f2=open("human_pfam.counts","r")
f3=open(sys.argv[3],"r")#f3=open("sumlist.bed","r")
table=f1.readline()
count=f2.readline()
totsum=f3.readline()
sys.stdout.write("\t".join(["domain","gene","nonsynct","synct","totalvarct","domcount","totalbp","mmaf","dn/ds","var/bp ratio"])+"\n")
while table and count and totsum:
	a=RecordA(table.split("\t"))
	b=RecordB(count.split("\t"))
	c=RecordC(totsum.split(";"))
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
		c=RecordC(totsum.split(";"))
	if b.domain<c.domain:
		count=f2.readline()
		b=RecordB(count.split("\t"))
	if b.domain>c.domain:
		totsum=f3.readline()
		c=RecordC(totsum.split(";"))
	if a.domain==b.domain==c.domain:
		print "\t".join([a.domain,a.gene,a.nonsynct,a.synct,a.totalct,b.count,c.totalbp,a.mmaf,a.dnds,str(round(float(a.totalct)/float(c.totalbp),2))])
		table=f1.readline()
		count=f2.readline()
		totsum=f3.readline()