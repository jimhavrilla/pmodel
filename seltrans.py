import sys
import string

class Appris(object):
	def __init__(self,fields):
		self.gene=fields[0]
		self.geneid=fields[1]
		self.transid=fields[2]
		self.ccdsid=fields[3]
		self.type=fields[4]

class Transcript(object):
	def __init__(self,fields):
		self.geneid=fields[0]
		self.transid=fields[1]
		self.gene=fields[2]
		self.length=fields[3]

f1=open(sys.argv[1],"r")
f2=open(sys.argv[2],"r")
appris=f1.readline().rstrip()
trans=f2.readline().rstrip()
principal=[]
maxlen=None

while appris:
	a=Appris(appris.split("\t"))
	t=Transcript(trans.split("\t"))
	old_a_gene=a.geneid
	old_t_gene=t.geneid
	if a.geneid == t.geneid:
		while a.geneid == old_a_gene and t.geneid == old_t_gene:
			if a.transid == t.transid:
				principal.append([a.geneid,a.transid,t.length])
				appris=f1.readline().rstrip()
				trans=f2.readline().rstrip()
				try:
					t=Transcript(trans.split("\t"))
				except IndexError:
					old_t_gene=None
				try:
					a=Appris(appris.split("\t"))
				except IndexError:
					old_a_gene=None
			elif a.transid > t.transid:
				trans=f2.readline().rstrip()
				try:
					t=Transcript(trans.split("\t"))
				except IndexError:
					old_t_gene=None
			elif a.transid < t.transid:
				appris=f1.readline().rstrip()
				try:
					a=Appris(appris.split("\t"))
				except IndexError:
					old_a_gene=None
		if principal!=[]:
			for i in range(len(principal)):
				if maxlen==None or principal[i][2]>maxlen[2]:
					maxlen=principal[i]
			sys.stdout.write("\t".join(maxlen)+"\n")
		elif principal==[] and a.geneid<t.geneid:
			sys.stdout.write("\t".join([a.geneid,a.transid,"-"])+"\n")
		maxlen=None
		principal=[]
	elif a.geneid > t.geneid:
		trans=f2.readline().rstrip()
	elif a.geneid < t.geneid:
		appris=f1.readline().rstrip()
	if appris and not trans:
		while appris:
			a=Appris(appris.split("\t"))
			sys.stdout.write("\t".join([a.geneid,a.transid,"-"])+"\n")
			appris=f1.readline().rstrip()
	if not appris and not trans:
		break