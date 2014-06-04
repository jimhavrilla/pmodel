import string
import sys

f1=open(sys.argv[1],"r")
f2=open(sys.argv[2],"r")

count=f1.readline()
totsum=f2.readline()

while count or totsum:
	c=count.split()
	t=totsum.split()
	if c[0]<t[0]:
		count=f1.readline()
	elif c[0]>t[0]:
		totsum=f2.readline()
	elif c[0]==t[0]:
		print c[0]+"\t"+str(float(c[1])/float(t[1]))
		count=f1.readline()
		totsum=f2.readline()