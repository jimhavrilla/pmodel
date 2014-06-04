import re
import string
import sys
import fileinput

for line in fileinput.input():
	foo=line.rstrip().split("\t")
	foo[8]=foo[8]+"; "
	a=foo[1]
	b=foo[2]
	if int(foo[3])>intq(foo[4]):
		foo[2]=foo[3]
		foo[1]=foo[4]
	else:
		foo[1]=foo[3]
		foo[2]=foo[4]
	foo[3]=a
	foo[4]=b
	m=re.search("pfamA_id.*?; ",foo[8])
	n=re.sub(" pfamA_id.*?;","",foo[8])
	foo[8]=m.group(0)+n
	if("ccds_id" in foo[8]):
		m=re.search("ccds_id.*?;",foo[8])
		n=re.sub("ccds_id.*?; ","",foo[8])
		foo[8]=n+m.group(0)
	foo=[foo[0],foo[1],foo[2],str(int(foo[2])-int(foo[1])),foo[3],foo[4],foo[5],foo[6],foo[7],foo[8]]
	line="\t".join(foo)
	sys.stdout.write(line+"\n")