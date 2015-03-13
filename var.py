import string
import sys
import re
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n","--nodom",dest="bool",help="default,nodoms",action="store_true")
parser.add_option("-d","--domain",dest="bool",help="domains",action="store_false")
(options,args)=parser.parse_args()
if options.bool==True:
	class Record1(object):
		def __init__(self, fields):
			self.chr = fields[0]
			self.start = fields[1]
			self.end = fields[2]
			self.gene = fields[3]
			self.uniqid = fields[4] 
			self.covct = fields[5]
			self.len = fields[6]
			self.covratio = fields[7]
			self.ref = fields[8]
			self.alt = fields[9]
			self.info = fields[10]

	for line in sys.stdin:
		fields=line.rstrip().split("\t")
		r_=Record1(fields)
		foo=r_.info.split(";")
		r_.maf=foo[11].lstrip("AF=").split(",")
		alt=r_.alt.split(",")
		ct=0
		for y in r_.maf:
			z=alt[ct]
			#Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE
			CSQ=foo[len(foo)-1]
			for impact in CSQ.lstrip("CSQ=").split(","):
				toks = impact.split("|")
				r_.csq="\t".join([x if x!="" else "." for x in toks]) # adds nodom domain field
				print "\t".join([r_.chr,str(int(r_.start)-1),r_.end,r_.ref,z,".",r_.uniqid,r_.covct,r_.len,r_.covratio,r_.gene,y,r_.csq])
			ct=ct+1

if options.bool==False:
	class Record1(object):
		def __init__(self, fields):
			self.chr = fields[0]
			self.start = fields[1]
			self.end = fields[2]
			self.domain = fields[3]
			self.gene = fields[4]
			self.uniqid = fields[5]
			self.covct = fields[6]
			self.len = fields[7]
			self.covratio = fields[8]
			self.ref = fields[9]
			self.alt = fields[10]
			self.info = fields[11]

	for line in sys.stdin:
		fields=line.rstrip().split("\t")
		r_=Record1(fields)
		foo=r_.info.split(";")
		r_.maf=foo[11].lstrip("AF=").split(",")
		alt=r_.alt.split(",")
		ct=0
		for y in r_.maf:
			z=alt[ct]
			#Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE
			CSQ=foo[len(foo)-1]
			for impact in CSQ.lstrip("CSQ=").split(","):
				toks = impact.split("|")
				r_.csq="\t".join([x if x!="" else "." for x in toks])
				print "\t".join([r_.chr,str(int(r_.start)-1),r_.end,r_.ref,z,r_.domain,r_.uniqid,r_.covct,r_.len,r_.covratio,r_.gene,y,r_.csq])
			ct=ct+1