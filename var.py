import string
import fileinput
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
			self.uniqid = fields[3].rstrip(";")
			self.gene = fields[4].rstrip(";").strip("\"")
			self.ref = fields[5]
			self.alt = fields[6]
			self.info = fields[7]

	for line in sys.stdin:
		fields=line.rstrip().split("\t")
		r_=Record1(fields)
		foo=r_.info.split(";")
		r_.maf=foo[4].lstrip("MAF=").split(",")
		r_.eamaf=str(float(r_.maf[0])/100);r_.aamaf=str(float(r_.maf[1])/100);r_.maf=str(float(r_.maf[2])/100)
		#Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE
		CSQ=foo[len(foo)-1]
		for impact in CSQ.lstrip("CSQ=").split(","):
			toks = impact.split("|")
			r_.csq="\t".join([x if x!="" else "." for x in toks])
			print "\t".join([r_.chr,str(int(r_.start)-1),r_.end,r_.ref,r_.alt,".",r_.uniqid,r_.gene,r_.eamaf,r_.aamaf,r_.maf,r_.csq])
if options.bool==False:
	class Record1(object):
		def __init__(self, fields):
			self.chr = fields[0]
			self.start = fields[1]
			self.end = fields[2]
			self.domain = fields[3]
			self.gene = fields[4].rstrip(";").strip("\"")
			self.uniqid = fields[5].rstrip(";")
			self.ref = fields[6]
			self.alt = fields[7]
			self.info = fields[8]

	for line in sys.stdin:
		fields=line.rstrip().split("\t")
		r_=Record1(fields)
		foo=r_.info.split(";")
		r_.maf=foo[4].lstrip("MAF=").split(",")
		r_.eamaf=str(float(r_.maf[0])/100);r_.aamaf=str(float(r_.maf[1])/100);r_.maf=str(float(r_.maf[2])/100)
		#Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE
		CSQ=foo[len(foo)-1]
		for impact in CSQ.lstrip("CSQ=").split(","):
			toks = impact.split("|")
			r_.csq="\t".join([x if x!="" else "." for x in toks])
			print "\t".join([r_.chr,str(int(r_.start)-1),r_.end,r_.ref,r_.alt,r_.domain,r_.uniqid,r_.gene,r_.eamaf,r_.aamaf,r_.maf,r_.csq])