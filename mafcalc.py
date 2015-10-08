#mafcalc.py - changes table according to maf cutoff

import string
import sys
from csv import reader
from gzip import open
import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-m","--maf",help="default>=0 MAF",default=0,type=float)
parser.add_argument("-f","--files",help="the files to input/output",nargs='*')
args = parser.parse_args()


class Record(object):
	def __init__(self,fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.transid = fields[3]
		self.domain = fields[4]
		self.gene = fields[5]
		self.autoreg = fields[6]
		self.covratio = fields[7]
		self.length = fields[8]
		self.dn = fields[9]
		self.ds = fields[10]
		self.na = fields[11]
		self.dnds = fields[12]
		self.density = fields[13]
		self.fvrv = fields[14]
		self.prevalence = fields[15]
		self.mafs = fields[16]
		self.impacts = fields[17]
		self.types = fields[18]
		self.starts = fields[19]
		self.ends = fields[20]
		self.mod = fields[21]

with open(args.files[0],"r") as f:
	bed = reader(f,delimiter="\t")
	bed.next() #the header
	for record in bed:
		r_=Record(record)
		try:
			mafs=np.array([y for x in r_.mafs.split(",") for y in x.split("|")],dtype=float)
		except ValueError:
			sys.stdout.write("\t".join([r_.chr,r_.start,r_.end,r_.transid,r_.domain,r_.gene,r_.autoreg,r_.length,str(r_.dn),str(r_.ds),r_.dnds,r_.density,r_.fvrv,r_.mafs,r_.impacts,r_.types,r_.starts,r_.ends,r_.mod])+"\n")
			continue
		types=np.array([y for x in r_.types.split(",") for y in x.split("|")],dtype=str)
		impacts=np.array([y for x in r_.impacts.split(",") for y in x.split("|")],dtype=str)
		impacts=impacts[np.where(mafs>args.maf)[0]]
		types=types[np.where(mafs>args.maf)[0]]
		mafs=mafs[np.where(mafs>args.maf)[0]]
		r_.dn=float(len(types[np.where(types=='dn')[0]]))
		r_.ds=float(len(types[np.where(types=='ds')[0]]))
		if r_.ds!=0:
			r_.dnds=str(r_.dn/r_.ds)
		else:
			r_.dnds=str(r_.dn/1)
		r_.density=str(float(len(types))/float(r_.length))
		r_.mod=str(args.maf)
		sys.stdout.write("\t".join([r_.chr,r_.start,r_.end,r_.transid,r_.domain,r_.gene,r_.autoreg,r_.length,str(r_.dn),str(r_.ds),r_.dnds,r_.density,r_.fvrv,r_.mafs,r_.impacts,r_.types,r_.starts,r_.ends,r_.mod])+"\n")
