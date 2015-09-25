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
		self.length = fields[7]
		self.dn = fields[8]
		self.ds = fields[9]
		self.dnds = fields[10]
		self.density = fields[11]
		self.fvrv = fields[12]
		self.mafs = fields[13]
		self.impacts = fields[14]
		self.types = fields[15]
		self.starts = fields[16]
		self.ends = fields[17]
		self.mod = fields[18]

with open(args.files[0],"r") as f:
	bed = reader(f,delimiter="\t")
	bed.next() #the header
	for record in bed:
		r_=Record(record)
		try:
			mafs=np.array([float(x) for x in r_.mafs.split(",")])
		except ValueError:
			sys.stdout.write("\t".join([r_.chr,r_.start,r_.end,r_.transid,r_.domain,r_.gene,r_.autoreg,r_.length,str(r_.dn),str(r_.ds),r_.dnds,r_.density,r_.fvrv,r_.mafs,r_.impacts,r_.types,r_.starts,r_.ends,r_.mod])+"\n")
			continue
		types=np.array(r_.types.split(","))
		impacts=np.array(r_.impacts.split(","))
		r_.impacts=impacts[np.where(mafs>args.maf)[0]]
		r_.types=types[np.where(mafs>args.maf)[0]]
		r_.mafs=mafs[np.where(mafs>args.maf)[0]]
		r_.dn=float(len(r_.types[np.where(r_.types=='dn')[0]]))
		r_.ds=float(len(r_.types[np.where(r_.types=='ds')[0]]))
		if r_.ds!=0:
			r_.dnds=str(r_.dn/r_.ds)
		else:
			r_.dnds=str(r_.dn/1)
		r_.density=str(float(len(r_.types))/float(r_.length))
		r_.mod=str(args.maf)
		r_.mafs=[str(i) for i in r_.mafs]
		sys.stdout.write("\t".join([r_.chr,r_.start,r_.end,r_.transid,r_.domain,r_.gene,r_.autoreg,r_.length,str(r_.dn),str(r_.ds),r_.dnds,r_.density,r_.fvrv,",".join(r_.mafs),",".join(r_.impacts),",".join(r_.types),r_.starts,r_.ends,r_.mod])+"\n")