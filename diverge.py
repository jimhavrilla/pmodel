import string
import sys
import numpy
import collections
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--uniqids","-u",help="run the uniqids protocol",action="store_true")
parser.add_argument("--pairs","-p",help="run the pairs protocol",action="store_true")
parser.add_argument("--domains","-d",help="run the domain protocol",action="store_true")

args=parser.parse_args()

def mad(a, axis=None):

    med = numpy.median(a, axis=axis)                # Median along given axis
    if axis is None:
        umed = med                              # med is a scalar
    else:
        umed = numpy.expand_dims(med, axis)         # Bring back the vanished axis
    mad = numpy.median(numpy.absolute(a - umed), axis=axis)  # MAD along given axis

    return mad

if args.uniqids:
	class Record(object):
		def __init__(self, fields):
			self.domain = fields[0]
			self.gene = fields[1]
			self.autoreg = fields[2]
			self.uniqid = fields[3]
			self.coverage = int(fields[4])
			self.truelength = int(fields[5])
			self.covratio = float(fields[6])
			self.nct = int(fields[7])
			self.sct = int(fields[8])
			self.ct = int(fields[9])
			self.dnds = float(fields[10])
			self.z = float(fields[11])
			self.density = float(fields[12])
			self.maf = float(fields[13])
			self.domcount = int(fields[14])

	gene={}

	#gene, domain list, count of domains, bps, bp sum, domcounts, domcount sum, dnds, min/max/stdev dnds, z-score, rank metric
	for line in sys.stdin:
		fields=line.rstrip().split("\t")
		r_=Record(fields)
		x=r_.gene
		try:
			foo=gene[x]
			foo[0]=x
			foo[1][r_.autoreg]=r_.domain
			foo[2]=foo[2]+1
			foo[3].append(r_.truelength)
			foo[4].append(r_.domcount)
			foo[5].append(r_.dnds)
			foo[6]=min(foo[5])
			foo[7]=max(foo[5])
			foo[8]=mad(foo[5])
			foo[9].append(r_.z)
			foo[10]=max(foo[9])-min(foo[9]) #range z-score
			foo[11]=(foo[7]-foo[6])*foo[8]*foo[10] #(max-min)*mad*z-range
			foo[12].append(r_.coverage)
			foo[13].append(r_.covratio)
			foo[14].append(r_.density)
			gene[x]=foo
		except KeyError:
			gene[x]=[x, collections.OrderedDict([(r_.autoreg,r_.domain)]),1,[r_.truelength],[r_.domcount],[r_.dnds],0,0,0,[r_.z],0,0,[r_.coverage],[r_.covratio],[r_.density]]

	for x in gene:
		i=0
		for y in gene[x][1].keys():
			sys.stdout.write('\t'.join([gene[x][0],y,gene[x][1][y],str(gene[x][2]),str(gene[x][3][i]),str(gene[x][4][i]),str(gene[x][5][i]),str(gene[x][6]),str(gene[x][7]),str(gene[x][8]),str(gene[x][9][i]),str(gene[x][10]),str(gene[x][11]),str(gene[x][12][i]),str(gene[x][13][i]),str(gene[x][14][i])])+'\n')
			i=i+1
		# sys.stdout.write('\t'.join([gene[x][0],','.join(gene[x][1]),str(gene[x][2]),str(gene[x][4]),','.join([str(y) for y in gene[x][5]]),','.join([str(y) for y in gene[x][7]]),str(gene[x][8]),str(gene[x][9]),str(gene[x][10]),','.join([str(y) for y in gene[x][11]]),str(gene[x][12]),str(gene[x][13]),','.join([str(y) for y in gene[x][14]]),str(gene[x][15]),','.join([str(y) for y in gene[x][16]]),','.join([str(y) for y in gene[x][17]])])+'\n')
		# gene, domains, num_of_domains, bp_total, domcount_total, dnds_list, min_dnds, max_dnds, mad_dnds, z-scores, z-score_range, divergence metric

if args.domains:
	class Record1(object):
		def __init__(self, fields):
			self.domain = fields[0]
			self.gene = fields[1]
			self.nct = fields[2]
			self.sct = fields[3]
			self.ct = fields[4]
			self.domcount = fields[5]
			self.bp = fields[6]
			self.mmaf = fields[7]
			self.dnds = fields[8]
			self.vbpratio = fields[9]

	gene={}

	#domain list, count of domains, bps, bp sum, domcounts, domcount sum, dnds, min/max/stdev dnds, rank metric
	for line in sys.stdin:
		fields=line.rstrip().split("\t")
		r_=Record1(fields)
		for x in r_.gene.split(","):
			try:
				foo=gene[x]
				foo[0]=x
				foo[1][r_.domain]=None
				foo[2]=foo[2]+1
				foo[3].append(int(r_.bp))
				foo[4]=sum(foo[3])
				foo[5].append(int(r_.domcount))
				foo[6]=sum(foo[5])
				foo[7].append(float(r_.dnds))
				foo[8]=min(foo[7])
				foo[9]=max(foo[7])
				foo[10]=mad(foo[7])
				foo[11]=(foo[9]-foo[8])*foo[10] #(max-min)*mad
				gene[x]=foo
			except KeyError:
				gene[x]=[x,collections.OrderedDict.fromkeys([r_.domain]),1,[int(r_.bp)],0,[int(r_.domcount)],0,[float(r_.dnds)],0,0,0,0]

	for x in gene:
		sys.stdout.write('\t'.join([gene[x][0],','.join(gene[x][1]),str(gene[x][2]),str(gene[x][4]),str(gene[x][6]),','.join([str(y) for y in gene[x][7]]),str(gene[x][8]),str(gene[x][9]),str(gene[x][10]),str(gene[x][11])])+'\n')

if args.pairs:
	class Record1(object):
		def __init__(self, fields):
			self.domain = fields[0]
			self.gene = fields[1]
			self.ct = fields[2]
			self.nct = fields[3]
			self.sct = fields[4]
			self.dnds = float(fields[5])
			self.z = float(fields[6])
			self.mmaf = float(fields[7])
			self.domcount = int(fields[8])
			self.bp = int(fields[9])

	gene={}

	#gene, domain list, count of domains, bps, bp sum, domcounts, domcount sum, dnds, min/max/stdev dnds, z-score, rank metric
	for line in sys.stdin:
		fields=line.rstrip().split("\t")
		r_=Record1(fields)
		x=r_.gene
		try:
			foo=gene[x]
			foo[0]=x
			foo[1][r_.domain]=None
			foo[2]=foo[2]+1
			foo[3].append(r_.bp)
			foo[4]=sum(foo[3])
			foo[5].append(r_.domcount)
			foo[6]=sum(foo[5])
			foo[7].append(r_.dnds)
			foo[8]=min(foo[7])
			foo[9]=max(foo[7])
			foo[10]=mad(foo[7])
			foo[11].append(r_.z)
			foo[12]=max(foo[11])-min(foo[11]) #range z-score
			foo[13]=(foo[9]-foo[8])*foo[10]*foo[12] #(max-min)*mad
			gene[x]=foo
		except KeyError:
			gene[x]=[x, collections.OrderedDict.fromkeys([r_.domain]),1,[r_.bp],0,[r_.domcount],0,[r_.dnds],0,0,0,[r_.z],0,0]

	for x in gene:
		sys.stdout.write('\t'.join([gene[x][0],','.join(gene[x][1]),str(gene[x][2]),str(gene[x][4]),str(gene[x][6]),','.join([str(y) for y in gene[x][7]]),str(gene[x][8]),str(gene[x][9]),str(gene[x][10]),','.join([str(y) for y in gene[x][11]]),str(gene[x][12]),str(gene[x][13])])+'\n')
		# gene, domains, num_of_domains, bp_total, domcount_total, dnds_list, min_dnds, max_dnds, mad_dnds, z-scores, z-score_range, divergence metric
		