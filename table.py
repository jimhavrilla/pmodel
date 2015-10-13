import sys
import string

class Record(object):
	def __init__(self,fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.ref = fields[3]
		self.alt = fields[4]
		self.domain = fields[5]
		self.autoreg = fields[6] #uniqid for nodoms
		self.uniqid = fields[7]
		self.len = fields[8]
		self.covratio = fields[9]
		self.gene = fields[10]
		self.maf = fields[11]
		self.type = fields[12]
		self.impact = fields[13]
		self.codonchange = fields[14]
		self.aachange = fields[15]
		self.geneidcsq = fields[16]
		self.genecsq = fields[17] # can be different from the gene name from ensembl exon file
		self.transidcsq = fields[18] # only picked those that match our transcript list, but may be a versioning issue if gene names don't match
		self.exonnum = fields[19]
		self.polyphen = fields[20]
		self.sift = fields[21]
		self.position = fields[22]
		self.biotype = fields[23]

ct={} # will be an autoreg/gene key that refers to a list of three values: [0] dn, [1] ds, [2] na
oldauto=None
oldgene=None
for line in sys.stdin:
	try:
		oldauto=r_.autoreg
		oldgene=r_.gene
	except:
		pass
	r_=Record(line.rstrip().split("\t"))
	types=r_.type.split("|")
	if r_.autoreg+r_.gene in ct:
		if "dn" in types:
			ct[r_.autoreg+r_.gene][0]+=1
		if "ds" in types:
			ct[r_.autoreg+r_.gene][1]+=1
		if "na" in types:
			ct[r_.autoreg+r_.gene][2]+=1
	else:
		ct[r_.autoreg+r_.gene]=[float(0),float(0),float(0),float(r_.len),"\t".join([r_.domain,r_.gene,r_.autoreg,r_.uniqid,r_.covratio,r_.len])]
		if "dn" in types:
			ct[r_.autoreg+r_.gene][0]+=1
		if "ds" in types:
			ct[r_.autoreg+r_.gene][1]+=1
		if "na" in types:
			ct[r_.autoreg+r_.gene][2]+=1
	if oldauto!=None and oldauto!=r_.autoreg:
		if ct[oldauto+oldgene][1]==0:
			dnds=str(int(ct[oldauto+oldgene][0]/1))
		else:
			dnds=str(round(ct[oldauto+oldgene][0]/ct[oldauto+oldgene][1],4))
		density=str(round((ct[oldauto+oldgene][0]+ct[oldauto+oldgene][1]+ct[oldauto+oldgene][2])/ct[oldauto+oldgene][3],4))
		print "\t".join([ct[oldauto+oldgene][4],str(int(ct[oldauto+oldgene][0])),str(int(ct[oldauto+oldgene][1])),str(int(ct[oldauto+oldgene][2])),dnds,density])
