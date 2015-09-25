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
		self.covcount = fields[8]
		self.len = fields[9]
		self.covratio = fields[10]
		self.gene = fields[11]
		self.maf = fields[12]
		self.impact = fields[13]
		self.type = fields[14]
		self.codonchange = fields[15]
		self.aachange = fields[16]
		self.geneidcsq = fields[17]
		self.genecsq = fields[18] # can be different from the gene name from ensembl exon file
		self.transidcsq = fields[19] # only picked those that match our transcript list, but may be a versioning issue if gene names don't match
		self.exonnum = fields[20]
		self.polyphen = fields[21]
		self.sift = fields[22]
		self.position = fields[23]
		self.biotype = fields[24]

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
		ct[r_.autoreg+r_.gene]=[float(0),float(0),float(0),float(r_.len),"\t".join([r_.domain,r_.gene,r_.autoreg,r_.uniqid,r_.covcount,r_.covratio,r_.len])]
	if oldauto!=None and oldauto!=r_.autoreg:
		if ct[oldauto+oldgene][1]==0:
			dnds=str(int(ct[oldauto+oldgene][0]/1))
		else:
			dnds=str(round(ct[oldauto+oldgene][0]/ct[oldauto+oldgene][1],4))
		density=str(round((ct[oldauto+oldgene][0]+ct[oldauto+oldgene][1]+ct[oldauto+oldgene][2])/ct[oldauto+oldgene][3],4))
		print "\t".join([ct[oldauto+oldgene][4],str(int(ct[oldauto+oldgene][0])),str(int(ct[oldauto+oldgene][1])),str(int(ct[oldauto+oldgene][2])),dnds,density])
