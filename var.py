import string
import sys
import re
from optparse import OptionParser

def judge_impact(impact):
	types=[]
	dn = set(['stop_gained','stop_lost','start_lost','initiator_codon_variant','rare_amino_acid_variant','missense_variant','5_prime_UTR_premature_start_codon_gain_variant','protein_altering_variant'])
	ds = set(['stop_retained_variant','synonymous_variant','exon_variant','start_retained_variant','incomplete_terminal_codon_variant'])
	na = set(['splice_acceptor_variant','splice_donor_variant','frameshift_variant','exon_loss_variant','chromosomal_deletion','inframe_insertion','inframe_deletion','disruptive_inframe_deletion','disruptive_inframe_insertion','5_prime_UTR_truncation','3_prime_UTR_truncation','splice_region_variant','mature_miRNA_variant','regulatory_region_variant','TF_binding_site_variant','regulatory_region_ablation','regulatory_region_amplification','TFBS_ablation','TFBS_amplification','5_prime_UTR_variant','3_prime_UTR_variant','intron_variant','upstream_gene_variant','downstream_gene_variant','intergenic_variant','conserved_intergenic_variant','intragenic_variant','gene_variant','transcript_variant','conserved_intron_variant','nc_transcript_variant','NMD_transcript_variant','non_coding_exon_variant','transcript_ablation','transcript_amplification','feature_elongation','feature_truncation','coding_sequence_variant','non_coding_transcript_exon_variant','non_coding_transcript_variant']) #coding_sequence_variant is really like a multi-nucleotide truncation - hard to say what it does
	foo=impact.split('|')
	for x in foo:
		y=x.split('&')
		for z in y:
			if z in dn:
				types.append("dn")
				break
			if z in ds:
				types.append("ds")
				break
			if z in na:
				types.append("na")
				break
	
	return '|'.join(types)

def split_fields(fields,field_num):
	return "|".join([fields[k][field_num] for k in csqstring[j]])

parser = OptionParser()
parser.add_option("-n","--nodom",dest="bool",help="default,nodoms",action="store_true")
parser.add_option("-d","--domain",dest="bool",help="domains",action="store_false")
(options,args)=parser.parse_args()
if options.bool==True:
	class Record(object):
		def __init__(self, fields):
			self.chr = fields[0]
			self.start = fields[1]
			self.end = fields[2]
			self.transid = fields[3]
			self.gene = fields[4]
			self.uniqid = fields[5] 
			self.len = fields[6]
			self.covratio = fields[7]
			self.ref = fields[8]
			self.alt = fields[9]
			self.info = fields[10]

	for line in sys.stdin:
		fields=line.rstrip().split("\t")
		r_=Record(fields)
		foo=r_.info.split(";")
		if r_.chr=="X" or r_.chr=="Y": # chr X and Y have AC_Hemi, a field that other chromosomes do NOT have
			r_.maf="|".join(foo[12].lstrip("AF=").split(","))
		else:
			r_.maf="|".join(foo[11].lstrip("AF=").split(","))
		CSQ=foo[len(foo)-1]
		toks=[x.split("|") for x in CSQ.lstrip("CSQ=").split(",")]
		var={}
		csqstring={}
		for i in toks:
			#Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM
			trans=i[5]
			impact=i[0]
			allele=i[11]
			cs=i[1:5]
			cs.extend(i[6:11])
			cs=[x if x!="" else "." for x in cs]
			try:
				var[trans][allele]=impact
				csqstring[trans][allele]=cs
			except KeyError:
				var[trans]={allele:impact}
				csqstring[trans]={allele:cs}
		for j in var:
			if r_.transid == j:
				impacts="|".join(var[j].values())
				types=judge_impact(impacts)
				csq1="\t".join([split_fields(csqstring[j],i) for i in range(0,4)])
				csq2="\t".join([split_fields(csqstring[j],i) for i in range(4,9)])
				print "\t".join([r_.chr,str(int(r_.start)-1),r_.end,r_.ref,r_.alt,".",r_.uniqid,r_.uniqid,r_.len,r_.covratio,r_.gene,r_.maf,impacts,types,csq1,j,csq2])

if options.bool==False:
	class Record(object):
		def __init__(self, fields):
			self.chr = fields[0]
			self.start = fields[1]
			self.end = fields[2]
			self.domain = fields[3]
			self.transid = fields[4]
			self.gene = fields[5]
			self.autoreg = fields[6]
			self.uniqid = fields[7]
			self.covratio = fields[8]
			self.len = fields[9]
			self.ref = fields[10]
			self.alt = fields[11]
			self.info = fields[12]

	for line in sys.stdin:
		fields=line.rstrip().split("\t")
		r_=Record(fields)
		foo=r_.info.split(";")
		if r_.chr=="X" or r_.chr=="Y": # chr X and Y have AC_Hemi, a field that other chromosomes do NOT have
			r_.maf="|".join(foo[12].lstrip("AF=").split(","))
		else:
			r_.maf="|".join(foo[11].lstrip("AF=").split(","))
		CSQ=foo[len(foo)-1]
		toks=[x.split("|") for x in CSQ.lstrip("CSQ=").split(",")]
		var={}
		csqstring={}
		for i in toks:
			#Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM
			trans=i[5]
			impact=i[0]
			allele=i[11]
			cs=i[1:5]
			cs.extend(i[6:11])
			cs=[x if x!="" else "." for x in cs]
			try:
				var[trans][allele]=impact
				csqstring[trans][allele]=cs
			except KeyError:
				var[trans]={allele:impact}
				csqstring[trans]={allele:cs}
		for j in var:
			if r_.transid == j:
				impacts="|".join(var[j].values())
				types=judge_impact(impacts)
				csq1="\t".join([split_fields(csqstring[j],i) for i in range(0,4)])
				csq2="\t".join([split_fields(csqstring[j],i) for i in range(4,9)])
				print "\t".join([r_.chr,str(int(r_.start)-1),r_.end,r_.ref,r_.alt,r_.domain,r_.autoreg,r_.uniqid,r_.len,r_.covratio,r_.gene,r_.maf,impacts,types,csq1,j,csq2])
