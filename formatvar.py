import string
import sys
import re
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--maf","-m",help="default=>0 maf",default=0,type=float)
args=parser.parse_args()

mpar=args.maf

class Record1(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.ref = fields[3]
		self.alt = fields[4]
		self.domain = fields[5]
		self.uniqid = fields[6]
		self.gene = fields[7]
		self.maf = fields[8]
		self.impact = fields[9]
		self.info = fields[10:20]

class Csq(object):
	def __init__(self):
		self.impact=''
		self.info=''

def judge_impact(impact_list):
	impact={}
	for c in impact_list:	
		dn = set(['stop_gained','stop_lost','start_lost','initiator_codon_variant','rare_amino_acid_variant','missense_variant','5_prime_UTR_premature_start_codon_gain_variant'])
		ds = set(['stop_retained_variant','synonymous_variant','exon_variant','start_retained_variant','incomplete_terminal_codon_variant'])
		na = set(['splice_acceptor_variant','splice_donor_variant','frameshift_variant','exon_loss_variant','chromosomal_deletion','inframe_insertion','inframe_deletion','disruptive_inframe_deletion','disruptive_inframe_insertion','5_prime_UTR_truncation','3_prime_UTR_truncation','splice_region_variant','mature_miRNA_variant','regulatory_region_variant','TF_binding_site_variant','regulatory_region_ablation','regulatory_region_amplification','TFBS_ablation','TFBS_amplification','5_prime_UTR_variant','3_prime_UTR_variant','intron_variant','upstream_gene_variant','downstream_gene_variant','intergenic_variant','conserved_intergenic_variant','intragenic_variant','gene_variant','transcript_variant','conserved_intron_variant','nc_transcript_variant','NMD_transcript_variant','non_coding_exon_variant','transcript_ablation','transcript_amplification','feature_elongation','feature_truncation','coding_sequence_variant']) #coding_sequence_variant is really like a multi-nucleotide truncation - hard to say what it does
		foo=c.impact.split('&')
		for x in foo:
			if x in dn:
				impact['dn']=x+'\t'+'dn'+'\t'+'\t'.join(c.info)
			if x in ds:
				impact['ds']=x+'\t'+'ds'+'\t'+'\t'.join(c.info)
			if x in na:
				impact['na']=x+'\t'+'na'+'\t'+'\t'.join(c.info)

	return impact

old_r=None
old_c=None
impact_list=[]
for line in sys.stdin:
	fields=line.rstrip().split("|")
	r_=Record1(fields)
	if float(r_.maf) < mpar:
		continue
	c_=Csq();c_.impact=r_.impact;c_.info=r_.info
	if old_r!=None:
		impact_list.append(old_c)
	if old_r!=None and old_r.start+old_r.end != r_.start+r_.end:
		d_=judge_impact(impact_list)
		for x in d_.values():
			print "\t".join([old_r.chr,old_r.start,old_r.end,old_r.ref,old_r.alt,old_r.domain,old_r.uniqid,old_r.gene,old_r.maf])+'\t'+x
		impact_list=[]
	old_r=r_
	old_c=c_

if old_c!=None:
	impact_list.append(old_c)
	d=judge_impact(impact_list)
	for x in d_.values():
		print "\t".join([old_r.chr,old_r.start,old_r.end,old_r.ref,old_r.alt,old_r.domain,old_r.uniqid,old_r.gene,old_r.maf])+'\t'+x