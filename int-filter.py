import string
import sys
import re
from optparse import OptionParser

def judge_impact(impacts,
                dn = set(['stop_gained','stop_lost','start_lost','initiator_codon_variant','rare_amino_acid_variant','missense_variant','5_prime_UTR_premature_start_codon_gain_variant','protein_altering_variant']),
                ds = set(['stop_retained_variant','synonymous_variant','exon_variant','start_retained_variant','incomplete_terminal_codon_variant']),
                na = set(['splice_acceptor_variant','splice_donor_variant','frameshift_variant','exon_loss_variant','chromosomal_deletion','inframe_insertion','inframe_deletion','disruptive_inframe_deletion','disruptive_inframe_insertion','5_prime_UTR_truncation','3_prime_UTR_truncation','splice_region_variant','mature_miRNA_variant','regulatory_region_variant','TF_binding_site_variant','regulatory_region_ablation','regulatory_region_amplification','TFBS_ablation','TFBS_amplification','5_prime_UTR_variant','3_prime_UTR_variant','intron_variant','upstream_gene_variant','downstream_gene_variant','intergenic_variant','conserved_intergenic_variant','intragenic_variant','gene_variant','transcript_variant','conserved_intron_variant','nc_transcript_variant','NMD_transcript_variant','non_coding_exon_variant','transcript_ablation','transcript_amplification','feature_elongation','feature_truncation','coding_sequence_variant','non_coding_transcript_exon_variant','non_coding_transcript_variant'])):
    #coding_sequence_variant is really like a multi-nucleotide truncation - hard to say what it does):
    types=[]
    for z in impacts.split('&'):
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

keys = "Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM".split(",")

def concat(r):
    info = dict(tuple(x.split("=")) if "=" in x else (x, None) for x in r.info.split(";") )
    r.maf = "|".join(info["AF"].lstrip("AF=").split(","))

    csqs = [x.split("|") for x in info['CSQ'].split(",")]

    alld = {}
    for csq in csqs:

        d = dict(zip(keys, csq))
        if d['Feature'] != r.transid: continue

        d['impact'] = judge_impact(d['Consequence'])
        #d['csq'] = csq

        if alld == {}:
            for k in d:
                if d[k] == "": d[k] = "."
            alld = d
        else:
            for k in alld:
                alld[k] = alld[k] + "|" + (d[k] or ".")

    return alld

def split_fields(fields, field_num):
    return "|".join([fields[k][field_num] for k in fields])

from collections import OrderedDict, defaultdict

parser = OptionParser()
parser.add_option("-n","--nodom",dest="bool",help="default,nodoms",action="store_true")
parser.add_option("-d","--domain",dest="bool",help="domains",action="store_false")
(options,args)=parser.parse_args()
if options.bool==True:
    class NoDomRecord(object):
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

    for r in (NoDomRecord(l.rstrip().split("\t")) for l in sys.stdin):
        alld = concat(r)
        csq = "\t".join(alld[k] for k in keys)

        print "\t".join([r.chr, r.start,
            r.end,r.ref,r.alt,".",r.uniqid,r.uniqid,r.len,r.covratio,r.gene,r.maf,
            alld['impact'], csq])




if options.bool==False:
    class Record(object):
        def __init__(self, fields):
            self.chr = fields[0]
            self.start = fields[1]
            self.end = fields[2]
            self.domain = fields[3]
            self.transid = fields[4]
            self.autoreg = fields[5]
            self.gene = fields[6]
            self.uniqid = fields[7]
            self.len = fields[8]
            self.covratio = fields[9]
            self.ref = fields[10]
            self.alt = fields[11]
            self.info = fields[12]

    for line in sys.stdin:
        fields=line.rstrip().split("\t")
        r_=Record(fields)

        alld = concat(r_)
        csq = "\t".join(alld[k] for k in keys)
        print "\t".join([r_.chr,r_.start,r_.end,r_.ref,r_.alt,r_.domain,r_.autoreg,r_.uniqid,r_.len,r_.covratio,r_.gene,r_.maf,alld['impact'],csq])
