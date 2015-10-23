#!/usr/bin/env python
import sys
from toolshed import reader, header, nopen
from optparse import OptionParser

#parser = OptionParser()
#
#parser.add_option( "--headers",
#    dest="headers",
#    help="header name CSV")
#
#(options, args) = parser.parse_args()
#
#if not options.headers:
#    parser.error('Headers not given')
#
#headers = options.headers.split(',')
#
#if len(headers) <= 1:
#    sys.stderr.write("Must give more than one header name\n");
#    exit(1)
#
#s = set()

keys = "Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM".split(",")

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


def concat(r):
    info = dict(tuple(x.split("=")) if "=" in x else (x, None) for x in r['INFO'].split(";") )
    r['maf'] = "|".join(info["AF"].lstrip("AF=").split(","))

    csqs = [x.split("|") for x in info['CSQ'].split(",")]

    alld = {}
    for csq in csqs:

        d = dict(zip(keys, csq))
        if d['Feature'] != r['transcript_id']: continue

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


print '#' + '\t'.join(['chr', \
                       'start', \
                       'end', \
                       'POS', \
                       'REF', \
                       'ALT', \
                       'transcript_id', \
                       'pfamA_id', \
                       'pfamA_auto_reg', \
                       'uniq_id', \
                       'len', \
                       'ratio_of_cover_to_total', \
                       'gene_name', \
                       'maf', \
                       'impact'] + \
                        keys)

for d in reader('-',header='ordered') :
    if 'pfamA_id' not in d:
        d['pfamA_id'] = '.'
    if 'pfamA_auto_reg' not in d:
        d['pfamA_auto_reg'] = d['uniq_id']

    alld = concat(d)
    csq = "\t".join(alld[k] for k in keys)
    print "\t".join([d['chr'],
                     d['start'],
                     d['end'],
                     d['POS'],
                     d['REF'],
                     d['ALT'],
                     d['transcript_id'],
                     d['pfamA_id'],
                     d['pfamA_auto_reg'],
                     d['uniq_id'],
                     d['len'],
                     d['ratio_of_cover_to_total'],
                     d['gene_name'],
                     d['maf'],
                     alld['impact'],
                     csq])
