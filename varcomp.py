import string
import sys
import cyvcf2
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-c","--clinvar",dest="bool",help="default,clinvar",action="store_true")
parser.add_argument("-g","--1000genomes",dest="bool",help="",action="store_false")
parser.add_argument("-f","--files",help="the files to input/output",nargs='*')
args = parser.parse_args()
if args.bool == True:
	rdr = cyvcf2.VCF(args.files[0])
	wtr = cyvcf2.Writer("/dev/stdout", rdr) # "/dev/stdout" pipes to stdout
	for variant in rdr:
		clnsig = variant.INFO.get('CLNSIG').split('|')
		clnsrc = variant.INFO.get('CLNSRC')
		if ("5" in clnsig) and clnsrc != ".":
			wtr.write_record(variant)

if args.bool == False:
	transcripts=set()
	with open(args.files[1],"r") as f1:
		for line in f1:
			fields=line.split()
			transcripts.add(fields[1])
	rdr = cyvcf2.VCF(args.files[0])
	nuc=set(['A','C','G','T'])
	for variant in rdr:
		if not variant.FILTER and variant.REF in nuc and [i in variant.ALT for i in nuc][0]:
			csq = variant.INFO.get('CSQ').split(",")
			if variant.ID == "":
				variant.ID = "."
			# variant.FILTER="PASS" can't display NOTHING
			for impact in csq:
				toks=impact.split("|")
				trans=toks[2]
				if trans in transcripts:	
					gene=toks[1]
					consequence=toks[4]
					sys.stdout.write("\t".join([str(variant.CHROM),str(variant.POS),str(variant.ID),variant.REF,",".join(variant.ALT),str(variant.QUAL),"PASS",str(variant.INFO.get('AC')),str(variant.INFO.get('AF')),str(variant.INFO.get('AN')),str(variant.INFO.get('NS')),gene,trans,consequence])+"\n")