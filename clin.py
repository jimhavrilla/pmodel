from cyvcf2 import VCF, Writer
import re
import sys

patt = re.compile(',|\|')

def clinvar(v):
    return v.INFO.get("CLNSIG") == "5"
    #return [x in "45" for x in re.split(patt,v.INFO.get("CLNSIG"))][0]

def aaf(v, max_aaf):
    if v.INFO.get("max_aaf_all") != None:
        return float(v.INFO.get("max_aaf_all")) <= float(max_aaf)
    else:
        return True

vcf_path = sys.argv[1]
max_aaf = float(sys.argv[2])

viter = VCF(vcf_path)
w = Writer("-", viter)
pos = lambda v: (v.CHROM, v.start, v.end)
for v in viter:
    if clinvar(v) and aaf(v, max_aaf):
        w.write_record(v)
