import sys
from cyvcf2 import VCF, Writer
vcf = VCF(sys.argv[1])

w = Writer("-", vcf)

for v in vcf:
    try:
        if any(float(m) < 10 for m in v.INFO["MAF"].split(",")): continue
    except:
        continue
    w.write_record(v)


