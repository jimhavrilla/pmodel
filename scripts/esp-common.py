from cyvcf2 import VCF, Writer
vcf = VCF("../gemini_install/data/gemini_data/ESP6500SI.all.snps_indels.tidy.v2.vcf.gz")

w = Writer("-", vcf)

for v in vcf:
    try:
        if any(float(m) < 5 for m in v.INFO["MAF"].split(",")): continue
    except:
        continue
    w.write_record(v)


