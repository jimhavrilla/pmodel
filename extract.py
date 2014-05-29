import vcf
import sys
import string
reclist=[];
ESPlist = ["ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr2.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr3.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr4.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr5.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr6.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr7.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr8.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr9.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr10.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr11.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr12.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr13.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr14.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr15.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr16.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr17.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr18.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr19.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr20.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr21.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chr22.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chrX.snps_indels.vcf",
"ESP6500SI-V2-SSA137.updatedRsIds.chrY.snps_indels.vcf"]

for i in range(len(ESPlist)):
  vcf_reader = vcf.Reader(open(ESPlist[i]))
  for record in vcf_reader:
    minilist=[record.CHROM,record.INFO['GL'][0],str(record.ID),str(int(record.POS)-1),str(record.POS),str(record.REF),str(record.ALT),str(round(float(record.INFO['MAF'][2])/100,6)),str(round(float(record.INFO['MAF'][1])/100,6)),str(round(float(record.INFO['MAF'][0])/100,6))]
    if minilist[1] is None:
      continue
    for i in range(len(record.INFO['FG'])):
      minilist.append(record.INFO['FG'][i])
    reclist.append(minilist)

l=[l for l in reclist if l[1]=='TP53']

i=0;
for a in l:
	# l[i][0]="chr"+l[i][0]
	# b=l[i][1]
	# c=l[i][2]
	# l[i][1]=l[i][3]
	# l[i][2]=l[i][4]
	# l[i][3]=b
	# l[i][4]=c
	# i=i+1
	d=[]
	d.append("chr"+a[0])
	b=a[1]
	c=a[2]
	d.append(a[3])
	d.append(a[4])
	d.append(b)
	d.append(c)
	l[i]=d
	i=i+1

while True:
	try:
		with open('./pmodel/gene.bed','w+') as f:
			for y in l:
				print>>f,'\t'.join([x for x in y])
				i=i+1
	except TypeError:
		print l[i] 
		break
	break