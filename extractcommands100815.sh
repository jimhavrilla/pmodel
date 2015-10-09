#pipeline

export DATA=~/work/data/pmodeldata
export SOFTWARE=~/software
export SCRATCH=~/work/scratch/analysis

#bill's domain count file:

mysql -N --raw -h wrpxdb.its.virginia.edu -u web_user -pfasta_www pfam27 < count_human_pfam.sql > human_pfam.counts

#for vcfs:

perl $SOFTWARE/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/ExAC.r0.3.sites.vep.vcf.gz --cache --merged --sift b --polyphen b --symbol --numbers --biotype --total_length --allele_number -o $DATA/VEPEXAC3.vcf --vcf --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM --force_overwrite --port 3337 --dir ~/.vep/ --compress "gunzip -c" --offline
zcat $DATA/Homo_sapiens.GRCh37.75.gtf.gz | grep protein_coding$'\t'exon | perl -pe 's/protein_coding\texon\t//g' | grep -P '^\d*[0-9X-Y]\t' | perl gtf2bed.pl | sort -k10,10 > $DATA/exons.bed
zcat $DATA/Homo_sapiens.GRCh37.75.gtf.gz | grep protein_coding$'\t'UTR | perl -pe 's/protein_coding\tUTR\t//g' | grep -P '^\d*[0-9X-Y]\t' | perl gtf2bed.pl | sort -k10,10 > $DATA/utrs.bed
python filterutrs.py $DATA/exons.bed $DATA/utrs.bed > $DATA/GRCh37.bed

#for vcf coverage:
for chrom in {1..22} X Y ; do wget -P $DATA ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/coverage/Panel.chr$chrom.coverage.txt.gz; done
gzcat $DATA/Panel.chr*.coverage.txt.gz | awk '/^#/ {sub(/#.*/,"");getline;}1 {print $1,$2-1,$2,$3}' | tr -s " " "\t" > $DATA/coverage.bed

#bills gtfs:

for chrom in {1..22} X Y ; do wget -P $DATA http://fastademo.bioch.virginia.edu/pfam_dna/Homo_sapiens.chr$chrom.pfam.gtf; done
for chrom in {1..22} X Y
do
sort -k2,2n $DATA/Homo_sapiens.chr$chrom.pfam.gtf | python rearrange.py | sort -k1,1 -k2,2n  > $DATA/chr$chrom.bed
done
cat $DATA/chr*.bed | sort -k1,1 -k2,2n > $DATA/all.bed
rm $DATA/chr*

# remove utrs and introns

bedtools intersect -a $DATA/all.bed -b $DATA/GRCh37.bed -wb | awk {'print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$17'} FS='\t' OFS='\t' | perl uniq.pl | perl -pe 's/"|;//g' > $DATA/foo.bed

# get number one appris transcripts by length or randomness

python seltrans.py <(sort -k2,2 -k3,3 $DATA/appris_data.principal.txt | grep -v ALTERNATIVE) <(sed '1d' $DATA/transcriptlengths.txt | sort -k1,1 -k2,2) > $DATA/transcripts.txt

# use appris to remove non-canonical transcripts; sort by ENSL gene_id and Pfam autoreg to merge uniqids into single domain

awk 'NR==FNR{a[$2]}$29 in a{print $0}' $DATA/transcripts.txt $DATA/foo.bed | tr -s " " "\t" | sort -k15,15 -k25,25 -k1,1 -k2,2n -k3,3n > $DATA/blah; mv $DATA/blah $DATA/foo.bed

# domain coverage and rearranging: // based on histograms, used 5x as a filter

python rearrange2.py <(awk '{if ($15==$29) print}' $DATA/foo.bed | tr -s " " "\t" | cut -f -45) | sort -k1,1 -k2,2n > $DATA/bar

bedtools intersect -a <(cat $DATA/bar | tr -s "\t" " " | cut -d " " -f 1-45 | tr -s " " "\t") -b <(awk '{if ($4>=5) print}' $DATA/coverage.bed) -wa -wb -sorted \
| awk '{ct[$1 $2 $3 $45]++; len[$1 $2 $3 $45]=$4; row[$1 $2 $3 $45]=$0} END {for (i in ct) print row[i],(ct[i]==0 ? ct[i]=0: ct[i]), (len[i]==0 ? len[i]=1: len[i]),ct[i]/(len[i]==0 ? len[i]=1: len[i])}' \
| tr -s " " "\t" | cut -f 1-45,50- > $DATA/alluniq.bed; rm $DATA/bar

# nodoms
perl -pe 's/tag\s*\S*?(?=\n|\s)//g' $DATA/GRCh37.bed | perl -pe 's/ccds_id\s*\S*?(?=\n|\s)//g' | perl -pe 's/"|;//g' | awk 'NR==FNR{a[$2]}$10 in a{print}' $DATA/transcripts.txt - | tr -s " " "\t" | sort -k10,10 > $DATA/foo2
awk '{split($15,trans,","); for (i in trans) print $1,$2,$3,$13,$25,trans[i],$43}' $DATA/alluniq.bed | sort -k6,6 | tr -s " " "\t"  > $DATA/foo3
python nodom.py $DATA/foo2 $DATA/foo3 > $DATA/foo; rm $DATA/foo2 $DATA/foo3

bedtools intersect -a <(sort -k1,1 -k2,2n $DATA/foo) -b <(awk '{if ($4>=5) print}' $DATA/coverage.bed) -wa -wb -sorted \
| awk '{ct[$1 $2 $3 $8 $24]++; len[$1 $2 $3 $8 $24]=$4; row[$1 $2 $3 $8 $24]=$0} END {for (i in row) print row[i],(ct[i]==0 ? ct[i]=0: ct[i]),(len[i]==0 ? len[i]=1: len[i]),ct[i]/(len[i]==0 ? len[i]=1: len[i])}' \
| tr -s " " "\t" | cut -f -24,29- | sort -k1,1 -k2,2n > $DATA/nodom.bed; rm $DATA/foo

# all regions (not necessary for the pipeline, just in case someone wants it)

cat $DATA/alluniq.bed $DATA/nodom.bed | sort -k1,1 -k2,2n | cat <(printf "#header for nodoms:\n#chr,start,end,length,info\n#info field contains gene_id, transcript_id, exon_number, gene_name, gene_source, gene_biotype, transcript_name, transcript_source, exon_id, near_pfamA_id (if applicable, describes what domain it is near), uniq_id\n#header for domains separated by exon:\n#chr,start,end,length,pfam_database_ver,type_of_sequence,blank_field,strand,blank_field2,info\n#info field contains pfamA_id, gene_name, transcript_id, protein_id, pfamseq_acc, pfamseq_id, pfamA_acc, pfamA_auto_reg, gene_id, matched_transcript_id, expn_number, matched_gene_name, gene_source, gene_biotype, transcript_name, transcript_source, exon_id, uniq_id, ccds_id (if applicable)\n") - > $DATA/allregions.bed

# do intersections for variants, merge lengths for autoreg+gene combos (lencount.py), get MAF, variant type (FG), gene, domain, chr, start, end, impact, other info for analysis

cat <(zgrep "^#" $DATA/VEPEXAC3.vcf.gz) <(zgrep -v "^#" $DATA/VEPEXAC3.vcf.gz | sort -k1,1 -k2,2n) > $DATA/VEPEXAC3.vcf; bgzip $DATA/VEPEXAC3.vcf
# filter if AN_Adj >= 97129.6:
cat <(zgrep "^#" $DATA/VEPEXAC3.vcf.gz) <(zgrep -v "^#" $DATA/VEPEXAC3.vcf.gz | awk -F ';' '{if ($1 ~ /X|Y/) {t=$17} else t=$16} {sub(/\w+=/,"",t); if (t>=0.8*60706*2) print}' | grep -w PASS | sort -k1,1 -k2,2n) | bgzip -c > $DATA/VEPEXAC3filter.vcf.gz
# singletons removed:
cat  <(zgrep "^#" $DATA/VEPEXAC3filter.vcf.gz) <(zgrep -v "^#" $DATA/VEPEXAC3filter.vcf.gz | awk -F ';' '{t=$1; sub(/.+=/,"",t); if (t!=1 && t!="1,1" && t!="1,1,1" && t!="1,1,1,1" && t!="1,1,1,1,1") print}' | sort -k1,1 -k2,2n) | bgzip -c > $DATA/VEPEXAC3nosingle.vcf.gz

#intersect

bedtools intersect -a <(cut -f 1,2,3,11,25,27,33,43,45,47,48 $DATA/alluniq.bed | python lencount.py | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -sorted -wb | cut -f 1,2,3,4,5,6,7,8,9,13,14,17 | python var.py -d > $DATA/uniqintfilter.bed
bedtools intersect -a <(cut -f 1,2,3,12,24,26,27 $DATA/nodom.bed | awk '{t=$7;$7=$6;$6=t;print}' | tr -s " " "\t" | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -sorted -wb | cut -f 1,2,3,4,5,6,7,11,12,15 | python var.py -n > $DATA/nodomintfilter.bed

# assign types to impacts

cat $DATA/uniqintfilter.bed $DATA/nodomintfilter.bed | cat <(printf "#chr,start,end,ref,alt,pfamA_id,autoreg,uniqid,length_of_region,covratio,gene_symbol,maf,impact,type,codons,amino_acids,gene_id_csq,gene_symbol_csq,transcript_id_csq,exon_number_csq,polyphen,sift,protein_position,biotype\n") <(awk 'NR==FNR{a[$2]}$19 in a{print $0}' $DATA/transcripts.txt - | sort -k11,11 -k7,7) > $DATA/allintfilter.bed

# do variant analysis by gene and maf (filters out entries that are "na," i.e, not dn or ds)

GENE="FLG"; MOD=g; MAF1=0.01 MAF2=''
bash lollipop.sh $DATA $OUT $GENE $MOD $MAF1 $MAF2

# sort domain occurrence count from bill, remove weird carriage return characters that screw things up

sort -k2,2 $DATA/human_pfam.counts | perl -pe 's/\r$//g'  > $DATA/blah; mv $DATA/blah $DATA/human_pfam.counts

# by uniqid a table

sed '1d' $DATA/allintfilter.bed | python table.py > $DATA/regionstable.txt

cat <(bedtools intersect -a <(awk '{$13=$33; print $0}' OFS="\t" $DATA/alluniq.bed | cut -f 1,2,3,11,13,25,27,43,45,47,48 | python lencount.py | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -v -sorted | cut -f 1,2,3,4,5,6,7,8,9,10,14,15,18) \
<(bedtools intersect -a <(cut -f 1,2,3,12,24,26,27 $DATA/nodom.bed | awk '{t=$8;$8=$7;$7=t;print}' | tr -s " " "\t" | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf.gz -v -sorted | cut -f 1,2,3,4,5,6,7,8,12,13,16 | awk '{print $1,$2,$3,".",$4,$5,$5,$6,$7,$8}') > $DATA/nointregions.txt
python noints.py $DATA/regionstable.txt $DATA/nointregions.txt
awk 'NR==FNR{a[$2]=$3} NR!=FNR{if ($1 in a) print $0"\t"a[$1]; else print $0"\t"0}' $DATA/human_pfam.counts $DATA/regionstable.txt > $DATA/blah; mv $DATA/blah $DATA/regionstable.txt

# make MAF spanning file

cat <(gawk 'NR==FNR{a[$5 $6][$1 $2 $3]=$1 " " $2 " " $3 " " $4} NR!=FNR{if ($3 $2 in a) {for (i in a[$3 $2]) print a[$3 $2][i],$0}}' <(cut -f 1,2,3,15,25,33,48 $DATA/alluniq.bed) <(cut -d " " -f 1,2,3,7,8,9,10,11 $DATA/regionstable.txt)) \
<(gawk 'NR==FNR{a[$5 $6][$1 $2 $3]=$1 " " $2 " " $3 " " $4} NR!=FNR{if ($3 $2 in a) {for (i in a[$3 $2]) print a[$3 $2][i],$0}}' <(cut -f 1,2,3,8,12,24,27 $DATA/nodom.bed | awk '{t=$6; $6=$5; $5=t; print}' OFS='\t') <(cut -d " " -f 1,2,3,7,8,9,10,11 $DATA/regionstable.txt)) \
| tr -s " " "\t" | sort -k1,1 -k2,2n > $DATA/regioncoordsdnds.bed
sort -k1,1 -k2,2n $DATA/allintfilter.bed > $DATA/blah; mv $DATA/blah $DATA/allintfilter.bed
cat <(bedtools intersect -a $DATA/regioncoordsdnds.bed -b $DATA/allintfilter.bed -wa -wb -sorted | awk '{if ($7==$23) print}' | cut -f -16,28,29,30) <(bedtools intersect -a $DATA/regioncoordsdnds.bed -b $DATA/allintfilter.bed -v -sorted | awk '{print $0"\t.\t.\t."}') \
| bedtools groupby -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 -c 17,18,19 -o collapse,collapse,collapse | sort -k1,1 -k6,6 -k5,5 -k7,7 \
| bedtools groupby -i - -g 1,6,5,7 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,3 -o distinct,min,max,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,collapse,collapse,collapse,collapse,collapse | cut -f 5-26 \
| python recompute.py | cat <(printf "#chr\tstart\tend\ttranscript\tdomain\tgene\tautoregs(uniqid for non-domain regions)\tcov_ratio\tlength\tdn\tds\tna\tdn/ds\tdensity\tfvrv\tprevalence\tmafs\timpacts\ttype\tstarts\tends\tmaf_modifier\n") - | bgzip > $DATA/regionsmafsdnds.bed.gz #also adds noint regions

# by gene a table

cat <(awk '{gene[$5]+=$10} END {for (i in gene) print i,gene[i]}' $DATA/nointsnodom.txt) <(awk '{gene[$5]+=$10} END {for (i in gene) print i,gene[i]}' $DATA/nointsuniq.txt) | tr -s " " "\t" | sort -k1,1 | bedtools groupby -g 1 -c 2 -o sum > $DATA/nointlen
awk '{gene[$5]+=$10} END {for (i in gene) print i,gene[i]}' $DATA/nointsuniq.txt | sort > $DATA/nointlenuniq

python genetable.py $DATA/allintfilter.bed | sort -k1,1 > $DATA/gene
awk 'NR==FNR{a[$1]=$2} NR!=FNR{{if ($1 in a) $2=$2+a[$1]} {print}}' $DATA/nointlen $DATA/gene > $DATA/genetablefilter.txt

python genetable.py <(grep -v NoDom $DATA/allintfilter.bed) | sort -k1,1 > $DATA/uniqgene
awk 'NR==FNR{a[$1]=$2} NR!=FNR{{if ($1 in a) $2=$2+a[$1]} {print}}' $DATA/nointlenuniq $DATA/gene > $DATA/uniqgenetablefilter.txt

#get 1000 genomes data and get clinvar data

wget -P $DATA ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/clinical_vcf_set/clinvar.vcf.gz

for chrom in {1..22} X Y ; do wget -P $DATA ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/functional_annotation/unfiltered/ALL.chr$chrom.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz; done

for chrom in {1..22} X ; do grep -v "^#" <(zcat $DATA/ALL.chr$chrom.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz) >> $DATA/1000genomes.vcf; done
grep -v "^#" <(zcat $DATA/ALL.chrY.phase3_integrated_func_anno.20130502.sites.vcf.gz) >> $DATA/1000genomes.vcf
cat <(grep "^#" <(zcat $DATA/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz)) <(sort -k1,1 -k2,2n <(grep -v "^#" $DATA/1000genomes.vcf)) > $DATA/blah; mv $DATA/blah $DATA/1000genomes.vcf

# filtering clinvar data and filtering 1000 genomes data

python varcomp.py -c -f $DATA/clinvar.vcf.gz | cat <(grep "^#" <(zcat $DATA/clinvar.vcf.gz)) <(sort -k1,1 -k2,2n -) > $DATA/clinvardeleterious.vcf

python varcomp.py -g -f $DATA/1000genomes.vcf $DATA/transcripts.txt | cat <(grep "^#" <(zcat $DATA/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz)) - > $DATA/1kgenomesfilter.vcf

# getting histogram lists for 1000genomes and clinvar

bedtools intersect -a <(awk 'NR==FNR{a[$7]=$0} NR!=FNR{if ($7 in a) print a[$7],$12,$13}' $DATA/regioncoordsdnds.bed <(zcat $DATA/regionsmafsdnds.bed.gz) | tr -s " " "\t" | sed '1d' | sort -k1,1 -k2,2n) -b $DATA/clinvardeleterious.vcf -wb -sorted | cut -f 1-16 > $DATA/clinvarhist.bed

bedtools intersect -a <(awk 'NR==FNR{a[$7]=$0} NR!=FNR{if ($7 in a) print a[$7],$12,$13}' $DATA/regioncoordsdnds.bed <(zcat $DATA/regionsmafsdnds.bed.gz) | tr -s " " "\t" | sed '1d' | sort -k1,1 -k2,2n) -b $DATA/1kgenomesfilter.vcf -wb -sorted | cut -f 1-16 > $DATA/1khist.bed

# divergence by domain table

sed '1d' $DATA/dtablemaf.g.01-.txt | awk '{if ($6>0) print $0}' | python diverge.py -d | awk '{if ($3>3) print $0}' | awk '{gene[$1]=$10; row[$1]=$0} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]+=1}} END {for (i in rank) print row[i],rank[i]/NR}' | tr -s " " "\t" | sort -k11,11nr  > $DATA/diverge.01.txt

# creates gene-protein pair tables and compares z score/mad metric to rvis and omim
# uses omim genemap of all genes and phenotypes for those genes as well as a list of genes that hit certain keywords: recessive, haploinsufficient, de novo, dominant negative, autosomal dominant/heterozygous mut 

OMIM=~/work/omim

#python omim.py $OMIM/genemap2.txt $DATA/foo.txt

# awk 'NR==FNR{a[$1]} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/denovo.txt)) $DATA/foo.txt \
# | awk 'NR==FNR{a[$1]} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/domneg.txt)) - \
# | awk 'NR==FNR{a[$1]} {if ($1 in a) print $0,"y"; else print $0 "n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/haplo.txt)) - \
# | awk 'NR==FNR{a[$1]} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/recessive.txt)) - \
# | awk 'NR==FNR{a[$1]} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/autodom.txt)) - >$DATA/omim.txt


bash div.sh $DATA/var.db $DATA/omim.txt g 0.00001 

# to make limited dtable-based dn/ds distributions and files use:

bash limit.sh $DATA g 0.001

# high dn/ds regions: 
# test for hgtables repeats and generate dn/ds file with coords

awk '{if ($9 >= 10) print $0}' $DATA/domaincoordsdnds.bed | wc -l # how many domain regions for each dn/ds?
grep -v "_ND" $DATA/uniqtable.txt | awk '{if ($11 >= 10) print $0}' | wc -l # how many autoregs for each dn/ds?
dir=($(ls $DATA/hg*))
dir[6]=$DATA/LCR-hs37d5.bed
for ((i=0; i<${#dir[@]}; i++)); do bedtools intersect -a $DATA/domaincoordsdnds.bed -b <(perl -pe 's/^chr(.*)$/$1/g' ${dir[$i]} | sort -k1,1 -k2,2n | bedtools merge) -wo -sorted | awk '{print $0,$3-$2}' OFS='\t' | sort -k4,4 -k5,5 -k6,6 | bedtools groupby -g 4,5,6,7,8,9,10 -c 14,15 -o sum,sum | awk '{print $0,$8/$9}' OFS='\t' > $SCRATCH/$(echo ${dir[$i]} | perl -pe 's/.*\/(.*)\..*/$1/g')'.txt'; done #perl -pe 's/^chr(.*)$/$1/g' (clips off the chr part) other perl script pastes file name into .txt
dir=($(ls $SCRATCH/hg*))
dir[6]=$SCRATCH/LCR-hs37d5.txt
for ((i=0; i<${#dir[@]}; i++)); do awk '{if ($7>=10) print $0}' ${dir[$i]} | wc -l; done #in alpha-order, interrupted, micro, repeats, segmental, selfchain, simple, then LCR
# test for k-mer alignability/blacklist/excludable regions
# running on sbatch/slurm: for ((i=0; i<${#dir[@]}; i++)); do sbatch test.sh $i; done
dir=($(ls $SCRATCH/wg*.bed*))
for ((i=0; i<${#dir[@]}; i++)); do bedtools intersect -a $DATA/domaincoordsdnds.bed -b <(awk '{if ($4<=0.1) print $0}' ${dir[$i]} | perl -pe 's/^chr(.*)$/$1/g' | sort -k1,1 -k2,2n | bedtools merge) -wo -sorted | awk '{print $0,$3-$2}' OFS='\t' | sort -k4,4 -k5,5 -k6,6 | bedtools groupby -g 4,5,6,7,8,9,10 -c 14,15 -o sum,sum | awk '{print $0,$8/$9}' OFS='\t' > $SCRATCH/$(echo ${dir[$i]} | perl -pe 's/.*\/(.*)\..*/$1/g')'.txt'; done
dir=($(ls $SCRATCH/wg*txt))
for ((i=0; i<${#dir[@]}; i++)); do awk '{if ($7>=10) print $0}' ${dir[$i]} | wc -l; done #in alpha-order, 100mer, 24, 36, 40, 50, 75, ConsensusExcludable, RegionsExcludable, Uniq20bp, Uniq35bp

dir=($(ls $SCRATCH/wg*.bed*; ls $DATA/LCR*; ls $DATA/hg* | grep -v selfchain))
for ((i=0; i<${#dir[@]}; i++)); do cut -f 1,2,3 ${dir[$i]} | sort -k1,1 -k2,2n | bedtools merge >> $DATA/blah.bed; done
sort -k1,1 -k2,2n $DATA/blah.bed > $DATA/allrepeatsmaps.bed

#to get .2 none/all and unique all/none:
bedtools intersect -a $DATA/domaincoordsdnds.bed -b $DATA/allrepeatsmaps.bed -wo -sorted | awk '{print $0,$3-$2}' OFS='\t' | sort -k4,4 -k5,5 -k6,6 | bedtools groupby -g 4,5,6,7,8,9,10 -c 14,15 -o sum,sum | awk '{print $0,$8/$9}' OFS='\t' > $SCRATCH/allrepeatsmaps.txt
sort -k4,4 -k5,5 -k6,6 $DATA/domaincoordsdnds.bed | bedtools groupby -g 4,5,6,7,8,9,10 -c 10 -o distinct | cut -f 1,2,3,4,5,6,7 | sort -k1,1 -k2,2 -k3,3 > foo
awk '{if ($10>=.2) print}' $SCRATCH/allrepeatsmaps.txt | cut -f 1,2,3,4,5,6,7 | sort -k1,1 -k2,2 -k3,3 > bar
grep -F -f bar foo # gives you the "all" intersections
grep -v -F -f bar foo # gives you the "none" intersections

#to get an idea of Zn:
grep -i -E 'ZNF|zf' bar
grep -i -E 'ZNF|zf' foo