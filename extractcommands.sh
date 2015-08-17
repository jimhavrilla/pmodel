#python extract.py (old vcf extract prog)

export DATA=~/work/data/pmodeldata
export SOFTWARE=~/software
export $SCRATCH=~/work/scratch/analysis

#bill's domain count file:

mysql -N --raw -h wrpxdb.its.virginia.edu -u web_user -pfasta_www pfam27 < count_human_pfam.sql > human_pfam.counts

#for vcfs:

#grep ^# $DATA/ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf; grep -v ^# -h $DATA/ESP*.vcf) > $DATA/foo.vcf; mv $DATA/foo.vcf $DATA/ESPALL.vcf
#perl $SOFTWARE/ensembl-tools-release-76/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/ESPALL.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --total_length -o $DATA/VEPESPALL.vcf --vcf --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --force_overwrite --port 3337 --dir ~/.vep/ --compress "gunzip -c"
sudo perl $SOFTWARE/ensembl-tools-release-76/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/ExAC.r0.3.sites.vep.vcf.gz --cache --sift b --polyphen b --symbol --numbers --biotype --total_length -o $DATA/VEPEXAC3.vcf --vcf --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --force_overwrite --port 3337 --dir ~/.vep/ --compress "gunzip -c" --offline
gzcat $DATA/Homo_sapiens.GRCh37.75.gtf.gz | grep protein_coding$'\t'exon | perl -pe 's/protein_coding\texon\t//g' | grep -P '^\d*[0-9X-Y]\t' | perl gtf2bed.pl | sort -k10,10 > $DATA/exons.bed
gzcat $DATA/Homo_sapiens.GRCh37.75.gtf.gz | grep protein_coding$'\t'UTR | perl -pe 's/protein_coding\tUTR\t//g' | grep -P '^\d*[0-9X-Y]\t' | perl gtf2bed.pl | sort -k10,10 > $DATA/utrs.bed
python filterutrs.py $DATA/exons.bed $DATA/utrs.bed > $DATA/GRCh37.bed
#grep -v -E "C|G" $DATA/mart_export.bed | sort -k1,1 -k2,2n | sed '395352d' > foo.bed; mv foo.bed $DATA/mart_export.bed; # human genes from Ensembl

#for vcf coverage:
for chrom in {1..22} X Y ; do wget -P $DATA ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/coverage/Panel.chr$chrom.coverage.txt.gz; done
gzcat $DATA/Panel.chr*.coverage.txt.gz | awk '/^#/ {sub(/#.*/,"");getline;}1 {print $1,$2-1,$2,$3}' | tr -s " " "\t" > $DATA/coverage.bed
#bedtools intersect -a $DATA/VEPEXAC.vcf -b $DATA/coverage.bed -wa -wb | awk '{if ($12 >= 5) print $1,$2,$3,$4,$5,$6,$7,$8,$12}' OFS="\t" > $DATA/covVEPEXAC.vcf

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

awk 'NR==FNR{a[$2];next}$29 in a{print $0}' $DATA/transcripts.txt $DATA/foo.bed | tr -s " " "\t" | sort -k27,27 -k25,25 -k1,1 -k2,2n -k3,3n > $DATA/blah.txt; mv $DATA/blah.txt $DATA/foo.bed

# domain coverage and rearranging:  // based on histograms, used 5x as a filter

python rearrange2.py -u $DATA/foo.bed | sort -k1,1 -k2,2n > $DATA/bar

bedtools intersect -a <(cat $DATA/bar | tr -s "\t" " " | cut -d " " -f 1-45 | tr -s " " "\t") -b <(awk '{if ($4>=5) print}' $DATA/coverage.bed) -wa -wb -sorted \
| awk '{ct[$1 $2 $3 $45]++; len[$1 $2 $3 $45]=$4; row[$1 $2 $3 $45]=$0} END {for (i in ct) print row[i],(ct[i]==0 ? ct[i]=0: ct[i]), (len[i]==0 ? len[i]=1: len[i]),ct[i]/(len[i]==0 ? len[i]=1: len[i])}' \
| tr -s " " "\t" | cut -f 1-45,50- > $DATA/alluniq.bed; rm $DATA/bar

#old dom code

sort -k11,11 -k1,1 -k2,2n $DATA/foo.bed | python rearrange2.py -d > $DATA/alldom.bed

# nodoms

python nodom.py <(perl -pe 's/tag\s*\S*?(?=\n|\s)//g' $DATA/GRCh37.bed | perl -pe 's/ccds_id\s*\S*?(?=\n|\s)//g' | perl -pe 's/"|;//g' | awk 'NR==FNR{a[$2];next}$10 in a{print}' $DATA/transcripts.txt - | tr -s " " "\t" | sort -k10,10 | uniq) <(awk '{split($15,trans,","); for (i in trans) print $1,$2,$3,$13,$25,trans[i],$43}' $DATA/alluniq.bed | sort -k6,6 | tr -s " " "\t" | uniq) > $DATA/foo

bedtools intersect -a <(sort -k1,1 -k2,2n $DATA/foo | uniq) -b <(awk '{if ($4>=5) print}' $DATA/coverage.bed) -wa -wb -sorted \
| awk '{ct[$1 $2 $3 $8 $24]++; len[$1 $2 $3 $8 $24]=$4; row[$1 $2 $3 $8 $24]=$0} END {for (i in ct) print row[i],(ct[i]==0 ? ct[i]=0: ct[i]),(len[i]==0 ? len[i]=1: len[i]),ct[i]/(len[i]==0 ? len[i]=1: len[i])}' \
| tr -s " " "\t" | cut -f -24,29- | sort -k1,1 -k2,2n > $DATA/nodom.bed; rm $DATA/foo

cat $DATA/alluniq.bed $DATA/nodom.bed | sort -k1,1 -k2,2n | cat <(printf "#header for nodoms:\n#chr,start,end,length,info\n#info field contains gene_id, transcript_id, exon_number, gene_name, gene_source, gene_biotype, transcript_name, transcript_source, exon_id, near_pfamA_id (if applicable, describes what domain it is near), uniq_id\n#header for domains separated by exon:\n#chr,start,end,length,pfam_database_ver,type_of_sequence,blank_field,strand,blank_field2,info\n#info field contains pfamA_id, gene_name, transcript_id, protein_id, pfamseq_acc, pfamseq_id, pfamA_acc, pfamA_auto_reg, gene_id, matched_transcript_id, expn_number, matched_gene_name, gene_source, gene_biotype, transcript_name, transcript_source, exon_id, uniq_id, ccds_id (if applicable)\n") - > $DATA/allregions.bed

# do intersections for variants, merge lengths for autoreg+gene combos, get MAF and convert from percent to fraction, variant type (FG), gene, domain, chr, start, end, impact, other info for analysis

cat <(grep "^#" $DATA/VEPEXAC3.vcf) <(grep -v "^#" $DATA/VEPEXAC3.vcf | sort -k1,1 -k2,2n) > $DATA/foo; mv $DATA/foo $DATA/VEPEXAC3.vcf
# filter:
cat <(grep "^#" $DATA/VEPEXAC3.vcf) <(grep -v "^#" $DATA/VEPEXAC3.vcf | awk -F ';' '{if ($1 ~ /X|Y/) {t=$17} else t=$16} {sub(/\w+=/,"",t); if (t>=0.8*60706*2) print}' | grep -w PASS | sort -k1,1 -k2,2n) > $DATA/VEPEXAC3filter.vcf
cat  <(grep "^#" $DATA/VEPEXAC3filter.vcf) <(grep -v "^#" $DATA/VEPEXAC3filter.vcf | awk -F ';' '{t=$1; sub(/.+=/,"",t); if (t!=1 && t!="1,1" && t!="1,1,1" && t!="1,1,1,1" && t!="1,1,1,1,1") print}' | sort -k1,1 -k2,2n) > $DATA/VEPEXAC3nosingle.vcf
#intersectttttt
bedtools intersect -a <(sort -k1,1 -k2,2n -k3,3n $DATA/alldom.bed | tr -s " " "\t" | cut -f 1,2,3,11,13,45 ) -b $DATA/VEPEXAC3filter.vcf -sorted -wb | cut -f 1,2,3,4,5,6,10,11,14 | python var.py -d > $DATA/domint.bed
#13=33 makes the pfam gene the ensembl gene because the latter is the correct one.
bedtools intersect -a <(awk '{$13=$33; print $0}' OFS="\t" $DATA/alluniq.bed | cut -f 1,2,3,11,13,25,27,43,45,46,47,48 | python lencount.py | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf -sorted -wb | cut -f 1,2,3,4,5,6,7,8,9,10,14,15,18 | python var.py -d > $DATA/uniqintfilter.bed
bedtools intersect -a <(cut -f 1,2,3,12,24,25,26,27 $DATA/nodom.bed | awk '{t=$8;$8=$7;$7=t;print}' | tr -s " " "\t" | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3filter.vcf -sorted -wb | cut -f 1,2,3,4,5,6,7,8,12,13,16 | python var.py -n > $DATA/nodomintfilter.bed

cat $DATA/uniqintfilter.bed $DATA/nodomintfilter.bed | cat <(printf "#chr,start,end,ref,alt,pfamA_id,autoreg,uniqid,covct,length_of_region,covratio,gene_symbol,maf,impact,codons,amino_acids,gene_id_csq,gene_symbol_csq,transcript_id_csq,exon_number_csq,polyphen,sift,protein_position,biotype\n") - > $DATA/allintfilter.bed

# sort domain occurrence count from bill

sort -k2,2 $DATA/human_pfam.counts > $DATA/blah.bed; mv $DATA/blah.bed $DATA/human_pfam.counts

# get total bp for each domain, counts

cat $DATA/alldom.bed | tr -s " " "\t" | cut -f 4,11 | awk '{arr[$2]+=$1} END {for (i in arr) {print i,arr[i]}}' | sort -k1,1 > $DATA/sumlist.bed

# make db for queries AND filter variants by canonical transcripts

bash makedb.sh $DATA/varfilter.db <(awk 'NR==FNR{a[$2];next}$19 in a{print $0}' $DATA/transcripts.txt <(sed '1d' $DATA/allint.bed)) 

# pick one impact per variant for domains

bash query.sh $DATA/var.db '%' g 0 | sort -k6,6 -k8,8 | python formatvar.py > $DATA/allint2dom.bed

# for uniqids

bash query.sh $DATA/varfilter.db '%' g 0 | sort -t '|' -k7,7 -k1,1 -k2,2 -k3,3 | python formatvar.py > $DATA/allint2uniqfilter.bed

# do variant analysis by gene and maf (filters out entries that are "na," i.e, not dn or ds)

GENE="FLG"; MOD=g; MAF1=0.01 MAF2=''
bash lollipop.sh $DATA $OUT $GENE $MOD $MAF1 $MAF2

# make table of counts, non-syn, syn, total var, total bp/exome per domain

python maketable.py $DATA/allint2dom.bed > $DATA/foo.bed
python mergetable.py -f $DATA/foo.bed $DATA/human_pfam.counts $DATA/sumlist.bed > $DATA/dtable.txt; rm $DATA/foo.bed

# by uniqid a table

gawk 'function median(v) {c=asort(v,j); if (c % 2) return j[(c+1)/2]; else return (j[c/2+1]+j[c/2])/2.0} {{if ($15=="ds") {sct[$7 $12]++; ct[$7 $12]++} else if ($15=="dn") nct[$7 $12]++; ct[$7 $12]++} len[$7 $12]=$11; maf[$7 $12][$1 $2 $3 $15]=$13; row[$7 $12]=$6 " " $12 " " $7 " " $8 " " $9 " " $10 " " $11} END {for (i in ct) print row[i],(nct[i]==0 ? nct[i]=0: nct[i]),(sct[i]==0 ? sct[i]=0: sct[i]),(nct[i]+sct[i]),nct[i]/(sct[i]==0 ? sct[i]+1: sct[i]),(nct[i]+sct[i])/len[i],median(maf[i])}' <(grep -v NoDom $DATA/allint2uniqfilter.bed) > $DATA/uniqtablefilter.txt

gawk 'function median(v) {c=asort(v,j); if (c % 2) return j[(c+1)/2]; else return (j[c/2+1]+j[c/2])/2.0} {{if ($15=="ds") {sct[$8 $12]++; ct[$8 $12]++} else if ($15=="dn") nct[$8 $12]++; ct[$8 $12]++} len[$8 $12]=$11; maf[$8 $12][$1 $2 $3 $15]=$13; row[$8 $12]=$6 " " $12 " " $7 " " $8 " " $9 " " $10 " " $11} END {for (i in ct) print row[i],(nct[i]==0 ? nct[i]=0: nct[i]),(sct[i]==0 ? sct[i]=0: sct[i]),(nct[i]+sct[i]),nct[i]/(sct[i]==0 ? sct[i]+1: sct[i]),(nct[i]+sct[i])/len[i],median(maf[i])}' <(grep NoDom $DATA/allint2uniqfilter.bed) > $DATA/nodomtablefilter.txt

bedtools intersect -a <(awk '{$13=$33; print $0}' OFS="\t" $DATA/alluniq.bed | cut -f 1,2,3,11,13,25,27,43,45,46,47,48 | python lencount.py | sort -k1,1 -k2,2n -k3,3n) -b $DATA/VEPEXAC3.vcf -v -sorted | cut -f 1,2,3,4,5,6,7,8,9,10,14,15,18 > $DATA/noints.txt
python noints.py $DATA/uniqtablefilter.txt $DATA/noints.txt

# by gene a table

gawk 'function median(v) {c=asort(v,k); if (c % 2) return k[(c+1)/2]; else return (k[c/2+1]+k[c/2])/2.0} {{if ($15=="ds") {sct[$12]++; ct[$12]++} else if ($15=="dn") nct[$12]++; ct[$12]++} auto[$7]=$7; len[$12 $7]=$11; maf[$12][$1 $2 $3 $15]=$13; row[$12]=$12}
END {
        for (l in row) {
                for (j in auto) {
                        leng[l]+=len[l j]
                }
        }
}
END {
        for (i in row) {
                print row[i],leng[i],(nct[i]==0 ? nct[i]=0: nct[i]),(sct[i]==0 ? sct[i]=0: sct[i]),(nct[i]+sct[i]),nct[i]/(sct[i]==0 ? sct[i]+1: sct[i]),(nct[i]+sct[i])/(leng[i]==0? leng[i]+1: leng[i]),median(maf[i])
        }
}' $DATA/allint2uniqfilter.bed | tr -s " " "\t" | sort -k11,11nr > $DATA/genetablefilter.txt

# divergence by domain table

sed '1d' $DATA/dtablemaf.g.01-.txt | awk '{if ($6>0) print $0}' | python diverge.py -d | awk '{if ($3>3) print $0}' | awk '{gene[$1]=$10; row[$1]=$0} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]+=1}} END {for (i in rank) print row[i],rank[i]/NR}' | tr -s " " "\t" | sort -k11,11nr  > $DATA/diverge.01.txt

# creates gene-protein pair tables and compares z score/mad metric to rvis and omim
# uses omim genemap of all genes and phenotypes for those genes as well as a list of genes that hit certain keywords: recessive, haploinsufficient, de novo, dominant negative, autosomal dominant/heterozygous mut 

OMIM=~/work/omim

#python omim.py $OMIM/genemap2.txt $DATA/foo.txt

# awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/denovo.txt)) $DATA/foo.txt \
# | awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/domneg.txt)) - \
# | awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0 "n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/haplo.txt)) - \
# | awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/recessive.txt)) - \
# | awk 'NR==FNR{a[$1];next} {if ($1 in a) print $0,"y"; else print $0,"n"}' <(perl -pe 's/.*?;\s(\w*)\s.*/$1/g' <(sed '1,5d' $OMIM/autodom.txt)) - >$DATA/omim.txt


bash div.sh $DATA/var.db $DATA/omim.txt g 0.00001 

# to make limited dtable-based dn/ds distributions and files use:

bash limit.sh $DATA g 0.001

# high dn/ds regions: 
# test for hgtables repeats and generate dn/ds file with coords
gawk 'NR==FNR{a[$4 $5][$1 $2 $3]=$1 " " $2 " " $3} NR!=FNR{if ($3 $2 in a) {for (i in a[$3 $2]) print a[$3 $2][i],$0}}' <(cut -f 1,2,3,25,33 $DATA/alluniq.bed) <(cut -d " " -f 1,2,3,7,8,9,11 $DATA/uniqtablefilter.txt | grep -v -w '\.') | tr -s " " "\t" | sort -k1,1 -k2,2n > $DATA/domaincoordsdnds.bed
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
grep -i -E 'ZNF|zf'

# clans crap

# alldomint<-read.delim(paste(DATA,"/domint.bed",sep=""),header=FALSE)
# alldomint$V6=gsub(";","",alldomint$V6)
# m<-merge(alldomint,clans,by.x="V6",by.y="pfamA_id",all.x=TRUE)
# write.table(m,paste(DATA,"/alldomsclansint.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
# cat <(printf "domain\tchr\tstart\tend\tref\talt\tuniqid\tgene\teamaf\taamaf\tmaf\timpact\tcodon\taminoacid\tgeneid\tgenecsq\ttransid\texonno\tpolyphen\tsift\tproteinpos\tbiotype\tpfamacc\tclanacc\tclanid\tclandomcount\n") <(sed '1d' $DATA/alldomsclansint.txt) > $DATA/foo.txt; mv $DATA/foo.txt $DATA/alldomsclansint.txt

# perl -pe 's/ /\t/g' $DATA/alldom.bed | cut -f -45 > $DATA/foo.bed
# alldom<-read.delim(paste(DATA,"/foo.bed",sep=""),header=FALSE)
# alldom$V11=gsub(";","",alldom$V11)
# m<-merge(alldom,clans,by.x="V11",by.y="pfamA_id",all.x=TRUE)
# m<-m[,c(2:11,1,12:49)]
# write.table(m,paste(DATA,"/alldomsclans.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
# cat <(printf "chr\tstart\tend\tbplength\tdb\tseqtype\tblank\tstrand\tblank\tpfamA_id\tgene_name\ttranscript_id\tprotein_id\tpfamseq_acc\tpfamseq_id\tpfamA_acc\tpfamA_auto_reg\tgene_id\ttranscript_id2\texon_number\tgene_name2\tgene_source\tgene_biotype\ttranscript_name\ttranscript_source\texon_id\tuniq_id\tpfamA_acc\tclan_acc\tclan_id\tdom_cnt\n") <(sed '1d' $DATA/alldomsclans.txt) > $DATA/foo.txt; mv $DATA/foo.txt $DATA/alldomsclans.txt