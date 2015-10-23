import sys
import string
from toolshed import reader, header, nopen

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-p",
    dest="pfam_counts_file_name",
    help="Pfam counts")

parser.add_option("-r",
    dest="fvrv_threshold",
    type=float,
    default=0.001,
    help="Fraction of rare variants threshold (default 0.001)")

parser.add_option("-m",
    dest="maf_cutoff",
    type=float,
    default=0.0,
    help="MAF threshold for dn/ds (default 0.0)")


(options, args) = parser.parse_args()

if not options.pfam_counts_file_name:
    parser.error('pfam count file not given')

f = open(options.pfam_counts_file_name,'r')

pfam = {}
for l in f:
    A = l.rstrip().split('\t')
    pfam[A[1]] = A[2]
f.close()



transcripts = {}

for d in reader('-'):
    if d['transcript_id'] not in transcripts:
        transcripts[d['transcript_id']] = {}
    if d['pfamA_auto_reg'] not in transcripts[d['transcript_id']]:
        transcripts[d['transcript_id']][d['pfamA_auto_reg']] = []

    transcripts[d['transcript_id']][d['pfamA_auto_reg']].append(d)

print '\t'.join(['#chr',
                 'start',
                 'end',
                 'transcript',
                 'domain',
                 'gene',
                 'autoregs',
                 'cov_ratio',
                 'length',
                 'dn',
                 'ds',
                 'na',
                 'dn_ds',
                 'density',
                 'fvrv',
                 'prevalence',
                 'mafs',
                 'impacts',
                 'type',
                 'starts',
                 'ends',
                 'maf_modifier',
                 'pos'])

for t in transcripts:
    for a in transcripts[t]:
        dn = 0
        ds = 0
        na = 0
        impacts = []
        poss = []
        mafs = []
        consqs = []
        starts = []
        starts_ends = set()
        ends = []
        for d in sorted(transcripts[t][a], key=lambda d: (int(d['start']),int(d['POS']) if d['POS'] != '.' else '.')):
            starts.append(int(d['start']))
            ends.append(int(d['end']))
            starts_ends.add((int(d['start']),int(d['end'])))

            if d['impact'] != '.':
                impacts.append(d['impact'])

                for impact in d['impact'].split('|'):
                    if impact == 'dn':
                        dn += 1
                    elif impact == 'ds':
                        ds += 1
                    elif impact == 'na':
                        na += 1
                    else:
                        sys.stderr.write('Unknown impact value "'+ d['impact'] + '"\n')
                        exit(1)

            if d['POS'] != '.':
                poss.append(d['POS'])

            if d['maf'] != '.':
                mafs.append(d['maf'])

            if d['Consequence'] != '.':
                consqs.append(d['Consequence'])

        d['dn'] = dn
        d['ds'] = ds
        d['na'] = na
        d['dn/ds'] = float(dn)/(float(ds) or 1.0 )
        d['density'] = (float(dn) + float(ds) + float(na)) / (float(d['len']) or 1.0)

        num_rv = 0
        tot_v = 0
        for maf in mafs:
            tf = 0
            for f in [float(x) for x in maf.split('|')]:
                tf+=f
            if tf <= options.fvrv_threshold:
                num_rv +=1
            tot_v +=1
        f_rv = float(num_rv)/(float(tot_v) or 1.0)

        pfam_prevl = "1"
        if d['pfamA_id'] in pfam:
            pfam_prevl = pfam[d['pfamA_id']]

        starts_a = []
        ends_a = []
        for start_end in sorted(starts_ends,key=lambda tup:tup[1]):
            starts_a.append(start_end[0])
            ends_a.append(start_end[1])

        print '\t'.join([\
                d['chr'],\
                str(min(starts)),\
                str(max(ends)),\
                d['transcript_id'],
                d['pfamA_id'],
                d['gene_name'],
                d['pfamA_auto_reg'],
                d['ratio_of_cover_to_total'],
                d['len'],
                str(d['dn']),
                str(d['ds']),
                str(d['na']),
                str(d['dn/ds']),
                str(d['density']),
                str(f_rv),
                pfam_prevl,
                ','.join(mafs) if len(mafs) > 0 else '.',
                ','.join(impacts) if len(impacts) > 0 else '.',
                ','.join(consqs) if len(consqs) > 0 else '.',
                ','.join([str(x) for x in starts_a]),
                ','.join([str(x) for x in ends_a]),
                str(options.maf_cutoff),
                ','.join(poss) if len(poss) > 0 else '.',
                ])



#class Record(object):
#	def __init__(self,fields):
#		self.chr = fields[0]
#		self.start = fields[1]
#		self.end = fields[2]
#		self.ref = fields[3]
#		self.alt = fields[4]
#		self.domain = fields[5]
#		self.autoreg = fields[6] #uniqid for nodoms
#		self.uniqid = fields[7]
#		self.len = fields[8]
#		self.covratio = fields[9]
#		self.gene = fields[10]
#		self.maf = fields[11]
#		self.type = fields[12]
#		self.impact = fields[13]
#		self.codonchange = fields[14]
#		self.aachange = fields[15]
#		self.geneidcsq = fields[16]
#		self.genecsq = fields[17] # can be different from the gene name from ensembl exon file
#		self.transidcsq = fields[18] # only picked those that match our transcript list, but may be a versioning issue if gene names don't match
#		self.exonnum = fields[19]
#		self.polyphen = fields[20]
#		self.sift = fields[21]
#		self.position = fields[22]
#		self.biotype = fields[23]
#
#ct={} # will be an autoreg/gene key that refers to a list of three values: [0] dn, [1] ds, [2] na
#oldauto=None
#oldgene=None
#for line in sys.stdin:
#	try:
#		oldauto=r_.autoreg
#		oldgene=r_.gene
#	except:
#		pass
#	r_=Record(line.rstrip().split("\t"))
#	types=r_.type.split("|")
#	if r_.autoreg+r_.gene in ct:
#		if "dn" in types:
#			ct[r_.autoreg+r_.gene][0]+=1
#		if "ds" in types:
#			ct[r_.autoreg+r_.gene][1]+=1
#		if "na" in types:
#			ct[r_.autoreg+r_.gene][2]+=1
#	else:
#		ct[r_.autoreg+r_.gene]=[float(0),float(0),float(0),float(r_.len),"\t".join([r_.domain,r_.gene,r_.autoreg,r_.uniqid,r_.covratio,r_.len])]
#		if "dn" in types:
#			ct[r_.autoreg+r_.gene][0]+=1
#		if "ds" in types:
#			ct[r_.autoreg+r_.gene][1]+=1
#		if "na" in types:
#			ct[r_.autoreg+r_.gene][2]+=1
#	if oldauto!=None and oldauto!=r_.autoreg:
#		if ct[oldauto+oldgene][1]==0:
#			dnds=str(int(ct[oldauto+oldgene][0]/1))
#		else:
#			dnds=str(round(ct[oldauto+oldgene][0]/ct[oldauto+oldgene][1],4))
#		density=str(round((ct[oldauto+oldgene][0]+ct[oldauto+oldgene][1]+ct[oldauto+oldgene][2])/ct[oldauto+oldgene][3],4))
#		print "\t".join([ct[oldauto+oldgene][4],str(int(ct[oldauto+oldgene][0])),str(int(ct[oldauto+oldgene][1])),str(int(ct[oldauto+oldgene][2])),dnds,density])
