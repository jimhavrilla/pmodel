import string
import sys
import pybedtools
from collections import defaultdict
from itertools import groupby
import re

class Exon(object):
    def __init__(self, fields):
        self.chr = fields[0]
        self.start = fields[1]
        self.end = fields[2]
        self.info = "\t".join(fields[6:24]).rstrip()
        self.uniqid = fields[7].strip(";").strip("\"")+"_"+fields[9].strip(";").strip("\"")+"_"+fields[23].rstrip().strip(";").strip("\"")+"_NoDom"
        self.transid = fields[9].strip(";").strip("\"")

    def __str__(self):
        return "\t".join([self.chr, self.start, self.end, self.info, self.uniqid])

    def __len__(self):
        return int(self.end) - int(self.start)

class Uniq(object):
    def __init__(self, fields):
        self.chr = fields[0]
        self.start = fields[1]
        self.end = fields[2]
        self.gene = fields[3]
        self.autoreg = fields[4]
        self.transid = fields[5]
        self.exonid = fields[6]

    def __str__(self):
        return "\t".join([self.chr, self.start, self.end])

exons_by_transcript_id = defaultdict(list)
for ex in (Exon(x.rstrip().split("\t")) for x in open(sys.argv[1])):
    exons_by_transcript_id[ex.transid].append(ex)

domains_by_transcript_id = defaultdict(list)
for dom in (Uniq(x.rstrip().split("\t")) for x in open(sys.argv[2])):
    domains_by_transcript_id[dom.transid].append(dom)

grp_u = lambda b: (b['uniqid'])
grp_start = lambda b: (b['uniqid'], b['start'])


def genuniq(a, b):
    for i in (x.split("\t") for x in str(a.subtract(b)).splitlines()):
        chrom, start, end = i[:3]
        info = "\t".join(i[3:len(i)-1])
        uniqid = i[len(i)-1]
        yield dict(chrom=chrom, start=int(start), end=int(end), info=info,
                   uniqid=uniqid)

for transcript, exons in exons_by_transcript_id.iteritems():
    doms = domains_by_transcript_id[transcript]

    if len(doms) == 0:
        ct = "1"
        for e in exons:
            #print len(e), e.start, e.end, e.chr, ct
            print "\t".join([e.chr, e.start, e.end, str(len(e)), e.info + " uniqid "+e.uniqid + ct])

        continue

    a = pybedtools.BedTool("\n".join(str(e) for e in exons), from_string=True)
    b = pybedtools.BedTool("\n".join(str(d) for d in doms), from_string=True)

    # sort by start, uniq then group by uniq and make sur ethey all have uniq ct
    for (grp, li) in groupby(sorted(genuniq(a, b), key=grp_start), grp_u):
        for ct, exon_piece in enumerate(li, start=1):
            e = exon_piece
            print("\t".join([e['chrom'] , str(e['start']), str(e['end']), str(e['end'] - e['start']),
                                        e['info'] + " uniqid " + e['uniqid'] + str(ct)]))
