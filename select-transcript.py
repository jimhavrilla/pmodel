import sys
import string

class Appris(object):
    def __init__(self,fields):
        self.gene=fields[0]
        self.geneid=fields[1]
        self.transid=fields[2]
        self.ccdsid=fields[3]
        self.type=fields[4]
        self.length = 0

class Transcript(object):
    def __init__(self,fields):
        self.geneid=fields[0]
        self.transid=fields[1]
        self.gene=fields[2]
        self.length=int(fields[3])

import collections
from operator import attrgetter

appris_objs = collections.defaultdict(dict)
for obj in (Appris(line.rstrip().split("\t")) for line in open(sys.argv[1])):
    appris_objs[obj.geneid][obj.transid] = obj

for trans in (Transcript(line.rstrip().split("\t")) for line in open(sys.argv[2])):
    if not trans.geneid in appris_objs:
        continue
    if not trans.transid in appris_objs[trans.geneid]: continue
    appris_objs[trans.geneid][trans.transid].length = trans.length

for geneid, transcript_dict in appris_objs.iteritems():
    aprs = transcript_dict.values()
    if len(aprs) != 1:
        aprs.sort(key=attrgetter('length'), reverse=True)

    print "\t".join([aprs[0].geneid, aprs[0].transid, str(aprs[0].length)])
