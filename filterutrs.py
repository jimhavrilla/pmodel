import string
import sys
import pybedtools
import subprocess 
import operator 

class Record(object):
    def __init__(self, fields):
        self.chr = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.crapandstrand = "\t".join(fields[3:6])
        self.info=fields[6]
        self.transid = self.info.replace("\"","").replace(";","").split()[3]
    def __str__(self):
        return '\t'.join([self.chr, \
                          str(self.start),\
                          str(self.end), \
                          self.crapandstrand,\
                          self.info])

E = {}
U = {}

f = open(sys.argv[1],"r")
for l in f:
    e = Record(l.rstrip().split("\t"))
    if not e.transid in E:
        E[e.transid] = []
    E[e.transid].append(e)
f.close()

f = open(sys.argv[2],"r")
for l in f:
    u = Record(l.rstrip().split("\t"))
    if not u.transid in U:
        U[u.transid] = []
    U[u.transid].append(u)
f.close()

total_dropped = 0

for e_transid in E:
    if not e_transid in U:
        x=1
        for e in E[e_transid]:
            print e
    else:
        f = open('u.tmp', 'w')
        f.write('\n'.join( \
                [str(u) for u in sorted(U[e_transid], \
                                        key=lambda x:(x.chr,x.start,x.end))]))
        f.close()

        f = open('e.tmp', 'w')
        f.write('\n'.join( \
                [str(e) for e in sorted(E[e_transid], \
                                        key=lambda x:(x.chr,x.start,x.end))]))

        f.close()

        f = open('r.tmp', 'w')
        p = subprocess.Popen('bedtools subtract -sorted -a e.tmp -b u.tmp',
                             shell=True,
                             stdout=f)
        p.wait()
        f.close()

        f = open('r.tmp', 'r')
        num_lines = sum(1 for line in f)
        f.close()

        total_dropped += len(E[e_transid]) - num_lines

#        if (len(E[e_transid]) != num_lines):
#            sys.stderr.write('Warning: ' + \
#                             str(len(E[e_transid])) + \
#                             ' exons in transcript and ' + \
#                             str(num_lines) + \
#                             ' exons after UTR removal for ' + \
#                             e_transid + '\n')

        f = open('r.tmp', 'r')
        for l in f:
            print l.rstrip()
        f.close()

sys.stderr.write('A total of ' + str(total_dropped) + 'exons removed' + '\n')
