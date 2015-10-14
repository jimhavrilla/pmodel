from __future__ import division
from collections import namedtuple
import string
import math
from collections import Counter, defaultdict
import toolshed as ts

def ientropy(vals):
    p, lns = Counter(vals), float(len(vals))
    return -sum( count/lns * math.log(count/lns) for count in p.values())

def entropy(domains):
    prevalence = len(domains) 
    vals = [len(d.mafs.split(",")) for d in domains]
    e = ientropy(vals) / (math.log(prevalence) or 1)
    return e
   
if __name__ == "__main__":

    by_domain = defaultdict(list)
    for region in ts.reader('$DATA/regionsmafsdnds.bed.gz', header=namedtuple):
        by_domain[region.domain if region.domain != "." else region.autoregs].append(region)
            
    for dname, domains in by_domain.iteritems():
        print dname, entropy(domains)
