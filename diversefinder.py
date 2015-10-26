import string
import sys
import numpy as np
import toolshed as ts
from collections import defaultdict
from math import sqrt

nodoms = defaultdict(dict)
domains = defaultdict(dict)
it = ts.reader(sys.argv[1])
for iv in it:
    if iv['domain'] == '.':
        nodoms[iv['transcript']][iv['autoregs']+iv['start']] = float(iv['dn_ds'])
    else:
        domains[iv['transcript']][iv['autoregs']+iv['start']] = float(iv['dn_ds'])

nstats = defaultdict(dict);
dstats = defaultdict(dict);
for trans in nodoms:
    c = [nodoms[trans][x] for x in nodoms[trans]]
    nstats[trans]['mean'] = np.mean(c)
    nstats[trans]['std'] = np.std(c)
    nstats[trans]['length'] = len(c)
for trans in domains:
    c = [domains[trans][x] for x in domains[trans]]
    dstats[trans]['mean'] = np.mean(c)
    dstats[trans]['std'] = np.std(c)
    dstats[trans]['length'] = len(c)

tstats = defaultdict(dict)
for trans in domains:
    N = dstats[trans]['length'] + nstats[trans]['length']
    tstats[trans]['N'] = N
    GM = (dstats[trans]['mean'] * dstats[trans]['length'] + \
    nstats[trans]['mean'] * nstats[trans]['length']) / N
    ESS = dstats[trans]['std']**2 * (dstats[trans]['length'] - 1) + \
    nstats[trans]['std']**2 * (nstats[trans]['length'] - 1)
    TGSS = (dstats[trans]['mean']-GM)**2 * nstats[trans]['length'] + \
    (nstats[trans]['mean']-GM)**2 * nstats[trans]['length']
    tstats[trans]['cstd'] = sqrt((ESS + TGSS) / (N - 1))
    print trans, tstats[trans]['N'], tstats[trans]['cstd']

    

    
