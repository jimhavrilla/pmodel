
from swindow import JimFile, windower
from math import log10
import numpy as np
import sys

myfile = sys.argv[1]
iterator = JimFile(myfile) #/uufs/chpc.utah.edu/common/home/u6000294/lustre/u6000294/pmodel/y.sort.bed.gz


def FRV(intervals, maf_cutoff=1.0/10000):
    return sum(1.0 for iv in intervals if iv.aaf <= maf_cutoff) / len(intervals)

def IAFI(intervals, n_samples=65000):
    minaf = 1.0 / (2 * n_samples + 1)
    val = sum(1.0/max(iv.aaf, minaf) for iv in intervals)
    return log10(val / len(intervals))

def dnds_ratio(intervals):
    dn, ds, na = 0, 0, 0
    values = defaultdict(list)  
    for iv in intervals:
        if iv.dnds is None:
            continue
        dnds = iv.dnds.split('|')
        dn += dnds.count('dn')
        ds += dnds.count('ds')
        na += dnds.count('na')
    values['dn'] = dn
    values['ds'] = ds
    values['na'] = na
    values['dnds'] = dn / float(ds or 1.0)
    return values


def bytranscriptdist(grp, inext):
    """ group by transcript and split at gaps > 50 bases"""
    return inext.transcript != grp[0].transcript \
            or inext.start - grp[-1].end > 50

def smallchunk(grp, inext):
    return len(grp) > 20 or inext.transcript != grp[0].transcript \
        or inext.start - grp[-1].end > 20

def rescale(vals):
    minv, maxv = min(vals), max(vals)
    return [float(v - minv) / ((maxv - minv) or 1) for v in vals]

# allow a gap of at most 50 bases.
from collections import defaultdict
saved = defaultdict(list)
for chunk in windower(iterator, smallchunk):

    frv = FRV(chunk)
    iafi = IAFI(chunk)
    values = dnds_ratio(chunk)
    
    saved["chrom"].append(chunk[0].chrom)
    saved["start"].append(chunk[0].start)
    saved["end"].append(chunk[-1].end)
    saved["frv"].append(frv)
    saved["iafi"].append(iafi)
    saved["dnds"].append(values["dnds"])
    saved["dn"].append(values["dn"])
    saved["ds"].append(values["ds"])
    saved["na"].append(values["na"])
    saved["trans"].append(chunk[0].transcript)

saved["iafi"] = rescale(saved["iafi"])
saved["frv"] = rescale(saved["frv"])

for i in range(len(saved["chrom"])):
    print "%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%s" % (
            saved["chrom"][i],
            saved["start"][i],
            saved["end"][i],
            saved["iafi"][i],
            saved["frv"][i],
            saved["dnds"][i],
            saved["dn"][i],
            saved["ds"][i],
            saved["na"][i],
            saved["trans"][i])
