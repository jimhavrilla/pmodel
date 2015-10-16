
m swindow import JimFile, windower
from math import log10

iterator = JimFile('/uufs/chpc.utah.edu/common/home/u6000294/lustre/u6000294/pmodel/y.sort.bed.gz')


def FRV(intervals, maf_cutoff=1.0/10000):
    return sum(1.0 for iv in intervals if iv.aaf <= maf_cutoff) / len(intervals)

def IAFI(intervals, n_samples=65000):
    minaf = 1.0 / (2 * n_samples + 1)
    val = sum(1.0/max(iv.aaf, minaf) for iv in intervals)
    return log10(val / len(intervals))


def bytranscriptdist(grp, inext):
    """ group by transcript and split at gaps > 50 bases"""
    return inext.transcript != grp[0].transcript \
            or inext.start - grp[-1].end > 50

# allow a gap of at most 50 bases.
for chunk in windower(iterator, bytranscriptdist):
    frv = FRV(chunk)
    iafi = IAFI(chunk)
    print "%s\t%d\t%d\t%.3f\t%.3f\t%s" % (chunk[0].chrom, chunk[0].start,
                                          chunk[-1].end, frv, iafi,
                                          chunk[0].transcript)
