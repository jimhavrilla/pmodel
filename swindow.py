from collections import namedtuple, defaultdict

interval = namedtuple('interval', ['chrom', 'start', 'end'])

def size_grouper(n):
    return lambda grp, inext: len(grp) >= n

def windower(iterable, grouper=size_grouper(20)):
    """
    windower takes an iterable of intervals and yields groups
    defined by grouper.

    grouper has signature (chunk, next_item) and returns true if a
    new chunk should be formed *without* the next_item

    >>> iterable = [interval('2', 22, 24), interval('2', 24, 25), interval('2', 33, 34)]
    >>> list(windower(iterable, lambda grp, iv: iv.start - grp[-1].end > 5))
    [[interval(chrom='2', start=22, end=24), interval(chrom='2', start=24, end=25)], [interval(chrom='2', start=33, end=34)]]
    >>> list(windower(iterable, lambda grp, iv: iv.start - grp[-1].end > 50))
    [[interval(chrom='2', start=22, end=24), interval(chrom='2', start=24, end=25), interval(chrom='2', start=33, end=34)]]

    """
    iterable = iter(iterable)

    chunk = [next(iterable)]
    for iv in iterable:
        if grouper(chunk, iv):
            yield chunk
            chunk = [iv]
        else:
            chunk.append(iv)
    yield chunk


def IAFI(intervals, n_samples):
    from math import log10
    minaf = 1.0 / (2 * n_samples + 1)
    val = sum(1.0/max(iv.aaf, minaf) for iv in intervals)
    return log10(val / len(intervals))

def IAFI_inline(intervals, n_samples):
    from math import log10
    total_region_len = 0
    min_af = 1.0 / (2 * n_samples + 1)
    val = 0

    for interval in intervals:
        region_len = interval.end - interval.start
        afs = map(float, (x if x != "." else min_af for x in interval.mafs.split(",")))
        afs.extend([min_af] * (region_len - len(afs)))
        assert len(afs) == region_len
        val += sum(1.0/af for af in afs)
        total_region_len += region_len
    return val / total_region_len

def FRV(intervals, maf_cutoff):
    return sum(1.0 for iv in intervals if iv.aaf <= maf_cutoff) / len(intervals)

def FRV_inline(intervals, maf_cutoff):

    n, s = 0, 0
    for interval in intervals:
        afs = map(float, (x if x != "." else 0 for x in interval.mafs.split(",")))
        n += len(afs)
        s += sum(1.0 for af in afs if af <= maf_cutoff)
    return s / n

def contingent(intervals, domain_name, nodoms_only=False):
    """
    intervals should be all intervals in all genes that contain the domain
    """
    import fisher

    n_domain_variants = sum(len(i.mafs.split(",")) for i in intervals if i.domain == domain_name)
    if nodoms_only:
        n_gene_variants = sum(len(i.mafs.split(",")) for i in intervals if i.domain == ".")
    else:
        n_gene_variants = sum(len(i.mafs.split(",")) for i in intervals if i.domain != domain_name)

    n_domain_bases, n_gene_bases = 0, 0
    for iv in intervals:
        starts = map(int, iv.starts.split(","))
        ends = map(int, iv.ends.split(","))
        l = sum(e - s for s, e in zip(starts, ends))
        assert all(e > s for s, e in zip(starts, ends)), domain_name
        if iv.domain == domain_name:
            n_domain_bases += l
        elif nodoms_only and iv.domain == ".":
            n_gene_bases += l
        elif not nodoms_only and iv.domain != domain_name:
            n_gene_bases += l
    tbl = "gene:%d/%d,dom:%d/%d" % (n_gene_variants, n_gene_bases, n_domain_variants, n_domain_bases)

    p = fisher.pvalue(n_gene_bases, n_gene_variants, n_domain_bases, n_domain_variants)

    denom = float(n_gene_variants) / (n_gene_bases or 1) or 1
    return p.two_tail, (float(n_domain_variants) / (n_domain_bases or 1)) / denom, tbl

def runcontingent(path):
    import toolshed as ts
    it = ts.reader(path)
    iterable = (Interval(**iv) for iv in it)

    by_transcript = defaultdict(list)
    by_domain = defaultdict(list)
    for iv in iterable:
        by_domain[iv.domain].append(iv)
        by_transcript[iv.transcript].append(iv)

    print "domain\tpval\ttable\tratio\tn_intervals\tn_domain"
    for domain, ivs in by_domain.items():
        if len(ivs) < 2: continue
        if domain == ".": continue
        intervals = ivs[:]
        for iv in ivs:
            intervals.extend(by_transcript[iv.transcript])
        intervals = set(intervals)
        if len(intervals) < 3: continue
        pval, ratio, tbl = contingent(intervals, domain, nodoms_only=False)
        print "%s\t%.4g\t%s\t%.2f\t%d\t%d" % (domain, pval, tbl, ratio, len(intervals), len(ivs))


def slider(iterable, grouper, metric, **kwargs):
    """
    >>> interval = namedtuple('interval', ['chrom', 'start', 'end', 'value'])
    >>> iterable = [interval('2', 22, 24, 0.25), interval('2', 24, 25, 0.75), interval('2', 33, 34, 0.25)]
    >>> vals = list(slider(iterable, size_grouper(2), lambda vals: sum(v.value for v in vals) / len(vals)))
    >>> len(vals)
    2
    >>> vals[0]
    ([interval(chrom='2', start=22, end=24, value=0.25), interval(chrom='2', start=24, end=25, value=0.75)], 0.5)

    >>> vals[1]
    ([interval(chrom='2', start=33, end=34, value=0.25)], 0.25)

    >>> interval = namedtuple('interval', ['chrom', 'start', 'end', 'aaf'])
    >>> iterable = [interval('2', 22, 24, 0.25), interval('2', 24, 25, 0.02), interval('2', 33, 34, 0.000002)]
    >>> [x[1] for x in slider(iterable, size_grouper(2), IAFI, n_samples=3333333)]
    [1.4313637641589874, 5.698970004336019]

    >>> [x[1] for x in slider(iterable, size_grouper(2), FRV, maf_cutoff=0.05)]
    [0.5, 1.0]

    >>> [x[1] for x in slider(iterable, size_grouper(2), FRV, maf_cutoff=0.0005)]
    [0, 1.0]

    >>> interval = namedtuple('interval', ['chrom', 'start', 'end', 'mafs'])
    >>> iterable = [interval('2', 22, 23, '0.05,0.002,0.1,0.005')]
    >>> [x[1] for x in slider(iterable, size_grouper(1), FRV_inline, maf_cutoff=0.05)]
    [0.75]

    """
    for chunk in windower(iterable, grouper):
        yield chunk, metric(chunk, **kwargs)


class Interval(object):
    def __init__(self, **entries):
        entries['start'] = int(entries['start'])
        entries['end'] = int(entries['end'])
        self.__dict__.update(entries)
    def __repr__(self):
        return "interval('%s@%s:%d-%d')" % (self.autoregs, self.chr, self.start, self.end)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __hash__(self):
        return hash(self.autoregs)

def example():
    import toolshed as ts
    from collections import namedtuple

    it = ts.reader('/scratch/ucgd/serial/quinlan_lab/data/u1021864/regionsmafsdnds.bed.gz')
    iterable = (Interval(**iv) for iv in it)
    for gene, val in slider(iterable, size_grouper(1), FRV_inline, maf_cutoff=0.005):
        print "%s\t%.3f\t%.3f" % (gene[0].autoregs, val, IAFI_inline(gene, 65000))

def example2():
    runcontingent('/scratch/ucgd/serial/quinlan_lab/data/u1021864/regionsmafsdnds.bed.gz')


if __name__ == "__main__":
    import doctest
    import sys
    print >>sys.stderr, (doctest.testmod())

    example2()
