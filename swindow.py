import toolshed as ts
from collections import namedtuple, defaultdict, Counter
from operator import attrgetter

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

from math import log10

def IAFI_inline(intervals, n_samples):
    total_region_len = 0
    min_af = 1.0 / (2 * n_samples + 1)
    val = 0

    if hasattr(intervals, "start"):
        intervals = [intervals]

    for interval in intervals:
        region_len = interval.end - interval.start
        afs = map(float, (x if x != "." else min_af for x in interval.mafs.split(",")))
        afs.extend([min_af] * (region_len - len(afs)))
        
        # NOTE that sometimes we have more AFs than we have bases in the region. need to figure out why.
        # for now, we just add the len(afs) to the toal region len
        assert len(afs) >= region_len, (len(afs), region_len, interval.start, interval.end)
        val += sum(1.0/af for af in afs)
        total_region_len += len(afs)
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

def constraint(intervals, maf_cutoff):

    import re
    dn, ds, na = 0, 0, 0
    values = defaultdict(list)
    for iv in intervals:
        if iv.type is None:
            continue
        dnds = re.split(',|\|',iv.type)
        dn += dnds.count('dn')
        ds += dnds.count('ds')
        na += dnds.count('na')
    return float(dn) / (ds or 1)

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
    gene=set()
    n_domain_bases, n_gene_bases = 0, 0
    for iv in intervals:
        gene.add(iv.gene)
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
    return p.two_tail, (float(n_domain_variants) / (n_domain_bases or 1)) / denom, tbl, gene

def overlaps(a, b):
    return a[0] < b.end and a[1] > b.start

def evaldoms(iterable, vcf_path, is_pathogenic=lambda v:
                                                v.INFO.get("CLNSIG", 0) in "45"):
    """
    given a some chunks with a metric applied, do we see a difference in
    the values between pathogenic and non pathogenic variants?
    """
    from cyvcf2 import VCF

    iterable = list(iterable)

    tbl = {True: [], False: []}

    by_chrom = defaultdict(list)
    for it in iterable:
        by_chrom[it[0][0].chrom].append(it)


    for v in VCF(vcf_path):
        patho = is_pathogenic(v)

        cvars = by_chrom[v.CHROM]
        def chunkse(chunk):
            return chunk[0].start, chunk[-1].end

        vals = [it[1] for it in cvars if overlaps(chunkse(it[0]), v)]
        tbl[patho].extend(vals)

    return tbl


def runcontingent(path):
    from entropy import entropy
    import toolshed as ts
    it = ts.reader(path)
    iterable = (Interval(**iv) for iv in it)
    values = defaultdict(list)
    genes = set()
    by_transcript = defaultdict(list)
    by_domain = defaultdict(list)
    for iv in iterable:
        by_domain[iv.domain].append(iv)
        by_transcript[iv.transcript].append(iv)

    for domain, ivs in by_domain.items():
        if len(ivs) < 2: continue
        if sum(iv.mafs.count(',') for iv in ivs) < 3: continue
        if domain == ".": continue
        intervals = ivs[:]
        for iv in ivs:
            intervals.extend(by_transcript[iv.transcript])
        intervals = set(intervals)
        if len(intervals) < 3: continue
        pval, ratio, tbl, gene = contingent(intervals, domain, nodoms_only=False)
        ent = entropy(intervals)
        values['domain'].append(domain)
        values['pval'].append(pval)
        values['ent'].append(ent)
        values['tbl'].append(tbl)
        values['ratio'].append(ratio)
        values['num_intervals'].append(len(intervals))
        values['num_domains'].append(len(ivs))
        [genes.add(x) for x in gene]
        values['genes'].append(",".join(genes))
        genes=set()
    return values['domain'],values['pval'],values['ent'],values['tbl'],values['ratio'],values['num_intervals'],values['num_domains'],values['genes']

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

def metrics(trues, falses):
    from sklearn import metrics
    import matplotlib
    import seaborn as sns
    from matplotlib import pyplot as plt
    truth = ([1] * len(trues)) + ([0] * len(falses))

    obs = trues + falses
    dmetrics = {
     'auc': metrics.roc_auc_score(truth, obs),
     'precision': metrics.average_precision_score(truth, obs)
    }

    prec, rec, thresh = metrics.precision_recall_curve(truth, obs)
    fig, axes = plt.subplots(2)
    axes[0].plot(rec, prec)
    axes[0].set_xlabel("recall")
    axes[0].set_ylabel("precision")

    fpr, tpr, thresh = metrics.roc_curve(truth, obs)
    axes[1].plot(fpr, tpr, label="AUC: %.2f" % dmetrics['auc'])
    axes[1].set_xlabel('1 - specificity (FPR)')
    axes[1].set_ylabel('sensitivity (TPR)')
    axes[1].plot([0, 1], [0, 1], ls='--')
    axes[1].legend(loc="upper left")

    plt.show()

    return dmetrics

def tfloat(n):
    try:
        return float(n)
    except ValueError:
        return min(map(float, n.split("|")))


class Interval(object):
    def __init__(self, **entries):
        entries['start'] = int(entries['start'])
        entries['end'] = int(entries['end'])
        entries['chrom'] = entries['chr']
        self.__dict__.update(entries)

    def __repr__(self):
        return "Interval('%s@%s:%d-%d')" % (self.autoregs, self.chr, self.start, self.end)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __hash__(self):
        return hash(self.autoregs)

    @property
    def istarts(self):
        return map(int, self.starts.split(","))

    @property
    def iends(self):
        return map(int, self.ends.split(","))

    @property
    def positions(self):
        if self.pos == ".": return []
        return map(int, self.pos.split(","))

    @property
    def fmafs(self):
        if self.mafs == ".": return []
        return map(tfloat, self.mafs.split(","))

    @property
    def ftypes(self):
        if self.impacts == ".": return [] #impacts and types are switched, should fix later
        return self.impacts.split(",")

    def split(self, include_empties=False):
        posns, mafs, types = [x - 1 for x in self.positions], self.fmafs, self.ftypes
        if include_empties:
            starts, ends = self.istarts, self.iends
            for s, e in zip(starts, ends):
                for ip in range(s, e):
                    if not ip in posns:
                        posns.append(ip)
                        mafs.append(0.0)
                        types.append(None)

        pms = sorted(zip(posns, mafs, types))

        for p, maf, dntype in pms:
            I = Interval(**dict(self.__dict__.items()))
            I.start = p
            I.end = p + 1
            I.__dict__['mafs'] = str(maf)
            I.__dict__['pos'] = str(p + 1)
            I.__dict__['dnds'] = dntype
            I.aaf = maf
            yield I

class JimFile(object):
    def __init__(self, path, include_empties=True):
        self.path = path
        self.include_empties = include_empties

    def __iter__(self):
        cache = []
        for d in ts.reader(self.path):
            d = Interval(**d)
            # keep appending to the cache until we reach a different transcript
            # because nodoms occur in a different line from the domains but we
            # want everything to come out in order.
            if len(cache) > 0 and d.transcript != cache[0].transcript:
                for iv in sorted(cache, key=attrgetter('start')):
                    yield iv
                cache = []
            for iv in d.split(self.include_empties):
                cache.append(iv)
        for iv in cache:
            yield iv


def example():
    import toolshed as ts
    from collections import namedtuple

    it = ts.reader('/uufs/chpc.utah.edu/common/home/u6000294/lustre/u6000294/pmodel/y.sort.bed.gz')
    iterable = (Interval(**iv) for iv in it)
    for gene, val in slider(iterable, size_grouper(1), FRV_inline, maf_cutoff=0.005):
        print "%s\t%.3f\t%.3f" % (gene[0].autoregs, val, IAFI_inline(gene, 65000))

def domlimit(domain, pval, ent):
    label = []
    x = []
    y = []
    for i in range(0,len(domain)):
        if ent[i] < 0.4 and ent[i] > 0:
            label.append(domain[i])
            x.append(ent[i])
            y.append(pval[i])
        if pval[i] < 10 and pval[i] > 9.6 and ent[i] > 0.7:
            label.append(domain[i])
            x.append(ent[i])
            y.append(pval[i])
        if pval[i] < 0.05 and pval[i] > 0 and ent[i] > 0.8032 and ent[i] < 0.8035:
            label.append(domain[i])
            x.append(ent[i])
            y.append(pval[i])
    return label,x,y

def example2():
    domain, pval, ent, tbl, ratio, num_intervals, num_domains, genes = runcontingent(sys.argv[1]) #'/uufs/chpc.utah.edu/common/home/u6000294/lustre/u6000294/pmodel/y.sort.bed.gz'
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import statsmodels.stats.multitest as smm
    from math import log
    adj_p = pval
    #rej,adj_p=smm.multipletests(pval,method='bonferroni')[:2]
    adj_p = [-log(y,10) for y in adj_p]
    labels, x, y = domlimit(domain,adj_p,ent)
    print "domain\tgenes\tburden_pval\ttable\tratio\tn_intervals\tn_domain\tn_genes\tentropy"
    gd = []
    s = []
    for i in range(0,len(domain)):
        num_genes = len(genes[i].split(","))
        dom_muts = float(tbl[i].split(",")[1].split(":")[1].split("/")[0])
        print "%s\t%s\t%.4g\t%s\t%.2f\t%d\t%d\t%d\t%.4f" % (domain[i], genes[i], adj_p[i], tbl[i], ratio[i], num_intervals[i], num_domains[i], num_genes, ent[i])
        gd.append(log(num_genes,2))
        s.append(log(dom_muts,10)*10)

    #matplotlib.use('Agg')
    sc = plt.scatter(ent, adj_p, c = gd, s=s, edgecolors='none', cmap=cm.spectral)
    plt.xlim((0,1.05))
    plt.ylim((-1,11))
    plt.xlabel('Normalized Shannon entropy score')
    plt.ylabel('Mutational burden (-log10 p-value)')
    cb = plt.colorbar(sc, shrink = 0.5)
    cb.set_label("log2 (number of genes in domain family)")
    for label, x, y in zip(labels, x, y):
        plt.annotate(
            label, xy = (x,y), textcoords = 'offset points', xytext = (-20, 15),
            alpha = 1, rotation = 45, position = (x,y), size = 'x-small',
            weight = 'semibold')
    l1 = plt.scatter([],[], s=log(10,10)*8, edgecolors='none')
    l2 = plt.scatter([],[], s=log(100,10)*8, edgecolors='none')
    l3 = plt.scatter([],[], s=log(1000,10)*8, edgecolors='none')
    l4 = plt.scatter([],[], s=log(10000,10)*8, edgecolors='none')

    labels = ["10", "100", "1000", "10000"]

    leg = plt.legend([l1, l2, l3, l4], labels, ncol=1, frameon=False, fontsize=9,
    handlelength=2, loc='center', bbox_to_anchor=(1.1,0.9), borderpad = 2,
    handletextpad=1, title='number of mutations in domain', scatterpoints = 1)
    plt.setp(leg.get_title(), fontsize = 'small')
    plt.show()

def example3():
    import toolshed as ts
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    from scipy.stats import mannwhitneyu as mw

    it = ts.reader(sys.argv[1]) #'/scratch/ucgd/serial/quinlan_lab/data/u1021864/regionsmafsdnds.bed.gz'
    iterable = (Interval(**iv) for iv in it)

    results = defaultdict(list)
    maf_cutoff = 0.005

    for iv in windower(iterable, size_grouper(1)):
        results['constraint'].append(constraint(iv, maf_cutoff=maf_cutoff))
        results['iafi'].append(IAFI_inline(iv, n_samples=65000))
        results['frv'].append(FRV_inline(iv, maf_cutoff=maf_cutoff))
        # TODO: jim add a lot more metrics here... e.g.:
        # results['supermetric'].append(supermet(iv))

    fig, axes = plt.subplots(2)
    for metric in results: 
        print metric
        counts = evaldoms(results[metric], sys.argv[2]) # /uufs/chpc.utah.edu/common/home/u6000771/Projects/gemini_install/data/gemini/data/clinvar_20150305.tidy.vcf.gz
        axes[0].hist(counts[True])
        axes[0].set_xlabel("pathogenic")
        axes[1].hist(counts[False])
        axes[1].set_xlabel("not-pathogenic")
        plt.savefig(metric)
        plt.close()
        print metrics(counts[True],counts[False])
        plt.close()

if __name__ == "__main__":
    import doctest
    import sys
    print >>sys.stderr, (doctest.testmod())

    example3()
