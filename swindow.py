try:
    import matplotlib
    matplotlib.use("Agg")
except ImportError:
    pass
import toolshed as ts
from collections import namedtuple, defaultdict, Counter
from operator import attrgetter
import re
from pyfaidx import Fasta
from precision import remove_trailing_zeros as rtz
from argparse import ArgumentParser
from math import log, e
import scipy.stats as ss
from itertools import izip as zip


interval = namedtuple('interval', ['chrom', 'start', 'end'])
patt = re.compile(',|\|')

def size_grouper(n):
    return lambda grp, inext: len(grp) >= n

def byregiondist(grp, inext):
    """ group by region and split at gaps > 50 bases"""
    return inext.autoregs != grp[0].autoregs \
            or inext.start - grp[-1].end > 50

def bytranscriptdist(grp, inext):
    """ group by gene, split at gaps > 100 bases """
    return inext.transcript != grp[0].transcript \
            or inext.start - grp[-1].end > 100

def smallchunk(grp, inext, regionsize=15):
    """ group by chunk, input size, default is 50 """
    return len(grp) >= regionsize or inext.transcript != grp[0].transcript \
        or inext.start - grp[-1].end > 40

def frange(start, stop, step):
    r = start
    while (stop-r) > 1e-05: # allows a 1e-05 margin of error
        r += step # I put step before, to have an right closed range, instead of a left closed range
        yield float("%g" % r)


def windower(iterable, grouper=size_grouper(20), chunksize = ""):
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

    if chunksize == "":
        chunk = [next(iterable)]
        for iv in iterable:
            if grouper(chunk, iv):
                yield chunk
                chunk = [iv]
            else:
                chunk.append(iv)
        yield chunk
    else:
        chunk = [next(iterable)]
        for iv in iterable:
            if grouper(chunk, iv, chunksize):
                yield chunk
                chunk = [iv]
            else:
                chunk.append(iv)
        yield chunk

def RVIS(intervals, maf_cutoff = 1e-03, patt = patt):
    ct, l = 0.0, 0.0
    for iv in intervals:
        l += iv.end - iv.start
        if float(iv.mafs) >= maf_cutoff:
            dnds = patt.split(iv.type)
            ct += dnds.count('dn')
    return ct / l

def get_conservation(r):
        """
        Calculate the conservation for the interval based on the
        GERP scores weighted by the number of base pairs of overlap
        """
        start = float(r['start']); end = float(r['end']); overlap = float(r['overlap']); gerp = float(r['gerp'])
        overlap_fractions = (end - start) / overlap
        w_mu = gerp / overlap_fractions

        return w_mu

def baseline(intervals, maf_cutoff = 1e-05, exclude = None, comparison = "le", patt = patt):
    import operator
    def get_truth(inp, compare, cut):
        return compare(inp, cut)
    ct, l = 0.0, 0.0
    if comparison == "lt":
        comp = operator.lt
    elif comparison == "gt":
        comp = operator.gt
    elif comparison == "le":
        comp = operator.le
    elif comparison == "ge":
        comp = operator.ge
    for iv in intervals:
        l += iv.end - iv.start
        if get_truth(float(iv.mafs), comp, maf_cutoff):
            dnds = patt.split(iv.type)
            if exclude != None:
                if not dnds.count(exclude):
                    ct += 1
            else:
                ct += 1
    return intervals[0].chrom, intervals[0].start, intervals[-1].end, ct, l

def upton(base, baserate):
    obs = base[3] / base[4]
    exp = baserate
    if base[3]/base[4] < exp:
        upton = (base[0], base[1], base[2], 1)
    else:
        upton = (base[0], base[1], base[2], ss.binom_test(base[3], base[4], exp))
    return upton

def CpG(intervals, genes):
    n, l = 0.0, 0.0
    for iv in intervals:
        start = int(iv.start)
        end = int(iv.end)
        chrom = iv.chrom
        seq = genes[chrom][start:end].seq
        n += seq.count('CG')
        l += end - start
    return n*2.0/l

def IAFI(intervals, n_samples):
    minaf = 1.0 / (2 * n_samples + 1)
    val = sum(1.0/max(iv.aaf, minaf) for iv in intervals)
    return log(val / len(intervals), 10)

def IAFI_inline(intervals, n_samples, patt = patt):
    total_region_len = 0
    min_af = 1.0 / (2 * n_samples + 1)
    val = 0

    if hasattr(intervals, "start"):
        intervals = [intervals]

    for interval in intervals:
        region_len = interval.end - interval.start
        afs = map(float, (x if x != "." else min_af for x in patt.split(interval.mafs)))
        afs.extend([min_af] * (region_len - len(afs)))

        # NOTE that sometimes we have more AFs than we have bases in the region. need to figure out why.
        # for now, we just add the len(afs) to the total region len
        assert len(afs) >= region_len, (len(afs), region_len, interval.start, interval.end)
        val += sum(1.0/af for af in afs)
        total_region_len += len(afs)
    return val / total_region_len

def FRV(intervals, maf_cutoff): # default cutoff is 0.005
    return sum(1.0 for iv in intervals if iv.aaf <= maf_cutoff) / len(intervals)

def FRV_inline(intervals, maf_cutoff, patt = patt):

    n, s = 0, 0
    for interval in intervals:
        afs = map(float, (x if x != "." else 0 for x in patt.split(interval.mafs)))
        n += len(afs)
        s += sum(1.0 for af in afs if af <= maf_cutoff)
    return s / n

def dnds_density(intervals, maf_cutoff, patt = patt):
    dn, ds, na, l = 0.0, 0.0, 0.0, 0.0
    assert (x in set(['dn','ds','na','.']) for x in patt.split(intervals[0].type))
    for iv in intervals:
        dnds = patt.split(iv.type)
        dn += dnds.count('dn')
        ds += dnds.count('ds')
        na += dnds.count('na')
        l += iv.end - iv.start
    return float(dn)/l, float(dn) / (ds or 1), float(dn+ds+na) / l

def constraint(intervals, maf_cutoff, genes, upton):
    #cpg = CpG(intervals, genes)
    #upton = upton[3]
    nons, dnds, density = dnds_density(intervals, maf_cutoff)
    #base = baseline(intervals, maf_cutoff)[3]
    #iafi = IAFI_inline(intervals, n_samples=61000)
    #dn_density = count_nons(intervals)
    constraint = density*nons#upton#*float(1-cpg)

    return constraint

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

def uptest():
    b = "/scratch/ucgd/lustre/u1021864/serial/y.sort.bed.gz"
    it = ts.reader(b)
    iterable = (Interval(**iv) for iv in it)

def rvistest():
    vcf_path = "/scratch/ucgd/lustre/u1021864/serial/clinvar-anno.vcf.gz"
    bed = "rvis.bed"

    def genregions():
        for d in ts.reader("rvis.bed"):
            score = float(d['pct'])
            chunk = [interval(d['chrom'], int(d['start']), int(d['end']))]
            yield chunk, -score

    res = evaldoms(genregions(), vcf_path)
    print metrics(res[True], res[False], "x.auc.png")

def evaldoms(iterable, vcf_path, is_pathogenic=lambda v:
                                any(x == "5" and v.INFO.get("max_aaf_all", -1) < 0.001 for x in re.split(patt,v.INFO.get("CLNSIG"))),
                                not_pathogenic=lambda v: any(x == "2" for x in re.split(patt,v.INFO.get("CLNSIG")))):
    """
    given a some chunks with a metric applied, do we see a difference in
    the values between pathogenic and non pathogenic variants?
    """
    from cyvcf2 import VCF
    from interlap import InterLap

    tbl = {True: [], False: []}

    tree = {True: defaultdict(InterLap), False: defaultdict(InterLap)}
    n, p = 0, 0
    if vcf_path.endswith((".vcf", ".vcf.gz")):
        for v in VCF(vcf_path):
            path = is_pathogenic(v)
            nopath = not_pathogenic(v)
            if path == nopath: continue
            if path:
                p += 1
            else:
                n += 1
            # is it pathogenic
            tree[path][v.CHROM].add((v.start, v.end))
    else:
        for d in ts.reader(vcf_path):
            path = is_pathogenic(d)
            nopath = not_pathogenic(d)
            if path == nopath: continue
            if path:
                p += 1
            else:
                n += 1

            chrom = d.get('chrom', d['chr']).replace('chr', '')
            start, end = int(d['start']), int(d['end'])
            tree[path][chrom].add((start, end))

    print >>sys.stderr, "pathogenic variants: %d non: %d" % (p, n)
    counts = {True: 0, False: 0, "missing": 0}
    for reg in iterable:
        chrom = reg[0][0].chrom
        start, end = reg[0][0].start, reg[0][-1].end
        patho = len(list(tree[True][chrom].find((start, end)))) != 0
        nonpatho = len(list(tree[False][chrom].find((start, end)))) != 0
        if not (patho or nonpatho):
            counts["missing"] += 1

        if patho == nonpatho:
            continue
        tbl[patho].append(reg[1])
        counts[patho] += 1
    print >>sys.stderr, "region counts:", counts
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
    [1.4313637641589871, 5.698970004336018]

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

def metrics(trues, falses, figname=None, cutoff = None):
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
    fig.tight_layout()
    axes[0].plot(rec, prec, label = "precision: %.2f" % dmetrics['precision'],
            marker='.', lw=0.5)
    axes[0].legend(frameon = False, loc = 'upper left')
    axes[0].set_xlabel("recall")
    axes[0].set_ylabel("precision")
    props = dict(boxstyle = 'round', facecolor = 'whitesmoke', alpha = 0.5)
    if cutoff is not None:
        axes[0].text(.85, .8, "CpG frac:\n" + cutoff.replace("-"," - "), transform = axes[0].transAxes, bbox = props)

    fpr, tpr, thresh = metrics.roc_curve(truth, obs)
    axes[1].plot(fpr, tpr, label = "AUC: %.2f" % dmetrics['auc'], marker='.',
            lw=0.5)
    axes[1].set_xlabel('1 - specificity (FPR)')
    axes[1].set_ylabel('sensitivity (TPR)')
    axes[1].plot([0, 1], [0, 1], ls='--')
    axes[1].legend(loc = "upper left")
    axes[1].text(0.55, 0.2, "%d true+ : %d true- regions" % (len(trues),
        len(falses)))

    plt.savefig(figname, bbox_inches = 'tight')
    plt.close()

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
        if self.type == ".": return [] #impacts and types are switched in ryans original file, fixed in mine and in code but need to make a note
        return self.type.split(",")

    def split(self, include_empties=False):
        posns, mafs, types = [x - 1 for x in self.positions], self.fmafs, self.ftypes
        if include_empties:
            starts, ends = self.istarts, self.iends
            sposns = set(posns)
            for s, e in zip(starts, ends):
                diff = set(range(s, e)) - sposns
                posns.extend(sorted(diff))
                mafs.extend([0.0] * len(diff))
                types.extend(['.'] * len(diff))

        pms = sorted(zip(posns, mafs, types))

        for p, maf, dntype in pms:
            I = Interval(**dict(self.__dict__.items()))
            I.start = p
            I.end = p + 1
            I.__dict__['mafs'] = str(maf)
            I.__dict__['pos'] = str(p + 1)
            I.__dict__['type'] = dntype
            I.aaf = maf
            yield I

    def __str__(self):
        return "%s\t%d\t%d\t%s\t%s\t%s" % (
                self.chrom, self.start, self.end, self.mafs, self.transcript, self.domain)

class JimFile(object):
    def __init__(self, path, regions_excluded=None, include_empties=True):
        self.path = path
        self.include_empties = include_empties
        self.re = regions_excluded

    def iterchunks(self):

        interval = namedtuple('interval', ['chrom', 'start', 'end', 'mafs'])
        for d in ts.reader(self.path):
            starts = map(int, d['starts'].split(","))
            ends = map(int, d['ends'].split(","))
            fmafs = map(tfloat, d['mafs'].split(",")) if d['mafs'] != "." else []
            if len(starts) == len(ends) == 1:
                fmafs.extend([0.0] * end)
                yield Interval(start=starts, end=ends, mafs=fmafs)
            print d['mafs']
            print starts,ends,fmafs


    def __iter__(self):
        cache = []
        for d in ts.reader(self.path):
            d = Interval(**d)
            if self.re == "nodoms":
                if "NoDom" not in d.autoregs:
                    continue
            elif self.re == "domains":
                if "NoDom" in d.autoregs:
                    continue
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
    domain, pval, ent, tbl, ratio, num_intervals, num_domains, genes = runcontingent(args.input) #'/uufs/chpc.utah.edu/common/home/u6000294/lustre/u6000294/pmodel/y.sort.bed.gz'
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import statsmodels.stats.multitest as smm
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

def clinvar(v):
    return any(x in "45" for x in re.split(patt,v.INFO.get("CLNSIG")))

def pli(v):
    return float(v['pLI']) < 0.9

def tester():

    print "#chrom	start	end	maf	transcript	domain"
    input = "/scratch/ucgd/lustre/u1021864/serial/y.sort.bed.gz"
    iterator = JimFile(input, include_empties=True)
    for p in iterator:
        print str(p)
def gerprunner():
    
    import pyBigWig

    b = pyBigWig.open("/scratch/ucgd/lustre/u1021864/serial/hg19.gerp.bw")
   # x = list(range(1,23)); x.append("X"), x.append("Y")

    input = sys.argv[1]
    iterator = JimFile(input)
    iterable = windower(iterator, smallchunk(regionsize=1))
    cutoff = 1e-3

    def genchunks():
        nsmall = 0
        for i, chunk in enumerate(iterable):
            #if len(chunk) < 5:
            #    continue
            score = b.stats("chr"+chunk[0].chrom, chunk[0].start, chunk[-1].end)
            yield chunk, score[0]
            if i % 100000 == 0:
                print i, chunk[0].chrom, chunk[0].start, score
        print >>sys.stderr, nsmall, "removed for being too short"
        print >>sys.stderr, i, "total chunks"

    vcf_path = "/scratch/ucgd/lustre/u1021864/serial/clinvar-anno.vcf.gz"
    res = evaldoms(genchunks(), vcf_path)
    print metrics(res[True], res[False], "gerp.auc.png")

def uptonrunner():

    input = "/scratch/ucgd/lustre/u1021864/serial/y.sort.bed.gz"
    iterator = JimFile(input)
    iterable = windower(iterator, smallchunk)
    cutoff = 1e-3

    def genchunks():
        nsmall = 0
        for i, chunk in enumerate(iterable):
            if i % 100000 == 0:
                print i, chunk[0].chrom, chunk[0].start
            if len(chunk) < 5:
                continue
            mafs = (float(x.mafs) for x in chunk)
            score = sum(1.0 - m for m in mafs if m < cutoff) / float(len(chunk))
            if score == 1:
                nsmall += 1
                continue
            yield chunk, score
        print >>sys.stderr, nsmall, "removed for being too short"
        print >>sys.stderr, i, "total chunks"

    # NOTE: these are for humvar only. not neede for clinvar.
    def is_pathogenic(d):
        return d['class'] == "deleterious"
    def not_pathogenic(d):
        return d['class'] == "neutral"

    eval_path = "/scratch/ucgd/lustre/u1021864/serial/clinvar-anno.vcf.gz"
    eval_path = "/uufs/chpc.utah.edu/common/home/u1007787/pmodel/humvar.both.bed"
    res = evaldoms(genchunks(), eval_path, is_pathogenic=is_pathogenic,
            not_pathogenic=not_pathogenic)
    print metrics(res[True], res[False], "upton.auc.png")

def example3():
    import toolshed as ts
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import seaborn as sns
    from scipy.stats import mannwhitneyu as mw
    import numpy as np

    iterator = JimFile(args.input, args.regions)
    #it = ts.reader(args.input) #'/scratch/ucgd/serial/quinlan_lab/data/u1021864/regionsmafsdnds.bed.gz'
    #iterable = (Interval(**iv) for iv in it)

    results = defaultdict(lambda : defaultdict(list))
    ms = defaultdict(list)
    ff = args.genome
    cpg_cutoff = {}
    maf_cutoff = float(args.maf) if args.maf else 1e-05
    start = 0
    end = .2
    step = .025
    j = start
    #for i in frange(start, end, step):
    #    cpg_cutoff[str(j)+"-"+str(i)] = (j, i)
    #    j = i
    #cpg_cutoff['0.2-1'] = (.2, 1)
    cpg_cutoff['0-1'] = (0, 1)

    base = []
    cons = []
    genes = None
    #genes = Fasta(ff)
    if args.regions == "chunks":
        regioner = smallchunk
        chunksize = args.regionsize
    if args.regions in ["domains", "nodoms", "all"]:
        regioner = byregiondist
        chunksize = ""
    if args.regions == "genes":
        regioner = bytranscriptdist
        chunksize = ""
    y = list(windower(iterator, regioner, chunksize))
    comparison = args.comparison
    if args.exclude:
        exclude = args.exclude
        ex = "ex" + args.exclude + "."
    else:
        exclude = None
        ex = ""
    cv = []
    if args.conservation:
        for r in ts.reader(args.conservation):
            v = get_conservation(r)
            cv.append(v)
    cpg=1
    if y:
        for iv in y: # iterable, size_grouper(1)
            #cpg = CpG(iv, genes = genes)
            b = baseline(iv, maf_cutoff = maf_cutoff, exclude = exclude, comparison = comparison, patt = patt)
            ms['baseline'].append((iv,b[3]/b[4],cpg))
            base.append(b)
    count = 0.0
    totlen = 0.0
    if base:
        for b in base:
            count += b[3]
            totlen += b[4]
        baserate = count/totlen
    for iv, b in zip(y, base):
        u = upton(b, baserate)
        c = constraint(iv, maf_cutoff = maf_cutoff, genes = genes, upton = u)
        r = RVIS(iv, maf_cutoff = 1e-3, patt = patt)
        ct = (iv,
               c,
               cpg)
        if c != 0:
            ms['nzconstraint'].append(ct)
        ms['constraint'].append(ct)
        ct = (iv,
                u,
                cpg)
        ms['upton'].append((ct[0],ct[1][3],ct[2]))
        ct = (iv,
                r,
                cpg)
        ms['rvis'].append((ct[0],ct[1],ct[2]))
        cons.append((u[0],u[1],u[2],c))
       # results['iafi'].append((iv, IAFI_inline(iv, n_samples=61000)))
       # results['frv'].append((iv, FRV_inline(iv, maf_cutoff=maf_cutoff)))
       # results['count_nons'].append((iv, count_nons(iv)))
        # TODO: jim add a lot more metrics here... e.g.:
    bedname = "."+ rtz(maf_cutoff) + "." + comparison + "." + args.regions + str(chunksize) + "." + ex
    f1 = open("constraint" + bedname + ".bed","w")
    f2 = open("baseline" + bedname + ".bed","w")
    for b,c in zip(base,cons):
        f1.write("\t".join(map(str,c))+"\n")
        f2.write("\t".join(map(str,b))+"\n")
    f1.close()
    f2.close()

    cutoffs = set()
    for cutoff in cpg_cutoff:
        co = str(cpg_cutoff[cutoff][0])+'-'+str(cpg_cutoff[cutoff][1])
        cutoffs.add(co)
        for metric in ms:
            for ct in ms[metric]:
                if ct[2] >= cpg_cutoff[cutoff][0] and ct[2] <= cpg_cutoff[cutoff][1]:
                    results[metric][co].append(ct)

    option = args.truetype
    trusrc = ""
    if option == "clinvar" or option == "c":
        func = clinvar
        trusrc = "clinvar"
    if option == "pli" or option == "p":
        func = pli
        trusrc = "pli"
    for metric in results:
        for cutoff in cutoffs:
            imgname = metric + "." + trusrc + "." + comparison + "." + args.regions + str(chunksize) + "." + ex + cutoff + "." + rtz(maf_cutoff)
            print metric, cutoff
            fig, axes = plt.subplots(2)
            fig.tight_layout()
            counts = evaldoms(results[metric][cutoff],
                    args.pathogenic, # forweb_cleaned_exac_r03_march16_z_data_pLI.txt from ExAC ftp or clinvar_20150305.tidy.vcf.gz from clinvar src
                    func)
            imin, imax = np.percentile(counts[True] + counts[False], [0.01, 99.99])
            axes[0].hist(counts[True], bins=80) #,label = cutoff)
            axes[0].set_xlabel("pathogenic")
            axes[0].set_xlim(imin, imax)
            props = dict(boxstyle = 'round', facecolor = 'whitesmoke', alpha = 0.5)
            axes[0].text(.875, .8, "CpG frac:\n" + cutoff.replace("-"," - "), transform = axes[0].transAxes, bbox = props)
            #axes[0].legend(loc = 1, frameon = True)
            axes[1].hist(counts[False], bins=80)
            axes[1].set_xlabel("not-pathogenic")
            axes[1].set_xlim(imin, imax)
            plt.show()
            plt.savefig(imgname + ".dist.png", bbox_inches = 'tight')
            print metrics(counts[True], counts[False], imgname + ".auc.png", cutoff = cutoff)
            print mw(counts[True], counts[False])
            del fig
            plt.close()

if __name__ == "__main__":
    import doctest
    import sys
    print >>sys.stderr, (doctest.testmod())
    #rvistest()

    if len(sys.argv) > 1 and sys.argv[1] == "test":
        sys.exit(tester())
    #uptonrunner()
    gerprunner()
    1/0
    #1/0

    parser = ArgumentParser()
    parser.add_argument("--input", "-i", help = "file with regions (doms and nodoms) defined", type = str)
    parser.add_argument("--genome", "-g", help = "fasta genome file for CpG calculations", type = str)
    parser.add_argument("--pathogenic", "-p", help = "file with truth set for pathogenicity comparison", type = str)
    parser.add_argument("--truetype", "-t", help = "truth set type, specifies what function to use, e.g., clinvar vs pLI", type = str)
    parser.add_argument("--maf", "-m", help = "maf cutoff for baseline/uptoN metric and other metrics that utilize maf cutoffs", type = float)
    parser.add_argument("--exclude", "-e", help = "what to exclude from calculations regarding upton metric", type = str)
    parser.add_argument("--comparison", "-c", help = "which comparison to use for maf cutoff, greater than/equal to (le), ge, lt, gt", type = str, default = "le")
    parser.add_argument("--regions", "-r", help = "select regional model, whether to go by file or pick region size, must use region size option if selected, choices are domains, nodoms, chunks, genes, all (regions)", type = str, default = "all")
    parser.add_argument("--regionsize", "-s", help = "select region size for analysis", type = int)
    parser.add_argument("--conservation", "-v", help = "GERP conservation file", type = str)

    args = parser.parse_args()

    example3()
