import sys
from swindow import windower, JimFile, metrics
from collections import defaultdict
import re

patt = re.compile(',|\|')


def chunker(size=40):

    def fn(grp, inext):
        """ group by chunk, input size, default is 50 """
        return len(grp) > size or inext.transcript != grp[0].transcript \
            or inext.start - grp[-1].end > (min(size, 40))
    return fn

def is_pathogenic(v):
    max_aaf_all = v.INFO.get("max_aaf_all", -1.0)
    pos, neg = 0, 0
    for sig in patt.split(v.INFO['CLNSIG']):
        if sig == "5" and max_aaf_all < 0.001:
            pos += 1
        elif sig == "2":
            neg += 1
    if pos == neg: return None
    if pos > 0 and neg > 0: return None
    return pos > 0


def evaluate(iterable, vcf_path, is_pathogenic=is_pathogenic):
    """
    given a some chunks with a metric applied, do we see a difference in
    the values between pathogenic and non pathogenic variants?
    the is_pathogenic must return True for pathogenic, False for non-pathogenic
    and None for all other variants.
    """
    from cyvcf2 import VCF
    from interlap import InterLap

    tbl = {True: [], False: []}
    counts = {True: 0, False: 0}

    tree = {True: defaultdict(InterLap), False: defaultdict(InterLap)}

    if vcf_path.endswith((".vcf", ".vcf.gz")):
        for v in VCF(vcf_path):
            path = is_pathogenic(v)
            if path is None: continue
            counts[path] += 1
            # is it pathogenic
            tree[path][v.CHROM].add((v.start, v.end))
    else:
        import toolshed as ts
        for d in ts.reader(vcf_path):
            path = is_pathogenic(d)
            if path is None: continue
            counts[path] += 1

            chrom = d.get('chrom', d['chr']).replace('chr', '')
            start, end = int(d['start']), int(d['end'])
            tree[path][chrom].add((start, end))

    print >>sys.stderr, "pathogenic variants: %d non: %d" % (counts[True], counts[False])
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


def upton(iterable, truth_set, cutoff=1e-3, vmin=1/(3*60706.)):

    def genchunks():
        nsmall, ones = 0, 0
        for i, chunk in enumerate(iterable):
            if i % 100000 == 0:
                if i > 0:
                    print i, chunk[0].chrom, chunk[0].start
            if len(chunk) < 8:
                nsmall += 1
                continue
            mafs = (float(x.mafs) for x in chunk)
            score = sum((1.0 - m)**2.0 for m in mafs if m < cutoff) / float(len(chunk))
            if score == 1:
                ones += 1
                continue
            yield chunk, score
        sys.stderr.write("%d (%.2f%%) removed for being too short\n" % (nsmall,
                         100.0 * nsmall / float(i)))
        sys.stderr.write("%d (%.2f%%) removed for having score of 1\n" % (ones,
                         100.0 * ones /float(i)))
        print >>sys.stderr, i, "total chunks"

    # NOTE: this is for humvar only. not needed for clinvar.
    def humvar_pathogenic(d):
        return d['class'] == "deleterious"

    if "humvar" in truth_set:
        res = evaluate(genchunks(), truth_set, is_pathogenic=humvar_pathogenic)
    else:
        res = evaluate(genchunks(), truth_set)

    print metrics(res[True], res[False], "upton.auc.png")


def main(input, truth_set, aaf_cutoff, chunk_size):

    iterator = JimFile(input)
    iterable = windower(iterator, chunker(chunk_size))

    upton(iterable, truth_set, aaf_cutoff)


if __name__ == "__main__":

    import argparse
    p = argparse.ArgumentParser()

    p.add_argument("--input", "-i",
            default="/scratch/ucgd/lustre/u1021864/serial/y.sort.bed.gz")
    p.add_argument("--truth-set", "-t",
            default="/scratch/ucgd/lustre/u1021864/serial/clinvar-anno.vcf.gz")
    p.add_argument("--aaf-cutoff", type=float,
            default=5e-3)
    p.add_argument("--chunk-size", type=int,
            help="in # of bases.",
            default=20)

    a = p.parse_args()
    main(a.input, a.truth_set, a.aaf_cutoff, a.chunk_size)
