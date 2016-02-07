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

def evaluate(iterable, vcf_path, is_pathogenic=lambda v:
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

def upton(iterable, truth_set, cutoff=1e-3, vmin=1/(3*60706.)):

    def genchunks():
        nsmall = 0
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
                continue
            yield chunk, score
        print >>sys.stderr, nsmall, "removed for being too short"
        print >>sys.stderr, i, "total chunks"

    # NOTE: these are for humvar only. not neede for clinvar.
    def is_pathogenic(d):
        return d['class'] == "deleterious"
    def not_pathogenic(d):
        return d['class'] == "neutral"

    if "humvar" in truth_set :
        res = evaluate(genchunks(), truth_set, is_pathogenic=is_pathogenic,
            not_pathogenic=not_pathogenic)
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
