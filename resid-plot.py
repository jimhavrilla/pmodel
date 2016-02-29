import toolshed as ts
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import OLSInfluence
import scipy.stats as ss
from statsmodels.formula.api import ols
import pandas as pd


import csv
import sys
csv.field_size_limit(14365000)

cutoff = 0.3

def fn(scores, covs):
    c = np.array(covs)**2
    #c /= c.sum()
    return sum(sc * cov for sc, cov, oc in zip(scores, c, covs) if oc > cutoff)

sigma = []
cgs, zs, xs, ys, genes = [], [], [], [], []
try:
    for i, d in enumerate(ts.reader(1)):
        gerps = [float(x)+13 for x in d['gerp'].split(",")]
        genes.append((d['chrom'], str(d['start']), str(d['end']), d['gene'],
            d['exon'], str(len(gerps))))

        coverage = map(float, d['coverage'].split(","))
        score = fn(gerps, coverage)
        cgs.append(float(d['cg_content']))

        ys.append(score)
        xs.append(sum(c for c in coverage if c > cutoff))
        #sigma.append(np.std([sc * cov for sc, cov in zip(gerps, coverage)]) or 1)

        # just subtract scaled gerp from scaled coverage.
        zs.append(sum(cov - g / 25. for cov, g in zip(coverage, gerps) if cov >
            cutoff))
except:
    raise

X = pd.DataFrame({"CpG": cgs, "sum-coverage": xs})
#X = np.array(cgs)
results = sm.OLS(ys, X, hasconst=False).fit()

l = results.predict([np.zeros(X.shape[1] if len(X.shape) > 1 else 1), X.max(0)])

print >>sys.stderr, results.params

plt.plot(cgs, ys, marker='.', ls='none')
plt.plot([0, max(cgs)], l)
plt.title("split by missense")
plt.xlabel("sum(coverage)")
plt.ylabel("sum(cov * g for cov, g in zip(coverage, gerp))")
#plt.show()

resid = results.resid
resid = OLSInfluence(results).get_resid_studentized_external()

pctile = 100.0 * np.sort(resid).searchsorted(resid) / float(len(resid))

score_pctile = 100.0 * np.sort(xs).searchsorted(xs) / float(len(xs))
z_pctile = 100.0 * np.sort(zs).searchsorted(zs) / float(len(zs))

print "chrom\tstart\tend\tgene\texon\tn\tgerp_resid\tgerp_resid_pctile\tscore_pctile\tz_pctile"
for i, row in enumerate(genes):
    print "\t".join(list(row) + ["%.3f" % resid[i], "%.9f" % pctile[i],
                                 "%.9f" % score_pctile[i],
                                 "%.9f" % z_pctile[i],
                                 ])

#plt.close()
#fig = plt.figure(figsize=(12,8))
#fig = sm.graphics.plot_partregress_grid(results, fig=fig)
#plt.show()
