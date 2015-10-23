import sys

col = {}

for l in sys.stdin:
    A = l.rstrip().split('\t')
    for i in range(len(A)):
        if not i in col:
            col[i] = {}
        if not A[i] in col[i]:
            col[i][A[i]] = 0
        col[i][A[i]] += 1

for i in col:
    print len (col[i])
