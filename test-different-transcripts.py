import re
import sys

patt = re.compile("ENST[^\s_]*")

for line in sys.stdin:
    assert len(set(patt.findall(line))) == 2
