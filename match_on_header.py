#!/usr/bin/env python
import sys
from toolshed import reader, header, nopen
from optparse import OptionParser

parser = OptionParser()

parser.add_option( "--headers",
    dest="headers",
    help="header name CSV")

(options, args) = parser.parse_args()

if not options.headers:
    parser.error('Headers not given')

headers = options.headers.split(',')

if len(headers) <= 1:
    sys.stderr.write("Must give more than one header name\n");
    exit(1)

for d in reader('-',header='ordered') :
    s = set()
    for h in headers:
        s.add(d[h])
    if len(s) == 1:
        print '\t'.join([x for x in s] + d.values())
