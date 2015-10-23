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

s = set()
for d in reader('-',header='ordered') :
    for h in headers:
        print d[h],
    print
