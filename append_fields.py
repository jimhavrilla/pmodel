import sys
import string
from toolshed import reader, header, nopen

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-f",
    dest="field_set",
    help="CSV of needed fields")

(options, args) = parser.parse_args()

if not options.field_set:
    parser.error('Field set not given')

field_set = options.field_set.split(',')

print '#' + '\t'.join(field_set)

for d in reader('-'):
    O = []
    for field in field_set:
        if field in d:
            O.append(d[field])
        else:
            O.append('.')
    print '\t'.join(O)
