import string
import fileinput
import sys
import re

class Record1(object):
	def __init__(self, fields):
		self.chr = fields[0]
		self.start = fields[1]
		self.end = fields[2]
		self.uniqid = fields[3].rstrip(";")
		self.gene = fields[4].rstrip(";").strip("\"")
		self.ref = fields[5]
		self.alt = fields[6]
		self.info = fields[7]