#!/usr/bin/python

import sys

n = int(sys.argv[2])
ifilename = sys.argv[1]
ofilename = sys.argv[1] + '.data'

with open(ifilename,'r') as ifile:
	data = map(float,ifile.read().split())
	maxVal = max(data)
	minVal = min(data)
	
	freq = {}
	delta = (maxVal - minVal)/n
	for i in xrange(n - 1):
		left = minVal + i*delta
		right = minVal + (i+1)*delta if i < n - 1 else maxVal + 1
		freq[minVal + delta*(i + 0.5)] = len(filter(lambda x: left <= x < right,data))


with open(ofilename,'w') as ofile:
	for k in sorted(freq.keys()):
		ofile.write("%.16lf\t%d\n" % (k, freq[k]))

print ofilename
