#!/usr/bin/python
import sys
with open(sys.argv[1],"r") as f:
	l = map(float,f.read().split())
	print "%.16lf %.16lf %.16lf" % (min(l),sum(l)/len(l),max(l))
