#!/usr/bin/python

import random
import sys

if len(sys.argv) < 8:
    print "usage:"
    print sys.argv[0] + " <seed> <#servers> <#services> <server spec file>" +\
        " <service spec file>"
    sys.exit()

seed         = int(sys.argv[1])
numservers   = int(sys.argv[2])
numservices  = int(sys.argv[3])
server_fname = sys.argv[4]
service_fname = sys.argv[5]

random.seed(seed)

print 2
print numservers
print numservices

serverf = open(server_fname, 'r')
for i in range(numservers):
    print serverf.readline()
serverf.close()

servicef = open(service_fname, 'r')
for i in range(numservices):
    print servicef.readline()
servicef.close()
