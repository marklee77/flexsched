#!/usr/bin/python

import math
import sys

prob_fname = sys.argv[1]
probf = open(prob_fname, 'r')

probf.readline()
numhosts = int(probf.readline())
probf.readline()

machines = dict()

for i in range(numhosts):
    fields = probf.readline().split()
    cpuval = int(20.0 * float(fields[2]))
    memval = int(20.0 * float(fields[3]))
    machines[(cpuval, memval)] = i

probf.close()

print len(machines)
