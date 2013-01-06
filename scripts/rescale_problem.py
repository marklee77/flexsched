#!/usr/bin/python

import sys

prob_fname = sys.argv[1]
out_fname = sys.argv[2]
slacks = [float(x) for x in sys.argv[3].split()]
loads = [float(x) for x in sys.argv[4].split()]

probf = open(prob_fname, 'r')

num_resources = int(probf.readline())
num_servers = int(probf.readline())
num_services = int(probf.readline())

resource_totals = []
requirement_totals = []
need_totals = []
for i in range(num_resources):
    resource_totals.append(0.0)
    requirement_totals.append(0.0)
    need_totals.append(0.0)

for i in range(num_servers):
    line = probf.readline()
    fields = line.split()
    for j in range(num_resources):
        resource_value = float(fields[num_resources + j])
        resource_totals[j] += resource_value

for i in range(num_services):
    line = probf.readline()
    fields = line.split()
    for j in range(num_resources):
        requirement_totals[j] += float(fields[2*num_resources + j])
        need_totals[j] += float(fields[3*num_resources + j])

probf.close()

print "resource totals:", resource_totals
print "requirement totals:", requirement_totals
print "need totals:", need_totals

requirement_factors = []
need_factors = []
for j in range(num_resources):

    if (requirement_totals[j] > 0.0):
        requirement_factors.append((1.0 - slacks[j]) * resource_totals[j] / 
            requirement_totals[j])
    else:
        requirement_factors.append(0.0)

    if (need_totals[j] > 0.0):
        # FIXME: my brain is being a little melty right now
        if (requirement_totals[j] > 0.0):
            need_factors.append(loads[j] * slacks[j] * resource_totals[j] / need_totals[j])
        else:
            need_factors.append(loads[j] * resource_totals[j] / need_totals[j])
    else:
        need_factors.append(0.0)

print "requirement factors:",
for j in range(num_resources):
    print requirement_factors[j],
print

print "need factors:",
for j in range(num_resources):
    print need_factors[j],
print

probf = open(prob_fname, 'r')

# just want to get to the good stuff again...
num_resources = int(probf.readline())
num_servers = int(probf.readline())
num_services = int(probf.readline())

outf = open(out_fname, 'w')

outf.write("%d\n" % (num_resources))
outf.write("%d\n" % (num_servers))
outf.write("%d\n" % (num_services))

for j in range(num_servers):
    outf.write("%s" % (probf.readline()))

for line in probf.readlines():
    fields = [float(x) for x in line.split()]
    for j in range(num_resources):
        outf.write("%5.3f\t" % (min(fields[j], 
            requirement_factors[j] * fields[2*num_resources + j])))
    for j in range(num_resources):
        outf.write("%5.3f\t" % min(fields[num_resources + j], 
            need_factors[j] * fields[3*num_resources + j]))
    for j in range(num_resources):
        outf.write("%5.3f\t" % (requirement_factors[j] * 
            fields[2*num_resources + j]))
    for j in range(num_resources):
        outf.write("%5.3f\t" % (need_factors[j] * fields[3*num_resources + j]))
    outf.write("\n")

probf.close()
outf.close()

