#!/usr/bin/python

import random
import sys

if len(sys.argv) < 11:
    print "usage:"
    print sys.argv[0] + " <seed> <#rigid> <#fluid> <#servers> <#services> <rigid res. sigma> <fluid res. sigma> <slack> <fluid_mean> <need_cv> <prob of sla> <sla value>"
    sys.exit()

seed        = int(sys.argv[1])
numrigid    = int(sys.argv[2])
numfluid    = int(sys.argv[3])
numservers  = int(sys.argv[4])
numservices = int(sys.argv[5])
rigidressigma = float(sys.argv[6])
fluidressigma = float(sys.argv[7])
slack       = float(sys.argv[8])
fluidmean   = float(sys.argv[9])
needcv      = float(sys.argv[10])
probsla     = float(sys.argv[11])
slavalue    = float(sys.argv[12])

random.seed(seed)

if rigidressigma < 0.0:
    print "rigid resource variation must be larger than 0.0"
    sys.exit()

if fluidressigma < 0.0:
    print "fluid resource variation must be larger than 0.0"
    sys.exit()

if slack < 0.0 or slack > 1.0:
    print "slack must be between 0 and 1."
    sys.exit()

if fluidmean < 0.0 or fluidmean > 1.0:
    print "cpu mean must be between 0 and 1."
    sys.exit()

if probsla < 0.0 or probsla > 1.0:
    print "probsla must be between 0 and 1."
    sys.exit()

if slavalue < 0.0 or slavalue > 1.0:
    print "slavalue must be between 0 and 1."
    sys.exit()

if numservers >= numservices:
    print "there should be more services than servers"
    sys.exit()

rigidres = []
fluidres = []

for j in range(numrigid):
    rigidres.append([])
    for i in range(numservers):
        rigidres[j].append(-1.0)
        while rigidres[j][i] < 0.0:
            rigidres[j][i] = random.normalvariate(1.0, rigidressigma)

for j in range(numrigid):
    resmax = max(rigidres[j])
    for i in range(numservers):
        rigidres[j][i] = max(0.01, rigidres[j][i] / resmax)

for j in range(numfluid):
    fluidres.append([])
    for i in range(numservers):
        fluidres[j].append(-1.0)
        while fluidres[j][i] < 0.0:
            fluidres[j][i] = random.normalvariate(1.0, fluidressigma)

for j in range(numfluid):
    resmax = max(fluidres[j])
    for i in range(numservers):
        fluidres[j][i] = max(0.01, fluidres[j][i] / resmax)

rigidmean   = 0.5
rigidsigma  = rigidmean * needcv
fluidsigma  = fluidmean * needcv

rigidreq = []
fluidreq = []
sla = []
rigidneedtotal = []
fluidneedtotal = []

for i in range(numservices):

    rigidreq.append([])
    fluidreq.append([])

    sla.append(0.0)
    if random.random() <= probsla:
        sla[i] = slavalue

    # rigid needs
    for j in range(numrigid):
        rigidneedtotal.append(0.0)
        rigidreq[i].append(0.0)
        while rigidreq[i][j] < 0.01 or rigidreq[i][j] > 1.0:
            rigidreq[i][j] = random.normalvariate(rigidmean, rigidsigma)
            rigidneedtotal[j] += rigidreq[i][j]

    # fluid needs
    for j in range(numfluid):
        fluidneedtotal.append(0.0)
        fluidreq[i].append(0.0)
        while fluidreq[i][j] < 0.01 or fluidreq[i][j] > 1.0:
            fluidreq[i][j] = random.normalvariate(fluidmean, fluidsigma)
            fluidneedtotal[j] += sla[i] * fluidreq[i][j]

rigidscale = []
for j in range(numrigid):
	rigidscale.append(1.0)
	rigidscale[j] = (1.0 - slack) * sum(rigidres[j]) / rigidneedtotal[j];

# Computation of a fluid scale, but really doesn't make sense
fluidscale = []
for j in range(numrigid):
    fluidscale.append(1.0)
    if fluidneedtotal[j] > 0:
        fluidscale[j] = (1.0 - slack) * sum(fluidres[j]) / fluidneedtotal[j];

print numrigid
print numfluid
print numservers
print numservices
for i in range(numservers):
    for j in range(numrigid):
        print "%.3f" % (rigidres[j][i]),
    for j in range(numfluid):
        print "%.3f" % (fluidres[j][i]),
    print "\n",

for i in range(numservices):
    # print the SLA
    print "%.3f" % (sla[i]),
    # print the rigid needs
    for j in range(numrigid):
        print "%.3f" % (min(1.0,rigidreq[i][j]*rigidscale[j])),
	# print the fluid needs
    for j in range(numfluid):
        print "%.3f" % (min(1.0,fluidreq[i][j]*1.0)),
#       	print "%.3f" % (min(1.0,fluidreq[i][j]*fluidscale[j])),
    print "\n",

