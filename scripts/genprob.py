#!/usr/bin/python

import random
import sys

if len(sys.argv) < 11:
        print "usage:"
        print sys.argv[0] + \
                " <seed> <#rigid> <#fluid> <#servers> <#services> <slack> <fluid_mean> <need_cv> <prob of sla> <sla value>"
        sys.exit()

seed        = int(sys.argv[1])
numrigid    = int(sys.argv[2])
numfluid    = int(sys.argv[3])
numservers  = int(sys.argv[4])
numservices = int(sys.argv[5])
slack  = float(sys.argv[6])
fluidmean   = float(sys.argv[7])
needcv     = float(sys.argv[8])
probsla     = float(sys.argv[9])
slavalue    = float(sys.argv[10])

random.seed(seed)


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


rigidmean = 0.5

rigidsigma  = rigidmean * needcv
fluidsigma  = fluidmean * needcv

rigidreq = []
fluidreq = []
sla = []
rigidtotal = []
fluidtotal = []


for i in range(numservices):

	rigidreq.append([])
	fluidreq.append([])

	# Is there an SLA?
        sla.append(0.0)
	if random.random() <= probsla:
		sla[i] = slavalue

	# rigid needs
	for j in range(numrigid):
          rigidtotal.append(0.0)
          rigidreq[i].append(0.0)
          while rigidreq[i][j] < 0.01 or rigidreq[i][j] > 1.0:
                  rigidreq[i][j] = random.normalvariate(rigidmean, rigidsigma)
          rigidtotal[j] += rigidreq[i][j]

	# fluid needs
	for j in range(numfluid):
          fluidtotal.append(0.0)
          fluidreq[i].append(0.0)
          while fluidreq[i][j] < 0.01 or fluidreq[i][j] > 1.0:
                  fluidreq[i][j] = random.normalvariate(fluidmean, fluidsigma)
          fluidtotal[j] += sla[i] * fluidreq[i][j]

 
rigidscale = []
for j in range(numrigid):
	rigidscale.append(1.0)
	rigidscale[j] = (1.0 - slack) * (numservers ) / rigidtotal[j];

# Computation of a fluid scale, but really doesn't make sense
fluidscale = []
for j in range(numrigid):
	fluidscale.append(1.0)
        if fluidtotal[j] > 0:
		fluidscale[j] = (1.0 - slack) * (numservers ) / fluidtotal[j];


print numrigid
print numfluid
print numservers
print numservices
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


