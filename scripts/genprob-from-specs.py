#!/usr/bin/python

import random
import sys

if len(sys.argv) < 6:
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
serverspecs = []
max_server_select = int(serverf.readline())
for line in serverf.readlines():
    fields = line.split();
    serverspecs.append((int(fields[0]), fields[1:]))
serverf.close()

for i in range(numservers):
    server_select = random.randint(1, max_server_select)
    j = 0
    while (j < len(serverspecs) and serverspecs[j][0] < server_select):
        j += 1
    for val in serverspecs[j-1][1]:
        print val,
    print

servicef = open(service_fname, 'r')
servicespecs = []
max_service_select = int(servicef.readline())
for line in servicef.readlines():
    fields = line.split()
    servicespecs.append((int(fields[0]), fields[1:]))
servicef.close()

for i in range(numservices):
    service_select = random.randint(1, max_service_select)
    j = 0
    while (j < len(servicespecs) and servicespecs[j][0] < service_select):
        j += 1
    for val in servicespecs[j-1][1]:
        print val,
    print
