#!/usr/bin/env python3

import sys
import re

if len(sys.argv)==1:
    inputFile = sys.stdin
else:
    inputFile = open(sys.argv[1])

regex = re.compile('^oe-([0-9]{5}-[0-9]{5}): *Memory Used: *([0-9]+) *Mb *$')
mpiTiles = {}
for line in sys.stdin.readlines():
    match = regex.match(line)
    if match:
        mpi = match.group(1)
        memoryStr = match.group(2)
        memory= int(memoryStr)
        if mpi in mpiTiles:
            if mpiTiles[mpi]<memory:
                mpiTiles[mpi] = memory
        else:
            mpiTiles[mpi] = memory

print("Here is the maximum memory per tile")
mpiKeys = list(mpiTiles.keys())
mpiKeys.sort()
for mpi in mpiKeys:
    print("\t{}\t{}".format(mpi,mpiTiles[mpi]))
