#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import subprocess, os

minimumTimingLimit = 3.0

if len(sys.argv) < 2 or sys.argv[1] == "--help" or sys.argv[1] == "-h":
    print()
    print(" #***************************************************************************#")
    print(" #                                                                           #")
    print(" # timingTool.py                                                             #")
    print(" #                                                                           #")
    print(" # example of usage for comparing 2 listing files:                           #")
    print(" #                                                                           #")
    print(" # timingTool.py file1 file2 1 > timings1.dat  # extracts timings from file1 #")
    print(" # timingTool.py file2 file1 2 > timings2.dat  # extracts timings from file2 #")
    print(" # xxdiff timings1.dat timings2.dat                                          #")
    print(" #                                                                           #")
    print(" #***************************************************************************#")
    print()
    exit()


filename = sys.argv[1]

# this allows the user to specify a reference listing - timing information
# will not be generated from this file, but it will be used to ensure a 
# common list (and order) of the labels to facilitate an xxdiff of files
# from running the script twice, once for each listing file
if len(sys.argv) > 2:
    filename2 = sys.argv[2]
else:
    filename2 = ''

# this allows the user to specify in which order to construct the 
# list of labels when a reference listing is provided
if len(sys.argv) > 3:
    exptOrder = sys.argv[3]
else:
    exptOrder = "1"

print(f"Processing the file: {filename}")
if filename2:
    print(f"Reference file supplied: {filename2}")
    if exptOrder == "2":
        print("Reference file used first for constructing label list")
    else:
        print("Main file used first for constructing label list")

if not os.path.exists(filename):
    print("The file does not exist!")
    exit()

allLines = []
with open(filename, "r") as ins:
    for line in ins:
        strIndex = line.find(" TMG:")
        if line.find(" TMG:") == -1 or line.find("oe-") == -1: 
            continue
        allLines.append(line.rstrip('\n'))

allLines2 = []
if os.path.isfile(filename2):
    with open(filename2, "r") as ins2:
        for line in ins2:
            strIndex = line.find(" TMG:")
            if line.find(" TMG:") == -1 or line.find("oe-") == -1: 
                continue
            allLines2.append(line.rstrip('\n'))

# get unique and complete list of labels
labelListOrig = []
if exptOrder == "1":
    # get all labels from main file, first
    for line in allLines:
        label = line.partition("LABEL=")[2].partition(",")[0].strip()
        if label not in labelListOrig:
            labelListOrig.append(label)

    # add new labels from reference file, if supplied
    if os.path.isfile(filename2):
        for line in allLines2:
            label = line.partition("LABEL=")[2].partition(",")[0].strip()
            if label not in labelListOrig:
                labelListOrig.append(label)
else:
    # start with reference file
    if os.path.isfile(filename2):
        for line in allLines2:
            label = line.partition("LABEL=")[2].partition(",")[0].strip()
            if label not in labelListOrig:
                labelListOrig.append(label)

    # now add extra labels from main file
    for line in allLines:
        label = line.partition("LABEL=")[2].partition(",")[0].strip()
        if label not in labelListOrig:
            labelListOrig.append(label)

# remove labels from the list if the max timing is less that a minimum threshold
labelList = []
for label in labelListOrig:
    timingList = []
    for line in allLines:
        if label == line.partition("LABEL=")[2].partition(",")[0].strip():
            timing = float(line.partition("TIME=")[2].partition("seconds")[0].strip())
            timingList.append(timing)
    if allLines2:
        for line in allLines2:
            if label == line.partition("LABEL=")[2].partition(",")[0].strip():
                timing = float(line.partition("TIME=")[2].partition("seconds")[0].strip())
                timingList.append(timing)
    # check maximum timing
    if max(timingList) > minimumTimingLimit:
        labelList.append(label)

# get min/max/mean timing and counts for each label in labelList
minTimingList = []
maxTimingList = []
meanTimingList = []
minCountList = []
maxCountList = []
meanCountList = []
timingLists = []
countLists = []
for label in labelList:
    timingList = []
    countList = []
    for line in allLines:
        if label == line.partition("LABEL=")[2].partition(",")[0].strip():
            timing = float(line.partition("TIME=")[2].partition("seconds")[0].strip())
            count  = int(line.partition("COUNT=")[2].partition(",")[0].strip())
            timingList.append(timing)
            countList.append(count)
    # timing stats
    if timingList: 
        minTimingList.append(min(timingList))
        maxTimingList.append(max(timingList))
        meanTimingList.append(sum(timingList)/len(timingList))
    else:
        minTimingList.append(0.0)
        maxTimingList.append(0.0)
        meanTimingList.append(0.0)

    # count stats
    if countList:
        minCountList.append(min(countList))
        maxCountList.append(max(countList))
        meanCountList.append(sum(countList)/len(countList))
    else:
        minCountList.append(0)
        maxCountList.append(0)
        meanCountList.append(0)

    # complete info
    timingLists.append((timingList))
    countLists.append((countList))

print()
print('{:30}'.format('LABEL'), '{:>8}'.format('minTime'), '{:>8}'.format('maxTime'), '{:>8}'.format('meanTime'), '{:>10}'.format('minCount'), '{:>10}'.format('maxCount'), '{:>10}'.format('meanCount'))
print('{:30}'.format('====='), '{:>8}'.format('======='), '{:>8}'.format('======='), '{:>8}'.format('========'), '{:>10}'.format('========'), '{:>10}'.format('========'), '{:>10}'.format('========='))
for listIndex in range(0,len(labelList)):
    print('{:30}'.format(labelList[listIndex]), \
          '{:8.2f}'.format(minTimingList[listIndex]), \
          '{:8.2f}'.format(maxTimingList[listIndex]), \
          '{:8.2f}'.format(meanTimingList[listIndex]), \
          '{:10d}'.format(minCountList[listIndex]),  \
          '{:10d}'.format(maxCountList[listIndex]),  \
          '{:10d}'.format(round(meanCountList[listIndex])))
