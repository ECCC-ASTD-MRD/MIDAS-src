#! /usr/bin/env python3
'''
Analyse dependencies of individual absolute files to deduce implicit dependencies
contained in them that will need to be added at link time

Usage: 
   recursiveDep.py OBJ_DEP TARGETLIST 

Arguments:
    TARGETLIST  target absolute file list
    OBJ_DEP     prerequisite makefile include file
'''


import re
import sys


def recurseDep(target, depFile='dep.obj.inc', uniqDep=None, verbose=False):
    if verbose: print('<<< %s'%target)
    if uniqDep == None: uniqDep=set()
    pattern = '^[^ ]*%s :'%target
    with open(depFile) as fun:
        for line in fun.readlines():
            if re.match(pattern, line):
                dep = re.sub(pattern, '', line)
                elems = dep.split()
                for f in elems:
                    if f in uniqDep: 
                        continue
                    if verbose: print(f)
                    if re.match('^[^ ]*\.o', f):
                        uniqDep.add(f)
                        uniqDep.union(
                            recurseDep( f, depFile=depFile,
                                        uniqDep=uniqDep, verbose=verbose))
                return uniqDep
        raise Exception('%s not in %s'%(target, depFile))

###| command arguments reading |#####################################

depFile = sys.argv[1]
targetList = sys.argv[2]
targetList = targetList.split()
suffix = '.Abs'

depDict = dict()
for target in targetList:
    obj=target.replace(suffix, '.o')
    depDict[target] = recurseDep(   obj, depFile=depFile)

    ## output to stdout
    print('%s : %s'%(target, obj), end=' ')
    for f in depDict[target]:
        print(f, end=' ')
    print()

# vim: set ts=4 sw=4:
