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
    if verbose: print(f'<<< {target}')
    if uniqDep == None: uniqDep=set()
    pattern = f'^{target} :'
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
                        if verbose: print(f'--search {f}')
                        uniqDep.union(
                            recurseDep( f, depFile=depFile,
                                        uniqDep=uniqDep, verbose=verbose))
                    if verbose: print(f'--end {f}')
                return uniqDep
        raise Exception(f'{target} not in {depFile}')

###| command arguments reading |#####################################

if __name__ == '__main__':
    verbose=False
    depFile = sys.argv[1]
    targetList = sys.argv[2]
    if len(sys.argv) > 3:
        verbose=bool(sys.argv[3])
    targetList = targetList.split()
    suffix = '.Abs'
    
    depDict = dict()
    for target in targetList:
        obj=target.replace(suffix, '.o')
        depDict[target] = recurseDep(obj, depFile=depFile, verbose=verbose)
    
        ## output to stdout
        print(f'{target} : {obj}', end=' ')
        for f in depDict[target]:
            print(f, end=' ')
        print()

# vim: set ts=4 sw=4:
