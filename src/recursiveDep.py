#! /usr/bin/env python3
'''
Analyse dependencies of individual absolute files to deduce implicit dependencies
contained in them that will need to be added at link time

Usage: 
   recursiveDep.py TARGETLIST 
   [-v] [-s=<str>] [-d=<file>]
   recursiveDep.py -h | --help

Arguments:
    TARGETLIST  target absolute file list

Options:
    -h, --help  show this screen and exit.
    -v          toggle verbose mode
    -s=<str>    absolute suffix [default: .Abs]
    -d=<file>   prerequisite makefile include file [default: ./depend.inc]
'''


import re
import sys
from docopt import docopt, DocoptExit


def recurseDep(target, depFile='depend.inc', uniqDep=None, verbose=False):
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
try:
    opt = docopt(__doc__, sys.argv[1:])
except DocoptExit as inst:
    print('Invalid command!')
    print(inst)
    sys.exit(1)


targetList = (opt['TARGETLIST']).split()
depFile = opt['-d']
suffix = opt['-s']
verbose = opt['-v']

if verbose: 
    print('Evaluating linking prerequisite for:')
    print(targetList)

depDict = dict()
for target in targetList:
    obj=target.replace(suffix, '.o')
    depDict[target] = recurseDep(   obj, depFile=depFile, verbose=verbose)

    ## output to stdout
    print('%s : %s'%(target, obj), end=' ')
    for f in depDict[target]:
        print(f, end=' ')
    print()

# vim: set ts=4 sw=4:
