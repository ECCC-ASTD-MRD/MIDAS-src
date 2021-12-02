#! /usr/bin/env python3
'''
Analyse the dependency tree of modules and programs to determine which programs
are impacted by a new external dependency in a given module.


SYNOPSIS
from /src directory of a MIDAS project

```
./findDependentAbs.py ${moduleNameWithout_mod_dot_o} ${buildDir}
```

ARGUMENTS
    ${moduleNameWithout_mod_dot_o}  the module *name* (without '_mod.o') in
                                    which new external dependencies have been 
                                    added.
    ${buildDir}                     the current build directory on the frontend

EXAMPLE
```
./findDependentAbs.py varqc ../compiledir/midas_bld-$(../midas.version.sh)/ubuntu-18.04-skylake-64/intel-19.0.3.199/
```
'''

import re
import sys

patternDepTarget = r'^(.*) :.*$'

def recurseModDependsOnMod(module, depFile='dep.obj.inc', uniqDep=None, verbose=False):
    '''
    Find all dependent target modules of a given module, recurse with found
    dependent targets.
    '''
    if verbose: print(f'-- searching for {module}')
    if uniqDep == None: uniqDep=set()
    patternLook4Mod = f'^.*:.*{module}_mod.o.*$' 
    patternMod = r'_mod.o'
    with open(depFile) as fun:
        for line in fun.readlines():
            #-- check for module in dependencies
            if re.match(patternLook4Mod, line):
                #-- extract dependant target
                depObject = re.findall(patternDepTarget, line)[0]
                if verbose : print(depObject)
                #-- is it a module?
                if re.match(r'.*'+patternMod, depObject): 
                    depObject = re.sub(patternMod,'', depObject)
                    if depObject in uniqDep : continue
                    uniqDep.add(depObject)
                    uniqDep.union(recurseModDependsOnMod(  
                                        depObject, depFile=depFile,
                                        uniqDep=uniqDep, verbose=verbose))
    if verbose: print(f'-- end {module}')
        
    return uniqDep

def findAbsDependsOnMod(module, depFile='dep.abs.inc', verbose=False):
    if verbose: print(f'-- searching for absolute depending on {module}')
    patternLook4Mod = f'^.*:.*{module}_mod.o.*$' 
    absList = list()
    with open(depFile) as fun:
        for line in fun.readlines():
            #-- check for module in dependencies
            if re.match(patternLook4Mod, line):
                #-- extract dependant target
                depAbs = re.findall(patternDepTarget, line)[0]
                depAbs = re.sub(r'.o', '', depAbs)
                absList.append(depAbs)
    return absList

###| command arguments reading |#####################################

if __name__ == '__main__':
    from glob import glob
    module = sys.argv[1]
    buildDir = sys.argv[2]
    outputModules = False
    
    try:
        depFileObj = glob(buildDir+'/dep.obj.inc')[0]
        depFileAbs = glob(buildDir+'/dep.abs.inc')[0]
    except: 
        raise RuntimeError('Inexistant directory', buildDir)

    dependentModules = recurseModDependsOnMod(  module, depFile=depFileObj, 
                                                verbose=False)

    if dependentModules : 
        if outputModules: 
            print(f'The following modules depends on {module}:')
            for mod in dependentModules:
                print(f'  * {mod}')    

        setAbs = list()
        for mod in dependentModules:
            setAbs.extend(findAbsDependsOnMod(mod, depFile=depFileAbs))
        setAbs=set(setAbs)

        if setAbs: 
            print(f'The following absolutes depends on {module}:')
            for absolute in setAbs:
                print(f'  * {absolute}')
        else:
            print(f'No absolutes depends on {module}')
    else:
        print(f'No modules depends on {module}')

# vim: set ts=4 sw=4:
