#! /usr/bin/env python3
'''
Analyze the dependency tree of modules and programs.

Usage:
    findDependentAbs.py OBJECT -a [ -v ] [ --path=<str> ]
    findDependentAbs.py OBJECT -c [ -v ] [ --path=<str> ]
    findDependentAbs.py -h | --help

Arguments:
    OBJECT  object name (module or program and without `_mod` or suffix)

Options:
    -h, --help      show this help and exit
    -v              toggle verbose mode
    -a              list which programs are impacted by a new external
                        dependency in the module OBJECT
    -c              show compilation order to build OBJECT
    --path=<str>    explicit path to build directory
'''
from docopt import docopt, DocoptExit
import re
import sys
import os
from glob import glob
import subprocess

patternDepTarget = r'^(.*) :.*$'

def recurseModDependsOnMod( module, depFile='dep.obj.inc', uniqDep=None,
                            verbose=False):
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
            if re.match(patternLook4Mod, line, re.IGNORECASE):
                #-- extract dependant target
                depObject = re.findall( patternDepTarget, line, 
                                        re.IGNORECASE)[0]
                if verbose : print(depObject)
                #-- is it a module?
                if re.match(r'.*'+patternMod, depObject, re.IGNORECASE): 
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
            if re.match(patternLook4Mod, line, re.IGNORECASE):
                #-- extract dependant target
                depAbs = re.findall(patternDepTarget, line, re.IGNORECASE)[0]
                depAbs = re.sub(r'.o', '', depAbs)
                absList.append(depAbs)
    return absList


def recurseCompilationOrder(module, depFile='dep.obj.inc', order=None,
                            lines=None, verbose=False):
    '''
    Establish  compilation order by recursing into prerequisites
    '''
    module = module.lower()
    if order == None : order = list()
    if lines == None :
        fun = open(depFile)
        lines=fun.readlines()
        fun.close()
    patternLook4Mod = f'^{module}.* :.*$'
    patternMod = r'_mod.o'
    if verbose :
        print(f'-- listing {module} prerequisites')
    for line in lines:
        if re.match(patternLook4Mod, line, re.IGNORECASE):
            if verbose: print(line)
            prereqList=line.split(':')[1].split()[1:]
            moduleList = list()
            for prereq in prereqList:
                moduleList.append(re.sub(patternMod, '', prereq))
            for prereq in moduleList:
                if prereq not in order:
                    order=recurseCompilationOrder(  prereq, depFile=depFile,
                                                    order=order, lines=lines,
                                                    verbose=verbose)
                    order.append(prereq)
            break
    return order


###| command arguments reading |#####################################
try:                                                                            
    opt = docopt(__doc__, sys.argv[1:])                                         
except DocoptExit as inst:                                                      
    print('Invalid command!')                                                   
    print(inst)                                                                 
    sys.exit(1)

module = opt['OBJECT'].lower()

buildDir = opt['--path']
if not buildDir:
    encoding=sys.getfilesystemencoding()
    cfgFile='config.dot.sh'
    envLeafBldDir='EC_ARCH'
    print(f'... sourcing {cfgFile}, use explicit --path to make it faster')
    version=subprocess.check_output('../midas.version.sh').strip().decode(encoding)
    relCompDir='../compiledir/'
    if not os.path.isdir(relCompDir):
        print(f'{relCompDir} does not exist.  Dependency files not accessible')
    ## -- loading compiler from SSM to get initialized EC_ARCH
    proc=subprocess.Popen(
        ['bash', '-c', f'source {cfgFile} && echo ${envLeafBldDir}'], 
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    proc.wait()
    ec_arch=proc.stdout.readlines()[-1].strip().decode(encoding)
    buildDir = relCompDir+'midas_bld-'+version+'/'+ec_arch
    print(f'Next time call with\n\t--path {buildDir}')

    if not os.path.isdir(buildDir):
        print(f'{buildDir} does not exist')
        print(f'Dependencies probably not built. Attempting to build them')
        proc=subprocess.Popen( ['bash', '-c', f'source {cfgFile} && make depend']) 
        proc.wait()
try:
    depFileObj = glob(buildDir+'/dep.obj.inc')[0]
    depFileAbs = glob(buildDir+'/dep.abs.inc')[0]
except: 
    raise RuntimeError('Inexistant directory', buildDir)

verbose=opt['-v']


if opt['-a']:
    dependentModules = list(recurseModDependsOnMod( module,
                                                    depFile=depFileObj,
                                                    verbose=verbose))
    if dependentModules : 
        dependentModules.sort()
        if verbose: 
            print(f'The following modules depends on {module}:')
            for mod in dependentModules:
                print(f'  * {mod}')    
    
        setAbs = list()
        for mod in dependentModules:
            setAbs.extend(findAbsDependsOnMod(mod, depFile=depFileAbs))
        setAbs=list(set(setAbs))
        setAbs.sort()
    
        if setAbs: 
            print(f'The following absolutes depends on {module}:')
            for absolute in setAbs:
                print(f'  * {absolute}')
        else:
            print(f'No absolutes depends on {module}')
    else:
        print(f'No modules depends on {module}')

elif opt['-c']:
    
    order = recurseCompilationOrder(module, depFile=depFileObj, verbose=verbose) 
    if order:
        print(f'Building {module} will result in compiling these objects:')
        for i, prereq in enumerate(order):
            print(f'{i+1} : {prereq}')
    else:
        print(f'{module} has no prerequisites') 
        
# vim: set ts=4 sw=4:
