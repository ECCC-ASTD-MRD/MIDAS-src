#!/bin/bash

set -e

## Build the module list use by program 'obsIO'
program=${1:-obsIO.f90}
toplevel=$(git rev-parse --show-toplevel)
#uses=`grep -i '^ *use *.*_mod' ${toplevel}/src/programs/${program} | sed 's/, *only *:.*//Ig' | sed 's/!.*//Ig' | sed 's/use //Ig' | tr '[:upper:]' '[:lower:]' | sort -u`
function getMod {
    set -e
    __getMod_path__=${1}
    grep -i '^ *use *.*_mod' ${__getMod_path__} | sed 's/, *only *:.*//Ig' | sed 's/!.*//Ig' | sed 's/use //Ig' | tr '[:upper:]' '[:lower:]' | sort -u
    unset __getMod_path__
}
uses=$(getMod ${toplevel}/src/programs/${program})

function findModuleDefinition {
    set -e
    __findModuleDefinition_mod__=$1
    git grep -il --full-name "^ *module *${__findModuleDefinition_mod__}" ${toplevel}/src/modules || true
    unset __findModuleDefinition_mod__
}
iteration=1
function findAllModulesUsed {
    set -e
    __findAllModulesUsed_uses__=$*
    __findAllModulesUsed_thisuse__=
    for __findAllModulesUsed_mod__ in ${__findAllModulesUsed_uses__}; do
        __findAllModulesUsed_modpath__=$(findModuleDefinition ${__findAllModulesUsed_mod__})
        [ -z "${__findAllModulesUsed_modpath__}" ] && continue
        __findAllModulesUsed_moduse__=$(getMod ${toplevel}/${__findAllModulesUsed_modpath__})
        #echo "__findAllModulesUsed_moduse__='${__findAllModulesUsed_moduse__}'" >&2
        #echo ${__findAllModulesUsed_moduse__} | grep -sqi rmatrix && printf "appel a rmatrix: '${__findAllModulesUsed_mod__}' '${__findAllModulesUsed_moduse__}'\n\n" >&2
        if [ -n "${__findAllModulesUsed_thisuse__}" ]; then
            __findAllModulesUsed_thisuse__=$(echo ${__findAllModulesUsed_thisuse__} ${__findAllModulesUsed_mod__} ${__findAllModulesUsed_moduse__} | tr ' ' '\n' | sort -u | tr '\n' ' ' | sed 's/ *$//')
        else
            __findAllModulesUsed_thisuse__="${__findAllModulesUsed_moduse__}"
        fi
        unset __findAllModulesUsed_modpath__
    done
    if [ "${__findAllModulesUsed_uses__}" = "${__findAllModulesUsed_thisuse__}" ]; then
        echo ${__findAllModulesUsed_uses__}
    else
        #echo ${__findAllModulesUsed_uses__} >&2
        #echo ${__findAllModulesUsed_thisuse__} >&2
        [ "${iteration}" -ge 10 ] && exit 1
        let iteration=iteration+1
        findAllModulesUsed ${__findAllModulesUsed_thisuse__}
    fi
    unset __findAllModulesUsed_uses__ __findAllModulesUsed_thisuse__ __findAllModulesUsed_moduse__ __findAllModulesUsed_mod__
}
modList=$(findAllModulesUsed ${uses})
echo ${modList}
#echo ${modList} | tr ' ' '\n'

## echo
## echo "We should have this list for program 'obsIO'"
## echo burpfiles_mod.mod
## echo burpread_mod.mod
## echo clib_interfaces_mod.mod
## echo cmafiles_mod.mod
## echo codeprecision_mod.mod
## echo codtyp_mod.mod
## echo earthconstants_mod.mod
## echo indexlistdepot_mod.mod
## echo mathphysconstants_mod.mod
## echo mpi_mod.mod
## echo mpivar_mod.mod
## echo obscolumnnames_mod.mod
## echo obsdatacolumn_mod.mod
## echo obsfiles_mod.mod
## echo obsspacedata_mod.mod
## echo obssubspacedata_mod.mod
## echo obsutil_mod.mod
## echo ramdisk_mod.mod
## echo sqlitefiles_mod.mod
## echo sqliteread_mod.mod
## echo timecoord_mod.mod
## echo utilities_mod.mod
