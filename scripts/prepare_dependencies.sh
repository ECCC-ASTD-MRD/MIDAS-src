#!/bin/bash

verbose="no"

echo
echo "LOOP OVER ALL MODULES AND SAVE LIST OF USED MODULES"

make_index () {
  local index_name=$1
  shift
  local -a value_array=("$@")
  local i
  # -A means associative array, -g means create a global variable:
  declare -g -A ${index_name}
  for i in "${!value_array[@]}"; do
    eval ${index_name}["${value_array[$i]}"]=$i
  done
}

numModules=0

# initialize zero'th element as blank (used for modules external to MIDAS)
filenames[$numModules]=""
fullfilenames[$numModules]=""
modulenames[$numModules]="UNKNOWN"
modulenames_lc[$numModules]="UNKNOWN"
useslist[$numModules]=""
revmodlist[$numModules]=""
revpgmlist[$numModules]=""
numberuses[$numModules]=0
revnumberuses[$numModules]=0
categories[$numModules]=0
prefixes[$numModules]=""

for fullfilename in modules/*_mod.f*90; do
  filename=`basename $fullfilename`
  modulenamelist=`grep -i "^ *module *[a-zA-Z]" $fullfilename \
    | grep -iv "^ *module *procedure *" | sed 's/module //Ig'`
  for modulename in $modulenamelist; do
    modulename_lc=`echo ${modulename} | tr '[:upper:]' '[:lower:]'`
    usedbymod=''
    usedbypgm=''
    numModules=$((numModules + 1))
    #-- gather all forward dependencies of the module
    uses=`grep -A 10000000 -iE "^ *module *${modulename_lc}" $fullfilename \
      | grep -B 10000000 -iE "^ *end *module *${modulename_lc}" \
      | grep -i "^ *use *.*_mod\>" | sed 's/, *only *:.*//Ig' | sed 's/!.*//Ig' \
      | sed 's/use //Ig' | tr '[:upper:]' '[:lower:]' | sort -u`

    #-- gather all reverse dependencies (all modules and programs using the module)
    usedbymod_files=`grep -il "^ *use *${modulename_lc}" modules/*_mod.f*90`
    if [ ! -z "${usedbymod_files}" ]; then 
      for modfile in ${usedbymod_files}; do
        lines=`grep -in "^ *use *${modulename_lc}" ${modfile} | cut -d':' -f1`
        for line in ${lines}; do
          # find which module use it
          usedbymod="${usedbymod} `cat ${modfile} | head -n ${line} \
            | tr '[:upper:]' '[:lower:]' | grep -v 'module *procedure' \
            | grep '^ *module ' | tail -n 1 | sed -n 's/^ *module *\([^!]*\).*/\1/p'`"
        done
      done
    fi
    usedbypgm_files=`grep -il "^ *use *${modulename_lc}" programs/*.f90`
    if [ ! -z "${usedbypgm_files}" ]; then
      for pgmfile in ${usedbypgm_files}; do
        usedbypgm="${usedbypgm} `grep -i '^ *program .*$' ${pgmfile} \
          | tr '[:upper:]' '[:lower:]' | grep -v 'end program' \
          | sed -n 's/^ *program *\([^!]*\).*/\1/p'`"
      done
    fi
    # This assumes only 1 module per file (usually the case, except for obsspacedata_mod.f90)
    category=`grep "category=" $fullfilename | sed "s/.*category=['\"]\([0-9]*\).*/\1/"`
    prefix=`grep "prefix=" $fullfilename | sed "s/.*prefix=['\"]\([a-z0-9]*\).*/\1/"`
    # If multiple modules found in the file, then assume category and prefix info is inside each module
    if [ `echo $category |wc -w` -gt 1 ]; then
      category=`grep -A 10000000 -iE "^ *module *${modulename_lc}" $fullfilename \
        | grep -B 10000000 -iE "^ *end *module *${modulename_lc}" |grep "category=" \
        | sed "s/.*category=['\"]\([0-9]*\).*/\1/"`
      prefix=`grep -A 10000000 -iE "^ *module *${modulename_lc}" $fullfilename \
        | grep -B 10000000 -iE "^ *end *module *${modulename_lc}" |grep "prefix=" \
        | sed "s/.*prefix=['\"]\([a-z0-9]*\).*/\1/"`
    fi

    filenames[$numModules]=$filename
    fullfilenames[$numModules]=$fullfilename
    modulenames[$numModules]=$modulename
    modulenames_lc[$numModules]=$modulename_lc
    useslist[$numModules]=$uses
    revmodlist[$numModules]=${usedbymod}
    revpgmlist[$numModules]=${usedbypgm}
    numberuses[$numModules]=`echo "${useslist[$numModules]}" |wc -w`
    revnumberuses[$numModules]=$(( $(echo ${revmodlist[$numModules]} | wc -w)+ $(echo ${revpgmlist[$numModules]} | wc -w)))
    categories[$numModules]=$category
    prefixes[$numModules]=$prefix

    [ "${verbose}" = "yes" ] && echo "$numModules $fullfilename $filename $modulename_lc ${revnumberuses[$numModules]} ($(echo ${revmodlist[$numModules]} | wc -w)+$(echo ${revpgmlist[$numModules]} | wc -w))"
    [ "${verbose}" = "yes" ] && echo
  done
  
done

echo "Finished scanning $numModules modules"
echo

make_index modulename_index "${modulenames_lc[@]}"
