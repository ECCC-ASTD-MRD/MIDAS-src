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
useslist[$numModules]=""
numberuses[$numModules]=0
categories[$numModules]=0
prefixes[$numModules]=""

for fullfilename in modules/*_mod.f* modules/*/*_mod.f*; do
  filename=`basename $fullfilename`
  modulenamelist=`grep -i "^ *module *[a-zA-Z]" $fullfilename |grep -iv "^ *module *procedure *" | sed 's/module //Ig' | tr '[:upper:]' '[:lower:]'`
  for modulename in $modulenamelist; do
    numModules=$((numModules + 1))
    [ "${verbose}" = "yes" ] && echo $numModules $fullfilename $filename $modulename
    uses=`grep -A 10000000 -iE "^ *module *${modulename}" $fullfilename |grep -B 10000000 -iE "^ *end *module *${modulename}" |grep -i "^ *use *.*_mod\>" | sed 's/, *only *:.*//Ig' | sed 's/!.*//Ig' | sed 's/use //Ig' | tr '[:upper:]' '[:lower:]' | sort -u`
    [ "${verbose}" = "yes" ] && echo $uses
    # This assumes only 1 module per file (usually the case, except for obsspacedata_mod.f90)
    category=`grep "category=" $fullfilename | sed "s/.*category=['\"]\([0-9]*\).*/\1/"`
    prefix=`grep "prefix=" $fullfilename | sed "s/.*prefix=['\"]\([a-z0-9]*\).*/\1/"`
    # If multiple modules found in the file, then assume category and prefix info is inside each module
    if [ `echo $category |wc -w` -gt 1 ]; then
      category=`grep -A 10000000 -iE "^ *module *${modulename}" $fullfilename |grep -B 10000000 -iE "^ *end *module *${modulename}" |grep "category=" | sed "s/.*category=['\"]\([0-9]*\).*/\1/"`
      prefix=`grep -A 10000000 -iE "^ *module *${modulename}" $fullfilename |grep -B 10000000 -iE "^ *end *module *${modulename}" |grep "prefix=" | sed "s/.*prefix=['\"]\([a-z0-9]*\).*/\1/"`
    fi
    [ "${verbose}" = "yes" ] && echo

    filenames[$numModules]=$filename
    fullfilenames[$numModules]=$fullfilename
    modulenames[$numModules]=$modulename
    useslist[$numModules]=$uses
    numberuses[$numModules]=`echo "${useslist[$numModules]}" |wc -w`
    categories[$numModules]=$category
    prefixes[$numModules]=$prefix
  done
  
done

echo "Finished scanning $numModules modules"
echo

make_index modulename_index "${modulenames[@]}"
