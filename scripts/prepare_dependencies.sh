#!/bin/bash

echo
echo "LOOP OVER ALL MODULES AND SAVE LIST OF USED MODULES"
echo

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

for fullfilename in modules/*_mod.f* modules/*/*_mod.f*; do
  filename=`basename $fullfilename`
  modulenamelist=`grep -i "^ *module *[a-zA-Z]" $fullfilename |grep -iv "^ *module *procedure *" | sed 's/module //Ig' | tr '[:upper:]' '[:lower:]'`
  for modulename in $modulenamelist; do
    numModules=$((numModules + 1))
    echo $numModules $fullfilename $filename $modulename
    uses=`grep -A 10000000 -iE "^ *module *${modulename}" $fullfilename |grep -B 10000000 -iE "^ *end *module *${modulename}" |grep -i "^ *use *.*_mod\>" | sed 's/, *only *:.*//Ig' | sed 's/!.*//Ig' | sed 's/use //Ig' | tr '[:upper:]' '[:lower:]' | sort -u`
    echo $uses
    echo

    filenames[$numModules]=$filename
    fullfilenames[$numModules]=$fullfilename
    modulenames[$numModules]=$modulename
    useslist[$numModules]=$uses
    numberuses[$numModules]=`echo "${useslist[$numModules]}" |wc -w`
  done
  
done

echo "Finished scanning $numModules modules"

make_index modulename_index "${modulenames[@]}"
