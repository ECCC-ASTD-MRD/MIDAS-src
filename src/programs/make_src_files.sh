#!/bin/bash

inputarg=$1

if [ "${inputarg}" = "all" -o "${inputarg}" = "" ]; then
  programlist=`ls *.f90 | tr '\n' ' '`
else
  programlist=$inputarg
fi

echo
echo "Will now generate src_files for the following list of programs:"
echo "    $programlist "
echo

# switch to the main source directory
ORIG_PWD=$PWD
cd ../

. $ORIG_PWD/../../scripts/prepare_dependencies.sh


for program in ${programlist}; do

programname=`echo ${program} |cut -f1 -d'.'`

echo
echo "START BUILDING DEPENDENCIES FOR THE PROGRAM: $programname"
echo

declare -a levels=("")
declare -a levelsfilename=("")

# build a list of modules used by the main program
uses1=`grep -i '^ *use *.*_mod' programs/$program | sed 's/, *only *:.*//Ig' | sed 's/!.*//Ig' | sed 's/use //Ig' | tr '[:upper:]' '[:lower:]' | sort -u`
for use1 in $uses1; do 
  echo $use1
  level=1; levels[$level]="${levels[$level]} $use1"
  index1=${modulename_index[$use1]}
  for use2 in ${useslist[$index1]}; do
    level=2; levels[$level]="${levels[$level]} $use2"
    index2=${modulename_index[$use2]}
    for use3 in ${useslist[$index2]}; do
      level=3; levels[$level]="${levels[$level]} $use3"
      index3=${modulename_index[$use3]}
      for use4 in ${useslist[$index3]}; do
        level=4; levels[$level]="${levels[$level]} $use4"
        index4=${modulename_index[$use4]}
        for use5 in ${useslist[$index4]}; do
          level=5; levels[$level]="${levels[$level]} $use5"
          index5=${modulename_index[$use5]}
          for use6 in ${useslist[$index5]}; do
            level=6; levels[$level]="${levels[$level]} $use6"
            index6=${modulename_index[$use6]}
            for use7 in ${useslist[$index6]}; do
              level=7; levels[$level]="${levels[$level]} $use7"
              index7=${modulename_index[$use7]}
              for use8 in ${useslist[$index7]}; do
                level=8; levels[$level]="${levels[$level]} $use8"
                index8=${modulename_index[$use8]}
                for use9 in ${useslist[$index8]}; do
                  level=9; levels[$level]="${levels[$level]} $use9"
                  index9=${modulename_index[$use9]}
                  for use10 in ${useslist[$index9]}; do
                    level=10; levels[$level]="${levels[$level]} $use10"
                    index10=${modulename_index[$use10]}
                    for use11 in ${useslist[$index10]}; do
                      level=11; levels[$level]="${levels[$level]} $use11"
                      index11=${modulename_index[$use11]}
                      for use12 in ${useslist[$index11]}; do
                        level=12; levels[$level]="${levels[$level]} $use12"
                        index12=${modulename_index[$use12]}
                        for use13 in ${useslist[$index12]}; do
                          level=13; levels[$level]="${levels[$level]} $use13"
                          index13=${modulename_index[$use13]}
                          for use14 in ${useslist[$index13]}; do
                            level=14; levels[$level]="${levels[$level]} $use14"
                            index14=${modulename_index[$use14]}
                            for use15 in ${useslist[$index14]}; do
                              level=15; levels[$level]="${levels[$level]} $use15"
                              index15=${modulename_index[$use15]}
                              for use16 in ${useslist[$index15]}; do
                                level=16; levels[$level]="${levels[$level]} $use16"
                                index16=${modulename_index[$use16]}
                                for use17 in ${useslist[$index16]}; do
                                  level=17; levels[$level]="${levels[$level]} $use17"
                                  index17=${modulename_index[$use17]}
                                  for use18 in ${useslist[$index17]}; do
                                    level=18; levels[$level]="${levels[$level]} $use18"
                                    index18=${modulename_index[$use18]}
                                    for use19 in ${useslist[$index18]}; do
                                    level=19; levels[$level]="${levels[$level]} $use19"
                                    echo " ERROR: More than 19 levels of dependencies found in program ${program}"
                                    echo " ERROR: to generate src_files for this program, increase the"
                                    echo " ERROR: number of levels in the scripts."
                                    exit
				  done
                                 done
                                done
                              done
                            done
                          done
                        done
                      done
                    done
                  done
                done
              done
            done
          done
        done
      done
    done
  done
done

echo
echo "Used modules at each level before first pruning"
echo

# remove duplicate modules at each level
numlevels=${#levels[@]}
numlevels=$((numlevels - 1))
for level in `seq $numlevels -1 1`; do
  levels[$level]=`echo ${levels[level]} | tr ' ' '\n' | sort -u | tr '\n' ' '`
  echo "level $level = ${levels[$level]}"
done

# prune to ensure modules appearing at multiple levels only remain at the lowest level
for level1 in `seq $numlevels -1 1`; do
  for modulename in ${levels[$level1]}; do
    #echo "pruning ${modulename}"
    for level2 in `seq $((level1 - 1)) -1 1`; do
      levels[$level2]=`echo "${levels[$level2]}" | sed "s/\<${modulename}\>//g"`
    done
  done
done

echo
echo "Used modules at each level after first pruning"
echo

for level in `seq $numlevels -1 1`; do
  levels[$level]=`echo ${levels[level]} | tr ' ' '\n' | sort -u | tr '\n' ' '`
  echo "level $level = ${levels[$level]}"
done

echo
echo "Used module filenames at each level before second pruning"
echo

# generate levels for filenames
for level in `seq $numlevels -1 1`; do
  for modulename in ${levels[$level]}; do
    index=${modulename_index[$modulename]}
    levelsfilename[$level]="${levelsfilename[$level]} ${filenames[index]}"
  done
  echo "level $level = ${levelsfilename[$level]}"
done

# prune in filenames opposite direction (needed for files containing with multiple modules)
for level1 in `seq 1 $numlevels`; do
  for filename in ${levelsfilename[$level1]}; do
    #echo "pruning ${filename}"
    for level2 in `seq $((level1 + 1)) $numlevels`; do
      levelsfilename[$level2]=`echo "${levelsfilename[$level2]}" | sed "s/\<${filename}\>//g"`
    done
  done
done

echo
echo "Used module filenames at each level after second pruning"
echo

for level in `seq $numlevels -1 1`; do
  echo "level $level = ${levelsfilename[$level]}"
done

# generate src_files file
srcfilesname="src_files_${programname}.sh"
echo 'SRC_FILES=""' > $ORIG_PWD/src_files/$srcfilesname
for level in `seq $numlevels -1 1`; do
  levels[$level]=`echo ${levels[level]} | tr ' ' '\n' | sort -u | tr '\n' ' '`
  echotext='SRC_FILES="$SRC_FILES '
  echotext="$echotext ${levelsfilename[$level]}"
  echotext="$echotext\""
  echo $echotext >> $ORIG_PWD/src_files/$srcfilesname
done

echo "DONE BUILDING DEPENDENCIES FOR THE PROGRAM: $program"

done
