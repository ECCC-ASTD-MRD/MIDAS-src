#!/bin/bash

GRAPHDIR=${1:-${PWD}/graphs/}
echo
echo "make_graphs.sh called with argument = $GRAPHDIR "

# some configuration variables
make_modules="yes"
make_programs="yes"
verbose="no"


commonLabel="Red boxes indicate modules with dependencies not shown\nLower level dependencies can be shown by clicking on a red box\nShaded boxes indicate modules with no dependencies"
moduleLabel=${commonLabel}
programLabel=${commonLabel}

# switch to the main source directory
ORIG_PWD=$PWD
SRCDIR=$PWD/../../src
cd $SRCDIR

. $ORIG_PWD/../prepare_dependencies.sh

if [ "${make_modules}" = "yes" ]; then

# prepare graphviz file for graphical representation of dependencies of each module
mkdir -p $GRAPHDIR/modules
rm -fR $GRAPHDIR/modules/*

echo "GENERATING DEPENDENCY GRAPHS FOR ALL MODULES"
for index1 in `seq 1 ${numModules}`; do
  dependencies_done=""
  all_modules=""
  modulename=${modulenames[$index1]}
  [ "${verbose}" = "yes" ] && echo
  [ "${verbose}" = "yes" ] && echo "GENERATING DEPENDENCY GRAPH FOR THE MODULE ${modulename}"
  echo "strict digraph ${modulename} {" > $GRAPHDIR/modules/${modulename}.gv
  echo "node [shape=box];" >> $GRAPHDIR/modules/${modulename}.gv
  echo "${modulename};" >> $GRAPHDIR/modules/${modulename}.gv
  if [ "${numberuses[$index1]}" = "0" ]; then
    echo "${modulename} [style=filled];" >> $GRAPHDIR/modules/${modulename}.gv
  fi
  for use2 in ${useslist[$index1]}; do
    index2=${modulename_index[$use2]}
    all_modules=`echo "${all_modules} ${use2}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
    echo "${modulename}->${use2};" >> $GRAPHDIR/modules/${modulename}.gv
  done

  for use in ${all_modules}; do
    index=${modulename_index[$use]}
    if [ "${numberuses[$index]}" = "0" ]; then
      echo "${use} [style=filled];" >> $GRAPHDIR/modules/${modulename}.gv
    else
      if [[ ! "${dependencies_done}" =~ "${use}" ]]; then
        echo "${use} [color=red URL=\"../modules/${use}.svg\"];" >> $GRAPHDIR/modules/${modulename}.gv
      fi
    fi
  done

  # finish the graph viz file
  echo "overlap=false" >> $GRAPHDIR/modules/${modulename}.gv
  echo "label=\"${moduleLabel}\"" >> $GRAPHDIR/modules/${modulename}.gv
  echo "fontsize=14;" >> $GRAPHDIR/modules/${modulename}.gv
  echo "}" >> $GRAPHDIR/modules/${modulename}.gv
  unflatten -l 8 -f $GRAPHDIR/modules/${modulename}.gv > $GRAPHDIR/modules/${modulename}_2.gv
  dot -Tsvg $GRAPHDIR/modules/${modulename}_2.gv > $GRAPHDIR/modules/${modulename}.svg
done

fi # if make_modules


if [ "${make_programs}" = "yes" ]; then

cd ${SRCDIR}/programs
programfilelist=`ls *.f90 | tr '\n' ' '`
cd - >/dev/null

mkdir -p $GRAPHDIR/programs
rm -fR $GRAPHDIR/programs/*

echo "GENERATING DEPENDENCY GRAPHS FOR ALL PROGRAMS"
for program in ${programfilelist}; do
  programname=`grep -i "^ *program *[a-zA-Z]" ${SRCDIR}/programs/$program | sed 's/program //Ig' | tr '[:upper:]' '[:lower:]'`

  dependencies_done=""
  all_modules=""

  [ "${verbose}" = "yes" ] && echo
  [ "${verbose}" = "yes" ] && echo "START BUILDING DEPENDENCY GRAPH FOR THE PROGRAM: $programname"

  # prepare graphviz file for graphical representation of dependencies of each program
  echo "strict digraph ${programname} {" > $GRAPHDIR/programs/${programname}.gv
  echo "node [shape=box];" >> $GRAPHDIR/programs/${programname}.gv
  echo "${programname};" >> $GRAPHDIR/programs/${programname}.gv

  # build a list of modules used by the main program
  uses1=`grep -i '^ *use *.*_mod' ${SRCDIR}/programs/$program | sed 's/, *only *:.*//Ig' | sed 's/!.*//Ig' | sed 's/use //Ig' | tr '[:upper:]' '[:lower:]' | sort -u`
  for use1 in $uses1; do 
    index1=${modulename_index[$use1]}
    all_modules=`echo "${all_modules} ${use1}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
    echo "${programname}->${use1};" >> $GRAPHDIR/programs/${programname}.gv
  done

  for use in ${all_modules}; do
    index=${modulename_index[$use]}
    if [ "${numberuses[$index]}" = "0" ]; then
      echo "${use} [style=filled];" >> $GRAPHDIR/programs/${programname}.gv
    else
      if [[ ! "${dependencies_done}" =~ "${use}" ]]; then
        echo "${use} [color=red URL=\"../modules/${use}.svg\"];" >> $GRAPHDIR/programs/${programname}.gv
      fi
    fi
  done

  # finish the graph viz file
  echo "overlap=false" >> $GRAPHDIR/programs/${programname}.gv
  echo "label=\"${programLabel}\"" >> $GRAPHDIR/programs/${programname}.gv
  echo "fontsize=14;" >> $GRAPHDIR/programs/${programname}.gv
  echo "}" >> $GRAPHDIR/programs/${programname}.gv
  unflatten -l 8 -f $GRAPHDIR/programs/${programname}.gv > $GRAPHDIR/programs/${programname}_2.gv
  dot -Tsvg $GRAPHDIR/programs/${programname}_2.gv > $GRAPHDIR/programs/${programname}.svg

done

fi # if make_programs
