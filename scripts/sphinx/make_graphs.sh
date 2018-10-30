#!/bin/bash

limit_levels=${1:-3}
GRAPHDIR=${2:-${PWD}/graphs/}
echo "ARGUMENTS = $limit_levels $GRAPHDIR "

make_modules="yes"
make_programs="yes"

# switch to the main source directory
ORIG_PWD=$PWD
SRCDIR=$PWD/../../src
cd $SRCDIR

. $ORIG_PWD/../prepare_dependencies.sh

if [ "${make_modules}" = "yes" ]; then

# prepare graphviz file for graphical representation of dependencies of each module
mkdir -p $GRAPHDIR/modules
rm -fR $GRAPHDIR/modules/*

for index1 in `seq 1 ${numModules}`; do
  dependencies_done=""
  all_modules=""
  modulename=${modulenames[$index1]}
  echo
  echo "GENERATING DEPENDENCY GRAPH FOR THE MODULE ${modulename}"
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
    if [ "${limit_levels}" -gt "1" ]; then
      for use3 in ${useslist[$index2]}; do
        index3=${modulename_index[$use3]}
        dependencies_done=`echo "${dependencies_done} ${use2}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
        all_modules=`echo "${all_modules} ${use3}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
        echo "${use2}->${use3};" >> $GRAPHDIR/modules/${modulename}.gv
        if [ "${limit_levels}" -gt "2" ]; then
          for use4 in ${useslist[$index3]}; do
            index4=${modulename_index[$use4]}
            dependencies_done=`echo "${dependencies_done} ${use3}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
            all_modules=`echo "${all_modules} ${use4}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
            echo "${use3}->${use4};" >> $GRAPHDIR/modules/${modulename}.gv
            if [ "${limit_levels}" -gt "3" ]; then
              for use5 in ${useslist[$index4]}; do
                index5=${modulename_index[$use5]}
                dependencies_done=`echo "${dependencies_done} ${use4}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
                all_modules=`echo "${all_modules} ${use5}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
                echo "${use4}->${use5};" >> $GRAPHDIR/modules/${modulename}.gv
              done
            fi
          done
        fi
      done
    fi
  done

  for use in ${all_modules}; do
    index=${modulename_index[$use]}
    if [ "${numberuses[$index]}" = "0" ]; then
      echo "${use} [style=filled];" >> $GRAPHDIR/modules/${modulename}.gv
    else
      if [[ ! "${dependencies_done}" =~ "${use}" ]]; then
        echo "${use} [color=red];" >> $GRAPHDIR/modules/${modulename}.gv
      fi
    fi
  done

  # finish the graph viz file
  echo "overlap=false" >> $GRAPHDIR/modules/${modulename}.gv
  echo "label=\"Read boxes indicate modules with dependencies not shown\nShaded boxes indicate modules with no dependencies\"" >> $GRAPHDIR/modules/${modulename}.gv
  echo "fontsize=14;" >> $GRAPHDIR/modules/${modulename}.gv
  echo "}" >> $GRAPHDIR/modules/${modulename}.gv
  unflatten -l 8 -f $GRAPHDIR/modules/${modulename}.gv > $GRAPHDIR/modules/${modulename}_2.gv
  dot -Tpng $GRAPHDIR/modules/${modulename}_2.gv > $GRAPHDIR/modules/${modulename}.png
  image_width=`file $GRAPHDIR/modules/${modulename}.png |cut -f2 -d',' |cut -f1 -d'x'`
  echo "The width of the image = $image_width"
  convert $GRAPHDIR/modules/${modulename}.png -shave 1x1 -bordercolor black -border 1 $GRAPHDIR/output.png
  mv $GRAPHDIR/output.png $GRAPHDIR/modules/${modulename}.png
done

fi # if make_modules


if [ "${make_programs}" = "yes" ]; then

cd ${SRCDIR}/programs
programfilelist=`ls *.f90 | tr '\n' ' '`
cd - >/dev/null

mkdir -p $GRAPHDIR/programs
rm -fR $GRAPHDIR/programs/*

for program in ${programfilelist}; do
  programname=`grep -i "^ *program *[a-zA-Z]" ${SRCDIR}/programs/$program | sed 's/program //Ig' | tr '[:upper:]' '[:lower:]'`

  dependencies_done=""
  all_modules=""

  echo
  echo "START BUILDING DEPENDENCY GRAPH FOR THE PROGRAM: $programname"

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
    if [ "${limit_levels}" -gt "1" ]; then
      for use2 in ${useslist[$index1]}; do
        index2=${modulename_index[$use2]}
        dependencies_done=`echo "${dependencies_done} ${use1}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
        all_modules=`echo "${all_modules} ${use2}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
        echo "${use1}->${use2};" >> $GRAPHDIR/programs/${programname}.gv
        if [ "${limit_levels}" -gt "2" ]; then
          for use3 in ${useslist[$index2]}; do
            index3=${modulename_index[$use3]}
            dependencies_done=`echo "${dependencies_done} ${use2}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
            all_modules=`echo "${all_modules} ${use3}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
            echo "${use2}->${use3};" >> $GRAPHDIR/programs/${programname}.gv
            if [ "${limit_levels}" -gt "3" ]; then
              for use4 in ${useslist[$index3]}; do
                index4=${modulename_index[$use4]}
                dependencies_done=`echo "${dependencies_done} ${use3}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
                all_modules=`echo "${all_modules} ${use4}" | tr ' ' '\n' | sort -u | tr '\n' ' '`
                echo "${use3}->${use4};" >> $GRAPHDIR/programs/${programname}.gv
              done
            fi
          done
        fi
      done
    fi
  done

  for use in ${all_modules}; do
    index=${modulename_index[$use]}
    if [ "${numberuses[$index]}" = "0" ]; then
      echo "${use} [style=filled];" >> $GRAPHDIR/programs/${programname}.gv
    else
      if [[ ! "${dependencies_done}" =~ "${use}" ]]; then
        echo "${use} [color=red];" >> $GRAPHDIR/programs/${programname}.gv
      fi
    fi
  done

  # finish the graph viz file
  echo "overlap=false" >> $GRAPHDIR/programs/${programname}.gv
  echo "label=\"Read boxes indicate modules with dependencies not shown\nShaded boxes indicate modules with no dependencies\"" >> $GRAPHDIR/programs/${programname}.gv
  echo "fontsize=14;" >> $GRAPHDIR/programs/${programname}.gv
  echo "}" >> $GRAPHDIR/programs/${programname}.gv
  unflatten -l 8 -f $GRAPHDIR/programs/${programname}.gv > $GRAPHDIR/programs/${programname}_2.gv
  dot -Tpng $GRAPHDIR/programs/${programname}_2.gv > $GRAPHDIR/programs/${programname}.png
  image_width=`file $GRAPHDIR/programs/${programname}.png |cut -f2 -d',' |cut -f1 -d'x'`
  echo "The width of the image = $image_width"
  convert $GRAPHDIR/programs/${programname}.png -shave 1x1 -bordercolor black -border 1 $GRAPHDIR/output.png
  mv $GRAPHDIR/output.png $GRAPHDIR/programs/${programname}.png

done

fi # if make_programs

echo
echo "FINISHED"
echo
