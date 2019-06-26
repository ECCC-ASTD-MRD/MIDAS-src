#!/bin/bash

set -e

# ARGUMENTS SPECIFY DIRECTORIES OF MIDAS SOURCE CODE AND DESTINATION

codedir=${1:-../../src}
htmldir=${2:-~/public_html/midas_sphinx}

# CHOOSE WHETHER OR NOT TO GENERATE DEPENDENCY GRAPHS (COSTLY) AND NAMELIST INFORMATION

do_graphs=no
do_namelists=no

# PREPARE THE MODULE DEPENDENCY ARRAYS

ORIG_PWD=$PWD
SRCDIR=$PWD/../../src
cd $SRCDIR
. $ORIG_PWD/../prepare_dependencies.sh
cd $ORIG_PWD

# SELECT WHICH FORTRAN SOURCE FILES TO INCLUDE IN DOCUMENTATION

# GENERATE LIST OF ALL PROGRAMS

#program_filelist=`ls -dR -1 $codedir/programs/*.f*90`
program_filelist=""
numPrograms=0
for file in $program_filelist ; do  
  numPrograms=$((numPrograms + 1))
  program_name=`grep -i '^program ' $file | awk '{print $2}'`
  program_name_lc=`echo $program_name |tr '[:upper:]' '[:lower:]'`
  echo ADDING THIS PROGRAM TO THE LIST: $program_name_lc
  program_names[$numPrograms]="$program_name_lc"
  program_files[$numPrograms]="$file"
  program_bfiles[$numPrograms]=`basename "$file"`
done
echo "Number of programs = $numPrograms"

module_filelist=`ls -dR -1 $codedir/modules/*f*90`

# DEFINE THE MODULE CATEGORY NAMES FOR EACH NUMERICAL CODE

module_category[1]="High-level functionality"
module_category[2]="High-level data objects"
module_category[3]="High-level transformations"
module_category[4]="Observation operators"
module_category[5]="B and R matrices"
module_category[6]="Observation input/output"
module_category[7]="Low-level data objects and utilities"
module_category[8]="Global constants and interfaces"
num_categories=${#module_category[@]}

# CREATE LINKS TO F90 FILES

rm -fR _src_files
mkdir _src_files

# FIRST CREATE DUMMY FORTRAN MODULES FOR ALL THAT HAVE NOT YET BEEN MODIFIED
for index in `seq 1 $numModules`; do
cat >> _src_files/${filenames[$index]} <<EOF
module ${modulenames[$index]}
! MODULE ${modulenames[$index]} (prefix='${prefixes[$index]}' category='${categories[$index]}. ${module_category[${categories[$index]}]}')
!
! **Purpose:** Temporary place holder for modules that have not yet been modified
end module ${modulenames[$index]}
EOF
done

for filenum in `seq 1 $numPrograms` ; do  
  echo ADDING THIS FILE TO src_files: ${program_files[$filenum]}
  cd _src_files
  ln -s ../${program_files[$filenum]} ./
  cd ../
done
for file in $module_filelist ; do  
  echo ADDING THIS FILE TO src_files: $file
  cd _src_files
  bname=`basename $file`
  rm -f $bname
  ln -s ../$file ./
  cd ../
done

rm -fR programs
mkdir programs
rm -fR modules
mkdir modules

# GENERATE THE RST FILES FOR THE MAIN PROGRAMS
for filenum in `seq 1 $numPrograms` ; do

cat > ./programs/${program_names[$filenum]}_src.rst <<EOF
==========================================
${program_names[$filenum]} source
==========================================

    .. literalinclude:: ../_src_files/${program_bfiles[$filenum]}
       :language: fortran
       :linenos:

EOF


cat > ./programs/${program_names[$filenum]}.rst <<EOF
==========================================
${program_names[$filenum]}
==========================================

    :doc:\`link to source code <${program_names[$filenum]}_src>\`

EOF

if [ "${do_graphs}" = "yes" ]; then
cat >> ./programs/${program_names[$filenum]}.rst <<EOF
    **Dependency Diagrams:** \`1-Level <level1/${program_names[$filenum]}.png>\`_, \`2-Level <level2/${program_names[$filenum]}.png>\`_, \`3-Level <level3/${program_names[$filenum]}.png>\`_

EOF
fi

cat >> ./programs/${program_names[$filenum]}.rst <<EOF

    .. f:autoprogram:: ${program_names[$filenum]}

EOF

if [ "${do_graphs}" = "yes" ]; then
cat >> ./programs/${program_names[$filenum]}.rst <<EOF

    **Dependency Diagrams:**

    .. figure:: /level1/${program_names[$filenum]}.png
        :height: 100px

        1-Level Dependency Diagram

    .. figure:: /level2/${program_names[$filenum]}.png
        :height: 100px

        2-Level Dependency Diagram

    .. figure:: /level3/${program_names[$filenum]}.png
        :height: 100px

        3-Level Dependency Diagram

EOF
fi

done


# GENERATE THE RST FILES FOR THE MODULES

for index in `seq 1 $numModules`; do
module_name=${modulenames[$index]}

cat > ./modules/${module_name}_src.rst <<EOF
==========================================
${module_name} source
==========================================

    .. literalinclude:: ../_src_files/${filenames[$index]}
       :language: fortran
       :linenos:

EOF

cat > ./modules/${module_name}.rst <<EOF
==========================================
$module_name
==========================================

    :doc:\`link to source code <${module_name}_src>\`

EOF

if [ "${do_graphs}" = "yes" ]; then
cat >> ./modules/${module_name}.rst <<EOF

    **Dependency Diagrams:** \`1-Level <level1/${module_name}.png>\`_, \`2-Level <level2/${module_name}.png>\`_, \`3-Level <level3/${module_name}.png>\`_
EOF
fi

cat >> ./modules/${module_name}.rst <<EOF

    .. f:automodule:: ${module_name}

EOF

if [ "${do_graphs}" = "yes" ]; then
cat >> ./modules/${module_name}.rst <<EOF

    **Dependency Diagrams:**

    .. figure:: /level1/${module_name}.png
        :height: 100px

        1-Level Dependency Diagram

    .. figure:: /level2/${module_name}.png
        :height: 100px

        2-Level Dependency Diagram

    .. figure:: /level3/${module_name}.png
        :height: 100px

        3-Level Dependency Diagram

EOF
fi

done

revision=$(git describe --always --dirty=_M 2>/dev/null)

# GENERATE THE MAIN PAGE

cat > index.rst << EOF
.. MIDAS documentation master file

Welcome to MIDAS documentation
==============================

This is the documentation for version ${revision}

This is the automatically generated MIDAS documentation. Below you
will find a list of all fortran programs and modules that make up
the MIDAS software. Documentation in the fortran code that appears 
in comments immediately following the program or module or subroutine
statement will be included. It can be formatted using *reStructuredText*.
A primer on this markup language can be found here:

http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html

https://matplotlib.org/sampledoc/cheatsheet.html

EOF

cat >> index.rst << EOF

Programs
========

.. toctree::
   :maxdepth: 2

EOF

for filenum in `seq 1 $numPrograms` ; do
  echo "   programs/${program_names[filenum]}" >> index.rst
done


for category_index in `seq 1 $num_categories`; do

cat >> index.rst << EOF

Modules: ${category_index}. ${module_category[$category_index]}
===========================================================================

.. toctree::
   :maxdepth: 2

EOF

module_list=''
for index in `seq 1 $numModules`; do
  module_list="${module_list} ${modulenames[$index]}"
done
module_list_sort=`echo $module_list | tr ' ' '\n' | sort | tr '\n' ' '`

for module_name in ${module_list_sort} ; do
  index=${modulename_index[$module_name]}
  if [ "${categories[$index]}" = "${category_index}" ]; then
    echo "   modules/$module_name" >> index.rst
  fi
done

done # category_index

cat >> index.rst << EOF

Library
========

A library is published and can be accessed with::
   . ssmuse-sh -d eccc/mrd/rpn/anl/midas/$(git describe --abbrev=0 | sed 's/v_//')

This library is compiled with \`\`CODEPRECISION_OBS_REAL_SINGLE\`\` defined.

This library contains the following set of modules:
EOF

modList=$($(dirname $(true_path $0))/buildList.sh obsIO.f90)
for module_name in ${modList}; do
    modPath=$(git grep -li --full-name "^ *module *${module_name}" $(git rev-parse --show-cdup)/src/modules | grep -v unit_tests || true)
    if [ -n "${modPath}" ]; then
        str2parse="$(echo ${module_name} | sed 's/_/\\_/g') <https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/blob/$(git describe)/${modPath}/>"
        #echo " * \`$(echo ${str2parse} | sed 's/_/\\_/g')\`_" >> index.rst
        echo " * \`${str2parse}\`_" >> index.rst
    fi
done

cat >> index.rst << 'EOF'

Additional information
======================

* Namelists: Information is `here on the definition of all namelists. <namelists.html>`_

* TMG timing blocks: Information is :doc:`here on the numbering and labelling of all TMG timing blocks. <tmg_information>`

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

EOF

# GENERATE TMG TIMING BLOCK INFORMATION

cat > tmg_information.rst << 'EOF'

TMG Timing block information
============================

======================================== ============= =================================
Filename                                 Block Number  Block Label
======================================== ============= =================================
EOF

cd $SRCDIR
tmgstrings=`grep -ir 'tmg_start(' | sed "s/['\"]//g" | sed 's/\!//g' | sed 's/ //g' | sed 's/:.*(/,/g' | sed 's/)//g' | awk -F',' '{ print $2 "," $1 "," $3}' | sort -u | sort -n`
for tmgstring in $tmgstrings; do
  filename=`echo $tmgstring |cut -f2 -d','`
  filename=`basename $filename`
  if (( ${#filename} > 40 )); then
    filename2="${filename:0:40}"
  else
    filename2=$filename
  fi
  filename_fmt=`printf "%-40s" "$filename2"`
  tmgnumber=`echo $tmgstring |cut -f1 -d','`
  tmgnumber_fmt=`printf "%-13s" "$tmgnumber"`
  tmglabel=`echo $tmgstring |cut -f3 -d','`
  tmglabel_fmt=`printf "%-30s" "$tmglabel"`
  echo "$filename_fmt $tmgnumber_fmt $tmglabel_fmt" >> ${ORIG_PWD}/tmg_information.rst
done
cd $ORIG_PWD

echo "======================================== ============= =================================" >> tmg_information.rst


rm -fR _build
mkdir _build

# RUN SPHYNX TO BUILD HTML FILES

PYTHONPATH="$PYTHONPATH:${PWD}/lib/python2.7:${PWD}/lib/python2.7/sphinx_fortran-1.0.1-py2.7.egg"
#sphinx-build -b html ./ ./_build
make html

# PUBLISH HTML FILES

[ -d "${htmldir}" ] && rm -rf ${htmldir}
mkdir -p ${htmldir}
mv _build/html/* ${htmldir}

# GENERATE DEPENDENCY GRAPHS

if [ "${do_graphs}" = "yes" ]; then
  ./make_graphs.sh 1 $PWD/_build/html/graphs
  mkdir -p ${htmldir}/programs/level1
  mv $PWD/_build/html/graphs/programs/*.png ${htmldir}/programs/level1/
  mkdir -p ${htmldir}/modules/level1
  mv $PWD/_build/html/graphs/modules/*.png ${htmldir}/modules/level1/

  ./make_graphs.sh 2 $PWD/_build/html/graphs
  mkdir -p ${htmldir}/programs/level2
  mv $PWD/_build/html/graphs/programs/*.png ${htmldir}/programs/level2/
  mkdir -p ${htmldir}/modules/level2
  mv $PWD/_build/html/graphs/modules/*.png ${htmldir}/modules/level2/

  ./make_graphs.sh 3 $PWD/_build/html/graphs
  mkdir -p ${htmldir}/programs/level3
  mv $PWD/_build/html/graphs/programs/*.png ${htmldir}/programs/level3/
  mkdir -p ${htmldir}/modules/level3
  mv $PWD/_build/html/graphs/modules/*.png ${htmldir}/modules/level3/
fi

# GENERATE NAMELIST INFORMATION

if [ "${do_namelists}" = "yes" ]; then
  ./make_namelists.sh ${codedir} ${htmldir} > namelist_listing.txt
fi

echo
echo "The HTML are in ${htmldir}"

# CLEAN UP

rm -fR _src_files
rm -fR _build
rm -fR *.rst programs modules
rm -fR graphs
