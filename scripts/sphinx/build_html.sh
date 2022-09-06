#!/bin/bash

set -e

# ARGUMENTS SPECIFY DIRECTORIES OF MIDAS SOURCE CODE AND DESTINATION

codedir=${1:-../../src}
htmldir=${2:-~/public_html/midas_sphinx}

toplevel=$(git rev-parse --show-toplevel)

# CHOOSE WHETHER OR NOT TO GENERATE DEPENDENCY GRAPHS (COSTLY) AND NAMELIST INFORMATION

do_graphs=yes
do_namelists=yes

# PREPARE THE MODULE DEPENDENCY ARRAYS

ORIG_PWD=$PWD
SRCDIR=$PWD/../../src
cd $SRCDIR
. $ORIG_PWD/../prepare_dependencies.sh
cd $ORIG_PWD

# SELECT WHICH FORTRAN SOURCE FILES TO INCLUDE IN DOCUMENTATION

# GENERATE LIST OF ALL PROGRAMS

program_filelist=`ls -dR -1 $codedir/programs/*.f*90`
#program_filelist=""
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
module_category[2]="B and R matrices"
module_category[3]="Observation input/output"
module_category[4]="Data object transformations"
module_category[5]="Observation operators"
module_category[6]="High-level data objects"
module_category[7]="Low-level data objects"
module_category[8]="Low-level utilities and constants"
module_category[9]="Global interfaces"
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

    **Dependency Diagrams:**

    .. figure:: /${program_names[$filenum]}.svg
        :height: 100px

        Direct Dependency Diagram

EOF
fi

cat >> ./programs/${program_names[$filenum]}.rst <<EOF

    .. f:autoprogram:: ${program_names[$filenum]}

EOF

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

    **Dependency Diagrams:**

    .. figure:: /${module_name}.svg
        :height: 100px

        Direct Dependency Diagram

    .. figure:: /${module_name}_rev.svg
        :height: 100px

        Reverse Dependency Diagram

EOF

fi
cat >> ./modules/${module_name}.rst <<EOF

    .. f:automodule:: ${module_name}

EOF

done

revision=$(${toplevel}/midas.version.sh)

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

cat >> index.rst << 'EOF'

Additional information
======================

* Namelists: Information is :doc:`here on which namelists are used for each program. <namelists_in_each_program>`

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
tmgstrings=`grep -ir 'call *utl_tmg_start(' | sed "s/['\"]//g" | sed 's/\!//g' | sed 's/ //g' | sed 's/:.*(/,/g' | sed 's/)//g' | awk -F',' '{ print $2 "," $1 "," $3}' | sort -u | sort -n`
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

# GENERATE LIST OF NAMELISTS USED FOR EACH PROGRAM

# Build dependency tree
echo "Building dependency tree"
echo "Sourcing config"
cd ${toplevel}/src
source ./config.dot.sh
GITDESC=$(git describe --abbrev=7 --dirty=_M)
OBJBLD_PATH=${toplevel}/compiledir/midas_bld-${GITDESC}/${ARCH}/
make depend

cat > namelists_in_each_program.rst << 'EOF'

Namelists (possibly) used in each MIDAS program
===============================================

EOF

# switch to programs directory
cd ${SRCDIR}/programs
programList=`ls *.f90 | tr '\n' ' '`
# switch to the src_files directory
PGM_DIR=$PWD
cd ../modules
MOD_DIR=$PWD
cd ${OBJBLD_PATH}
echo "programList = ${programList}"
for programName in ${programList}; do
    echo 
    echo "Generating list of namelists for the program: ${programName}"
    programNameBase=$(basename ${programName} .f90)
    echo "========================================" >> ${ORIG_PWD}/namelists_in_each_program.rst
    echo "${programNameBase}" >> ${ORIG_PWD}/namelists_in_each_program.rst
    echo "========================================" >> ${ORIG_PWD}/namelists_in_each_program.rst
    src_files=`cat dep.obj.inc | grep "^${programNameBase}.o :" | sed "s/^${programNameBase}.o : ${programNameBase}\.[^ ]* //"`
    src_files2=`echo ${src_files} | sed "s/\(\w\)\.o/\1/g" | tr " " "\n"`
    # finding proper suffix (.f90 or .ftn90)
    src_files3=""
    for src_file in ${src_files2}; do
        src_files3="${src_files3} `ls ${MOD_DIR} | grep ${src_file}`"
    done
    nameListList=`grep -i 'namelist */ *[a-z0-9]* */' ${PGM_DIR}/${programName} \
                  | cut -f2 -d '/' | tr '[:upper:]' '[:lower:]' \
                  | grep -vi ptopo` || true
    for src_file in ${src_files3}; do
        nameListList="${nameListList} `grep -i 'namelist */ *[a-z0-9]* */' ${MOD_DIR}/${src_file} \
                                       | cut -f2 -d '/' | tr '[:upper:]' '[:lower:]' \
                                       | grep -vi ptopo`" || true
    done
    nameListList=`echo $nameListList | tr " " "\n" | sort -u`
    echo $nameListList
    for nameList in ${nameListList}; do
      echo "- ${nameList}" >> ${ORIG_PWD}/namelists_in_each_program.rst
    done
    echo >> ${ORIG_PWD}/namelists_in_each_program.rst
done
echo "========================================" >> ${ORIG_PWD}/namelists_in_each_program.rst
cd $ORIG_PWD

rm -fR _build
mkdir _build

# RUN SPHYNX TO BUILD HTML FILES

PYTHONPATH="$PYTHONPATH:${toplevel}/tools/sphinx-fortran"
#sphinx-build -b html ./ ./_build
make html

# PUBLISH HTML FILES

[ -d "${htmldir}" ] && rm -rf ${htmldir}
mkdir -p ${htmldir}
mv _build/html/* ${htmldir}

# GENERATE DEPENDENCY GRAPHS

if [ "${do_graphs}" = "yes" ]; then
  ./make_graphs.sh $PWD/_build/html/graphs
  mkdir -p ${htmldir}/programs
  mv $PWD/_build/html/graphs/programs/*.svg ${htmldir}/programs
  mkdir -p ${htmldir}/modules
  mv $PWD/_build/html/graphs/modules/*.svg ${htmldir}/modules
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
rm -fR *.rst programs modules namelist_listing.txt
rm -fR graphs
