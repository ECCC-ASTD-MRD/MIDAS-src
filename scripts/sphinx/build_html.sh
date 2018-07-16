#!/bin/bash

set -e

# SELECT WHICH FORTRAN SOURCE FILES TO INCLUDE IN DOCUMENTATION

codedir=${1:-../../src}
htmldir=${2:-~/public_html/midas_sphinx}

# ALL THE FILES (ONLY AFTER WORK IS DONE TO MODIFY ALL SOURCE FILES)
program_filelist=`ls -dR -1 $codedir/programs/*.f*90`
#module_filelist=`ls -dR -1 $codedir/modules/*.f*90 $codedir/modules/*/*.f*90`

# SMALL NUMBER OF FILES (JUST AS PROOF OF CONCEPT)
#program_filelist=`ls -dR -1 $codedir/programs/var.f90`
module_filelist=`ls -dR -1 $codedir/modules/utilities_mod.f90 $codedir/modules/mpi_mod.f90 $codedir/modules/mpivar_mod.f90`

# CREATE LINKS TO F90 FILES

rm -fR _src_files
mkdir _src_files

for file in $program_filelist ; do  
  echo ADDING THIS FILE TO src_files: $file
  cd _src_files
  ln -s ../$file ./
  cd ../
done
for file in $module_filelist ; do  
  echo ADDING THIS FILE TO src_files: $file
  cd _src_files
  ln -s ../$file ./
  cd ../
done

# GENERATE LIST OF ALL PROGRAMS

program_list=''
for file in $program_filelist ; do  
  program_name=`grep -i '^program ' $file | awk '{print $2}'`
  program_name_lc=`echo $program_name |tr '[:upper:]' '[:lower:]'`
  echo ADDING THIS PROGRAM TO THE LIST: $program_name_lc
  program_list="$program_list $program_name_lc"
done

# GENERATE LIST OF ALL MODULES

module_list=''
for file in $module_filelist ; do  
  module_name=`grep -i '^module ' $file | awk '{print $2}'`
  module_name_lc=`echo $module_name |tr '[:upper:]' '[:lower:]'`
  echo ADDING THIS MODULE TO THE LIST: $module_name_lc
  module_list="$module_list $module_name_lc"
done


# GENERATE THE RST FILES FOR THE MAIN PROGRAMS

rm -fR programs
mkdir programs
for program_name in $program_list ; do
cat > ./programs/${program_name}.rst <<EOF
==========================================
$program_name
==========================================

    .. f:autoprogram:: $program_name

EOF
done


# GENERATE THE RST FILES FOR THE MODULES

rm -fR modules
mkdir modules
for module_name in $module_list ; do
cat > ./modules/${module_name}.rst <<EOF
==========================================
$module_name
==========================================

    .. f:automodule:: $module_name

EOF
done


# GENERATE THE MAIN PAGE

cat > index.rst << 'EOF'
.. MIDAS documentation master file

Welcome to MIDAS documentation
==============================

This is the automatically generated MIDAS documentation. Below you
will find a list of all fortran programs and modules that make up
the MIDAS software. Documentation in the fortran code that appears 
in comments immediately following the program or module or subroutine
statement will be included. It can be formatted using *reStructuredText*.
A primer on this markup language can be found here:

http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html

.. note::  Currently only a few programs/modules are included since significant
           changes will be required in modifying existing comments to avoid sphinx
           from aborting. For example, since ``-------`` is used in *reStructuredText* 
           to indicate a heading, these must be removed from all comments.

Programs
========

.. toctree::
   :maxdepth: 2

EOF

for program_name in $program_list ; do
  echo "   programs/$program_name" >> index.rst
done

cat >> index.rst << 'EOF'

Modules
=======

.. toctree::
   :maxdepth: 2

EOF

for module_name in $module_list ; do
  echo "   modules/$module_name" >> index.rst
done

cat >> index.rst << 'EOF'

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

EOF


# RUN SPHYNX TO BUILD HTML FILES

rm -fR _build
mkdir _build
PYTHONPATH="$PYTHONPATH:${PWD}/lib/python2.7"
#sphinx-build -b html ./ ./_build
make html

# PUBLISH HTML FILES

[ -d "${htmldir}" ] && rm -rf ${htmldir}
mkdir -p ${htmldir}
mv _build/html/* ${htmldir}

echo "The HTML are in ${htmldir}"
# CLEAN UP

rm -fR _src_files
rm -fR _build
rm -fR *.rst programs modules
