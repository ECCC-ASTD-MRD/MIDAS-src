#!/bin/bash
#

for program in *.f90 ; do
    program_name=$(basename $program .f90)
    echo ''
    echo " ------------------------------------------------ "
    echo " >>> program $program_name "    
    read -rsp $'     press enter to compile... \n'
    ./compile_program.sh $program_name
done
