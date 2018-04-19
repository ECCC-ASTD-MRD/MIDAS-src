#!/bin/bash
#

for program in *.f90 ; do
    program_name=$(basename $program .f90)
    echo ''  1>> compile.${program_name}.listing.${ORDENV_PLAT} 2>&1
    echo " ------------------------------------------------ "  1>> compile.${program_name}.listing.${ORDENV_PLAT} 2>&1
    echo " >>> program $program_name "     1>> compile.${program_name}.listing.${ORDENV_PLAT} 2>&1
    echo " >>> program $program_name "
    ## Ask the question only if an interactive prompt
    if tty -s; then
        read -rsp $'     press enter to compile... \n'
    fi
    ./compile_program.sh $program_name 1>> compile.${program_name}.listing.${ORDENV_PLAT} 2>&1 &
done

wait

for program in *.f90 ; do
    program_name=$(basename $program .f90)
    cat compile.${program_name}.listing.${ORDENV_PLAT}
    rm compile.${program_name}.listing.${ORDENV_PLAT}
done
