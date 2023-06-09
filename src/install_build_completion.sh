#!/bin/bash

###############################################################################
##  `midas_build` auto-completion install
##
##  Provide bash auto-completion for `midas_build`.
##  Need to be installed once; modify
##      * `${HOME}/.bash_completion`
##      * `${HOME}/.bash_completion.d/`
##
##  @SYNOPSIS
##      ```
##      ./install_build_completion.sh
##      source ${HOME}/.bash_completion
##      ```
##
##  @REFERENCES 
##      * `./README.md`
##  
##  @AUTHORS
##      Martin Deshaies-Jacques (@mad001) - CMDA - December 2020
##
###############################################################################

set -x
typeset __here=$(pwd)
typeset __bc_dir=${HOME}/.bash_completion.d
typeset __check_line="## produced by MIDAS install_build_completion.sh"

## Check existing ~/.bash_completion
if [ -f ${HOME}/.bash_completion ]
then
    ## check if produced by this install
    if [ ! "$(cat ${HOME}/.bash_completion | head -n1)" == "${__check_line}" ]
    then
        echo "${HOME}/.bash_completion exist and was not produced by this install script"
        echo "You might want to move this script inside ${__bc_dir}."
        exit 1
    fi
fi
echo ${__check_line} > ${HOME}/.bash_completion
cat >>  ${HOME}/.bash_completion << EOF
for bcfile in ~/.bash_completion.d/*
do
    [ -f "\$bcfile" ] && . \$bcfile
done
EOF

## Create ~/.bash_completion.d and copy bc script
mkdir -p ${__bc_dir}
cp .build-completion.bash ${__bc_dir}/midas_build.bc

source ${HOME}/.bash_completion
set +x


echo
echo "Auto completion for midas_build installed"
echo 
echo "To use it directly (in the present shell):"
echo '   source ${HOME}/.bash_completion'
echo "(in any case, it will be automatically loaded on next shells)"
echo
