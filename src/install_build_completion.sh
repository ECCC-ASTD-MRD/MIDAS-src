#!/bin/bash
set -x
typeset __here=$(pwd)
typeset __bc_dir=${HOME}/.bash_completion.d
typeset __check_line="## produced by MIDAS install_build_completion.sh"
typeset __profile_post=${HOME}/.profile.d/interactive/post

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
cp .build-completion.bash ${__bc_dir}/build_midas.bc

## Adding build_midas to PATH
if [ -z "$(cat ~/.profile.d/interactive/post  | grep PATH | grep $(pwd))" ]
then 
    cat >> ${__profile_post} << EOF
## Adding ${__here} to PATH for build_midas auto-completion
export PATH=\${PATH}:${__here}
EOF
    ## does not propagate outside of the install shell
    #source ${__profile_post}
fi

source ${HOME}/.bash_completion
set +x


echo
echo "Auto completion for build_midas installed"
echo 
echo "To use it directly:  \`source ${__profile_post}\`"
echo "(it will be automatically loaded on next shells)"
echo
