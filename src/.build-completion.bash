#/usr/bin/env bash

_build_completions()
{
    COMPREPLY+=($(compgen -W "$(cat Makefile | grep '^[1-9a-z_.]*:' | cut -d':' -f1)" -- "${COMP_WORDS[COMP_CWORD]}"))
    COMPREPLY+=($(compgen -W "$(ls programs/*.f90 |sed 's/programs\/\(.*\).f90/\1.Abs/')" -- "${COMP_WORDS[COMP_CWORD]}"))
    COMPREPLY+=($(compgen -W "$(ls programs/*.f90 |sed 's/programs\/\(.*\).f90/\1.o/')" -- "${COMP_WORDS[COMP_CWORD]}"))
    COMPREPLY+=($(compgen -W "$(ls modules/*.f90 |sed 's/modules\/\(.*\).f90/\1.o/')" -- "${COMP_WORDS[COMP_CWORD]}"))

}
complete -F _build_completions ./midas_build
