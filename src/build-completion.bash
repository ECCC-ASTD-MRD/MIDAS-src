#/usr/bin/env bash
_build_completions()
{
    if [ "${#COMP_WORDS[@]}" != "2" ]; then
        return
    fi

    
    COMPREPLY+=($(compgen -W "$(cat Makefile | grep '^[1-9a-z_]*:' | cut -d':' -f1)" -- "${COMP_WORDS[1]}"))
    COMPREPLY+=($(compgen -W "$(ls programs/*.f90 |sed 's/programs\/\(.*\).f90/\1.Abs/')" -- "${COMP_WORDS[1]}"))
    COMPREPLY+=($(compgen -W "$(ls programs/*.f90 |sed 's/programs\/\(.*\).f90/\1.o/')" -- "${COMP_WORDS[1]}"))
    COMPREPLY+=($(compgen -W "$(ls modules/*.f90 |sed 's/modules\/\(.*\)_mod.f90/\1.o/')" -- "${COMP_WORDS[1]}"))

}
complete -F _build_completions ./build.sh
