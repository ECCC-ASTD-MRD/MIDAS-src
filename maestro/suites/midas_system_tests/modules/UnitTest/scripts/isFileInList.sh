#!/bin/bash

set -e

program_name=$0

print_usage () {
    set -e

    if [ $# -ne 0 ]; then
        echo "isFileInList.print_usage: This function takes no argument!" >&2
        return 1
    fi

    echo "isFileInList:"
    echo "       takes the 1 first argument and check if it is contained in the list"
    echo "         ${program_name} one in one two three"
    echo "       will return 0 since 'one' is contained in the list 'one two three'"
    echo "       and this command"
    echo "         ${program_name} four in one two three"
    echo "       will return 1 since 'four' is not contained in the list 'one two three'."
    echo "       "
    echo "       It can also check if it is not in the list:"
    echo "         ${program_name} four not in one two three"
    echo "       will return 0 since 'four' is not contained in the list 'one two three'"
    echo "       and this command"
    echo "         ${program_name} four not in one two three"
    echo "       will return 1 since 'four' is not contained in the list 'one two three'."
    echo
    echo "optional arguments:"
    echo "   -h, --help       show this help message and exit"
    echo "   -u, --unittests  run the unit tests"
}

isFileInList () {
    set -e

    if [ $# -lt 3 ]; then
        echo "isFileInList.isFileInList: need at least 3 arguments!" >&2
        echo >&2
        print_usage >&2
        return 1
    fi

    __isFileInList_file__=${1}

    if [ "${2}" = "in" ]; then
        shift 2
    elif [ "${2}" = "not" -a "${3}" = "in" ]; then
        shift 3
        ## If call with  'not in', the we call with 'in' and inverse the result
        __isFileInList_status__=0
        isFileInList "${__isFileInList_file__}" in $* || __isFileInList_status__=1
        unset __isFileInList_file__
        if [ "${__isFileInList_status__}" -eq 0 ]; then
            unset __isFileInList_status__
            return 1
        else
            unset __isFileInList_status__
            return 0
        fi
    else
        echo "isFileInList.isFileInList: Must be called with the form 'elem in \${list}' or 'elem not in \${list}'" >&2
        echo "isFileInList.isFileInList: You called with ${program_name} $*" >&2
        echo >&2
        print_usage >&2
        return 1
    fi

    __isFileInList_status__=0
    echo $* | grep -sq "\b${__isFileInList_file__}\b" || __isFileInList_status__=1
    unset __isFileInList_file__

    return ${__isFileInList_status__}
} ## End of function 'isFileInList'


unittests () {
    set -e

    global_status=0
    status=0
    isFileInList allo in allo bonjour salut || status=1
    [ "${status}" -eq 0 ] || echo "Test 1 unsuccessful"
    [ "${status}" -eq 0 ] || global_status=1

    status=0
    isFileInList bonjour in allo bonjour salut || status=1
    [ "${status}" -eq 0 ] || echo "Test 2 unsuccessful"
    [ "${status}" -eq 0 ] || global_status=1

    status=0
    isFileInList salut in allo bonjour salut || status=1
    [ "${status}" -eq 0 ] || echo "Test 3 unsuccessful"
    [ "${status}" -eq 0 ] || global_status=1

    status=0
    isFileInList hello in allo bonjour salut || status=1
    [ "${status}" -ne 0 ] || echo "Test 4 unsuccessful"
    [ "${status}" -eq 0 ] || global_status=1

    status=0
    isFileInList allo not in allo bonjour salut || status=1
    [ "${status}" -ne 0 ] || echo "Test 5 unsuccessful"
    [ "${status}" -eq 0 ] || global_status=1

    status=0
    isFileInList bonjour not in allo bonjour salut || status=1
    [ "${status}" -ne 0 ] || echo "Test 6 unsuccessful"
    [ "${status}" -eq 0 ] || global_status=1

    status=0
    isFileInList salut not in allo bonjour salut || status=1
    [ "${status}" -ne 0 ] || echo "Test 7 unsuccessful"
    [ "${status}" -eq 0 ] || global_status=1

    status=0
    isFileInList hello not in allo bonjour salut || status=1
    [ "${status}" -eq 0 ] || echo "Test 8 unsuccessful"
    [ "${status}" -eq 0 ] || global_status=1

    if [ "${global_status}" -eq 0 ]; then
        echo "All the tests are successfull!"
    else
        echo "Some of the tests do not pass!"
        exit 1
    fi
} ## End of function 'unittests'


if [ $# -eq 1 -a "${1}" = --unittests ]; then
    unittests
elif [ $# -eq 1 -a "${1}" = -u ]; then
    unittests
elif [ $# -eq 1 -a "${1}" = --help ]; then
    print_usage
elif [ $# -eq 1 -a "${1}" = -h ]; then
    print_usage
elif [ $# -lt 3 ]; then
    echo "isFileInList: need at least 3 arguments!" >&2
    echo >&2
    print_usage >&2
    exit 1
else
    isFileInList $*
fi
