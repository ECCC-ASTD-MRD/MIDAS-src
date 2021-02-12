#!/bin/bash

# get list of functions
functionList=" "
for fileName in */*.f90
do
    functionNames=`grep -v '^ *\!' $fileName | grep '^ *function ' | awk '{print $2}' | cut -f1 -d'('`
    for functionName in `echo $functionNames`
    do
	functionList="$functionList $functionName"
    done

    functionNames=`grep -v '^ *\!' $fileName | grep '^ *real.* *function ' | awk '{print $3}' | cut -f1 -d'('`
    for functionName in `echo $functionNames`
    do
	functionList="$functionList $functionName"
    done

    functionNames=`grep -v '^ *\!' $fileName | grep '^ *integer.* *function ' | awk '{print $3}' | cut -f1 -d'('`
    for functionName in `echo $functionNames`
    do
	functionList="$functionList $functionName"
    done

    functionNames=`grep -v '^ *\!' $fileName | grep '^ *logical.* *function ' | awk '{print $3}' | cut -f1 -d'('`
    for functionName in `echo $functionNames`
    do
	functionList="$functionList $functionName"
    done

    functionNames=`grep -v '^ *\!' $fileName | grep '^ *character.* *function ' | awk '{print $3}' | cut -f1 -d'('`
    for functionName in `echo $functionNames`
    do
	functionList="$functionList $functionName"
    done
done

echo "FUNCTION_LIST="
echo "$functionList"
echo ""
echo ""

for functionName in `echo $functionList`
do
    # one line write statements
    grep -ir '^ *write' |grep -i "$functionName *("
    # also try to detect write statements with line continuation
    grep -irA 2 '^ *write.*& *$' |grep -i "$functionName *(" |grep -iv 'write'
done

echo "DONE"


