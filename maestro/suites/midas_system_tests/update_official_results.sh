#!/bin/bash

set -ex

## Execute this script in the 'check/input' directory
## or specify this directory as the second argument

archive_path=$1

if [ $# -eq 2 ]; then
    cd $2
fi

archive=$(/bin/ls -1 *.ca)

echo ${archive_path}/${archive}

if [ -d "${archive_path}" ]; then
    echo "archive_path=${archive_path} exists!"
    exit 1
else
    mkdir -vp ${archive_path}
fi

cmcarc -f ${archive_path}/${archive} -v --dereference -a anal.anl.model.*

if [ -d anal.inc.hi.model.* ]; then
    for file in anal.inc.hi.model.*/*; do
        bfile=$(basename ${file})
        ext=$(echo ${bfile} | cut -d_ -f2)
	if [ "${ext}" = inc ]; then
            newfile=anal.inc.hi.model.$(echo ${bfile} | cut -d_ -f1)_000_$(echo ${bfile} | cut -d_ -f2-)
        else
            newfile=anal.inc.hi.model.${bfile}
        fi
        ln -svi ${file} ${newfile}
        cmcarc -f  ${archive_path}/${archive} -v --dereference -a ${newfile}
    done
fi

for file in postalt/*_*; do
    ln -svi ${file} banco.postalt.$(basename $file)
done

cmcarc -f ${archive_path}/${archive} -v --dereference -a  anal.inc.lo.model* banco.postalt.*
cmcarc -f ${archive_path}/${archive} -t -v

input=$(readlink ${archive})
input_source=$(basename $(dirname $(dirname ${input})))
ln -sv ../${input_source}/inputs $(dirname ${archive_path})
