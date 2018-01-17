#!/bin/bash

set -ex

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
        ln -svi ${file} anal.inc.hi.model.$(basename $file)
        cmcarc -f  ${archive_path}/${archive} -v --dereference -a anal.inc.hi.model.$(basename $file)
    done
fi

for file in postalt/*_*; do
    ln -svi ${file} banco.postalt.$(basename $file)
done

cmcarc -f ${archive_path}/${archive} -v --dereference -a  anal.inc.lo.model* banco.postalt.*
cmcarc -f ${archive_path}/${archive} -t -v
