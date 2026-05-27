#!/usr/bin/ksh
set -ex

SCRIPT=`readlink -f $0`
SCRIPT_PATH=`dirname $SCRIPT`
cd $SCRIPT_PATH

rm -rf airs.dup airs.dup2 ua.db

. ssmuse-sh -p cmda/base/master/libudfsqlite_1.6-${COMP_ARCH}

set -A files open_close query query2 insere update query_m1 query_m2 query_m3 query_m4 query_m5 query_m6 uastn
for file in ${files[@]}; do
    echo $file
    read "Enter"
    ./$file
done

