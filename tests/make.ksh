#!/usr/bin/ksh
set -ex

rm -rf *.mod core *.o arjen.db toto titi  rm *.a tit

# load appropriate compilers for each architecture
if [[ -z ${COMP_ARCH} ]]; then
    if [[ "${ORDENV_PLAT}" = "aix-7.1-ppc7-64" ]]; then
        . ssmuse-sh -d hpcs/201402/01/base -d hpcs/ext/xlf_13.1.0.10
    elif [[ "${ORDENV_PLAT}" = "ubuntu-10.04-amd64-64" || "${ORDENV_PLAT}" = "ubuntu-12.04-amd64-64" ]]; then
        . ssmuse-sh -d hpcs/201402/01/base -d hpcs/201402/01/intel13sp1u2
    else
       echo "Unsupported architecture: ${ORDENV_PLAT}"
       exit 1
    fi
fi
. ssmuse-sh -p cmda/base/master/libudfsqlite_1.6-${COMP_ARCH}

set -A files open_close.f90 query.f90 query2.f90 insere.f90 update.f90 query_m1.f90 query_m2.f90 query_m3.f90 query_m4.f90 query_m5.f90 query_m6.f90 uastn.f90 
for file in ${files[@]}; do
    echo $file
    binary=`echo $file | cut -f1 -d"."` 
    s.compile -o $binary -src $file -includes ../include -libpath ../lib -optf =-static -libappl f90sqlite sqlite3 udfsqlite -debug
done

#rm -f *.export *.so *.sl *.a *.mod airs /tmp/afsdhmd/airs.dup \
#./airs.dup ./airs.dup2 ua.db query_m1 query_m2 query_m3 query_m4 query_m5 query_m6 uastn \
#query query2 open_close update insere

