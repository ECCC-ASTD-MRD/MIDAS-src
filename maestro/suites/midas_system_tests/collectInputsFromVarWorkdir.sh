#!/bin/bash

set -e

## ajouter une option pour linker les fichiers plutot que les copier

source=${1}
destination=${2}

do_copy_file=${do_copy_file}

if [ -z "${source}" ]; then
    echo "You must provide a first argument which is the source directory."
    exit 1
fi
if [ -z "${destination}" ]; then
    echo "You must provide a second argument which is the target directory where files will be copied."
    exit 1
fi

if [ -d "${destination}" ]; then
    echo "The destination directory '${destination}' already exists."
    exit 1
fi

mkdir -v ${destination}

for file in ${source}/*; do
    if [[ "$(basename ${file})" = burpfiles* ]]; then
	continue
    elif [[ "$(basename ${file})" = sqlfiles* ]]; then
	continue
    elif [[ "$(basename ${file})" = obsfiles* ]]; then
	continue
    elif [[ "$(basename ${file})" = *splitobs* ]]; then
	continue
    elif [[ "$(basename ${file})" = *reunirobs* ]]; then
	continue
    elif [[ "$(basename ${file})" = bcif_* ]]; then
	continue
    elif [[ "$(basename ${file})" = coeff_file_* ]]; then
	continue
    elif [[ "$(basename ${file})" = atms ]]; then
	continue
    elif [[ "$(basename ${file})" = split_grid ]]; then
	continue
    elif [[ "$(basename ${file})" = *_mpi_tmpdir* ]]; then
	continue
    elif [[ "$(basename ${file})" = fort.* ]]; then
	continue
    elif [[ "$(basename ${file})" = rebm* ]]; then
	continue
    elif [[ "$(basename ${file})" = rehm* ]]; then
	continue
    elif [[ "$(basename ${file})" = anlm* ]]; then
	continue
    elif [[ "$(basename ${file})" = *.Abs ]]; then
	continue
    elif [[ "$(basename ${file})" = *hostname* ]]; then
	continue
    elif [[ "$(basename ${file})" = core* ]]; then
	continue
    elif [[ "$(basename ${file})" = ensemble* ]]; then
	continue
    elif [[ "$(basename ${file})" = *nml ]]; then
	continue
    elif [[ "$(basename ${file})" = *barrier* ]]; then
	continue
    elif [[ "$(basename ${file})" =  corvert_modular.fst ]]; then
	continue
    elif [[ "$(basename ${file})" =  VAR3D_STATUS.dot ]]; then
	continue
    fi
    echo "Linking '${file}' to '${destination}'"
    ${do_copy_file} ln -s ${file} ${destination}
done

cd ${destination}
cmcarc -f .inputs.ca -v --md5 --dereference -a *
chmod -v -wx .inputs.ca
rm *

for file in ${source}/{obs,sql,burp}files*; do
    if [[ "$(basename ${file})" = *files_*.updated ]]; then
	continue
    fi
    echo "Linking '${file}' to '${destination}'"
    if [ -d "${file}" ]; then
        ${do_copy_file} ln -s ${file} ${destination}
        for thisfile in $(basename ${file})/*; do
            [ -f "${thisfile}" ] && cmcarc -f .inputs_obsfiles.ca -v --md5 --dereference -a ${thisfile}
        done
        rm $(basename ${file})
    fi
done
chmod -v -wx .inputs_obsfiles.ca

for file in ${source}/{obs,sql,burp}files*.updated; do
    if [ -d "${file}" ]; then
        echo "Linking '${file}' to '${destination}'"
        ${do_copy_file} ln -s ${file} ${destination}
        for thisfile in $(basename ${file})/*; do
            [ -f "${thisfile}" ] && cmcarc -f .results_obsfiles.ca -v --md5 --dereference -a ${thisfile}
        done
        rm $(basename ${file})
    fi
done
if [ -f .results_obsfiles.ca ]; then
    chmod -v -wx .results_obsfiles.ca
fi

for bah in rebm rehm anlm; do
    copied_something=no
    for file in ${source}/${bah}* ${source}/../${bah}*; do
        echo "Linking '${file}' to '${destination}'"
        if [ -f "${file}" ]; then
            ${do_copy_file} ln -s ${file} ${destination}
            copied_something=yes
        fi
    done
    if [ "${copied_something}" = yes ]; then
        cmcarc -f .results_${bah}.ca -v --md5 --dereference -a ${bah}*
        chmod -v -wx .results_${bah}.ca
        rm *
    fi
done

mv -v .inputs.ca inputs.ca
mv -v .inputs_obsfiles.ca inputs_obsfiles.ca
if [ -f .results_obsfiles.ca ]; then
    mv -v .results_obsfiles.ca results_obsfiles.ca
fi
for bah in rebm rehm anlm; do
    if [ -f .results_${bah}.ca ]; then
        mv -v .results_${bah}.ca results_${bah}.ca
    fi
done

if [ -d ${source}/ensemble ]; then
    echo "Linking 'ensemble' to '${destination}'"
    ${do_copy_file} ln -s ${source}/ensemble ensemble
    cmcarc -f inputs_ensemble.ca -v --md5 --dereference -a ensemble/*
    chmod -v -wx inputs_ensemble.ca
    rm ensemble
fi
