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
    if [[ "$(basename ${file})" = burpfiles_*.updated ]]; then
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
    elif [[ "$(basename ${file})" = core* ]]; then
	continue
    fi
    echo "Copying '${file}' to '${destination}'"
    ${do_copy_file} cp -rL ${file} ${destination}
done

cd ${destination}
cmcarc -f .inputs.ca -v --md5 -a *
rm -rf *
mv -v .inputs.ca inputs.ca
chmod -v -wx inputs.ca
