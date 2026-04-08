#!/bin/bash

## Ce script sert a tester le programme 'midas.createLamTemplateGrids.Abs'.

set -eo pipefail

CREATELAMGRIDS_UNITTEST_DATABASE=${CREATELAMGRIDS_UNITTEST_DATABASE:-/home/sanl000/data/ppp7/UnitTests/midas/createLamTemplateGrids/0001}

createLAMgrids_PGM=${1:-./midas.createLamTemplateGrids.Abs}
inputFile=${2:-${CREATELAMGRIDS_UNITTEST_DATABASE}/2026033106_006}
referenceFile=${3:-${CREATELAMGRIDS_UNITTEST_DATABASE}/analysisgrid}

echo "Testing program '${createLAMgrids_PGM}'"

. r.load.dot eccc/mrd/rpn/utils/20260202
. r.load.dot rpn/code-tools/20251217/env/inteloneapi-2025.1.0
. ssmuse-sh -d main/opt/hdf5-netcdf4/parallel/intelmpi-2025.1.0/alllib/inteloneapi-2025.1.0/01

${createLAMgrids_PGM} $inputFile -grd_ext_x 280 -grd_ext_y 280

fstcomp -a ./analysisgrid -b $referenceFile -n 2>&1 | tee fstcomp.list

cat fstcomp.list | grep -vE '^[*] [*] [*] | FSTCOMP |^   \*|^1|^c_fstopl option REDUCTION32|^ Debug TG= T| SONT EGAUX$' | grep -v "^  NOM    ETIKET           IP1            IP2       IP3   E-REL-MAX   E-REL-MOY   VAR-A        C-COR        MOY-A        BIAIS       E-MAX       E-MOY" | grep -v "^  \*\*   SKIPPING RECORD \"!!\", CAN'T COMPARE  \*\*$" | grep -v "PAS DE COMPARAISON" | grep -v "No comparison possible:" | awk 'NF>0 && $(NF-1)!="0.0000E+00"' | tee fstcomp.list.grep

if [ -s fstcomp.list.grep ]; then
    echo "The 'midas.createLamTemplateGrids.Abs' test did not pass!!!"
    exit 1
else
    echo "The 'midas.createLamTemplateGrids.Abs' test is successull!!!"
    rm -f analysisgrid fstcomp.list fstcomp.list.grep
fi
