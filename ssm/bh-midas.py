#! /usr/bin/env python
# -*- coding: utf-8 -*-

from os import environ
import sys
import os
from bh import bhlib
from bh.actions import package as actions


def _init(b):
    global compiler

    environ["BH_PROJECT_NAME"]    = "midas"
    environ["CONTROL_DIR"]    = "%(BH_PACKAGE_CONTROL_DIR)s/%(BH_PACKAGE_NAMES)s/.ssm.d" % environ

def _pull(b):
    ## no need to pull MIDAS code since this script is contained directly in the MIDAS git depot
    pass

def _clean(b):
    pass

def _make(b):
    global compiler

    if b.platform == "all":
        b.shell("""
           (set -ex
            mkdir -p ${CONTROL_DIR}
            CONTROL_FILE=${CONTROL_DIR}/control.template
            echo \"Package:  x\"                                 > ${CONTROL_FILE}
            echo \"Version:  x\"                                >> ${CONTROL_FILE}
            echo \"Platform:  x\"                               >> ${CONTROL_FILE}
            echo \"Maintainer: Mark Buehner, Jean-Francois Caron, Ervig Lapalme)\"       >> ${CONTROL_FILE}
            echo \"BuildInfo: Compile script for MIDAS version ${BH_PULL_SOURCE_GIT_BRANCH}\" >> ${CONTROL_FILE}
            echo \"           with source ${BH_PULL_SOURCE}\"   >> ${CONTROL_FILE}
            echo \"Description:  Modular and Integrated Data Assimilation System (MIDAS)\">> ${CONTROL_FILE}

            cat > ${CONTROL_DIR}/control.json <<EOF
{
    \"package\": \"x\",
    \"version\": \"x\",
    \"platform\": \"x\",
    \"maintainer\": \"RPN-AD\",
    \"summary\": \"Modular and Integrated Data Assimilation System (MIDAS)\",
    \"build_info\": \"git clone -b ${BH_PULL_SOURCE_GIT_BRANCH} ${BH_PULL_SOURCE}; cd midas; cd src/programs; ./compile_all.sh\"
}
EOF
           )""")
    else:
        b.shell("""
           (set -ex
            mkdir -p ${CONTROL_DIR}
            CONTROL_FILE=${CONTROL_DIR}/control.template
            echo \"Package:  x\"                                 > ${CONTROL_FILE}
            echo \"Version:  x\"                                >> ${CONTROL_FILE}
            echo \"Platform:  x\"                               >> ${CONTROL_FILE}
            echo \"Maintainer: Mark Buehner, Jean-Francois Caron, Ervig Lapalme)\"       >> ${CONTROL_FILE}
            echo \"BuildInfo: Compiled with './compile_all.sh' version ${BH_PULL_SOURCE_GIT_BRANCH}\" >> ${CONTROL_FILE}
            echo \"           with source ${BH_PULL_SOURCE}\"   >> ${CONTROL_FILE}
            echo \"Description:  Modular and Integrated Data Assimilation System (MIDAS)\">> ${CONTROL_FILE}

            cat > ${CONTROL_DIR}/control.json <<EOF
{
    \"package\": \"x\",
    \"version\": \"x\",
    \"platform\": \"x\",
    \"maintainer\": \"RPN-AD\",
    \"summary\": \"Modular and Integrated Data Assimilation System (MIDAS)\",
    \"build_info\": \"git clone -b ${BH_PULL_SOURCE_GIT_BRANCH} ${BH_PULL_SOURCE}; cd midas; cd src/programs; ./compile_all.sh\"
}
EOF

            cd ${BH_MIDAS_TOP_LEVEL_DIR}/src/programs

            MIDAS_ABS=${BH_MIDAS_ABS} ./check_if_all_programs_compiled.sh ${ORDENV_PLAT}

            ## Compile 'obsIO' without erasing compiling directory.
            ## We will need some files to build the library.
            ./compile_program.sh obsIO full no
           )""")


def _install(b):
    if b.platform == "all":
        b.shell("""
        (set -ex
         INSTALL_DIR=${BH_INSTALL_DIR}/bin
         mkdir -p ${INSTALL_DIR}
         cp ${BH_MIDAS_TOP_LEVEL_DIR}/src/programs/commons/compile_setup.sh ${INSTALL_DIR}/midas.compile_commons.sh
         cat > ${INSTALL_DIR}/midas.compile.sh <<'EOF'
#!/bin/bash

if [ $1 = '-h' ]; then
    echo 'To compile a program with MIDAS libs, please call'
    echo '    midas.compile.sh -src ${src_files} -o ${program}'
    exit
fi

ORIGINAL_ORDENV_PLAT=${ORDENV_PLAT}

. midas.compile_commons.sh --no-rttov

if [ \"${ORDENV_PLAT}\" != \"${ORIGINAL_ORDENV_PLAT}\" ]; then
    ## load MIDAS libs
    cmpscr=$(true_path ${0})
    ssm=${cmpscr%*/*/bin/midas.compile.sh}
    echo \"loading MIDAS libs '${ssm}'\"
    . ssmuse-sh -d ${ssm}
fi

s.compile ${COMPF} -O ${FOPTMIZ} ${MPILIBDIR} -libappl midas burp_module descrip ${MPILIB} f90sqlite udfsqlite -libsys ${LIBSYS} -librmn ${LIBRMN} $*
EOF
         chmod +x ${INSTALL_DIR}/midas.compile.sh
        )""")
    else:
        b.shell("""
        (set -ex

         INSTALL_DIR=${BH_INSTALL_DIR}/bin
         mkdir -p ${INSTALL_DIR}
         for absname in ${BH_MIDAS_ABS}/midas-*_${ORDENV_PLAT}-${MIDAS_VERSION}.Abs; do
             cp ${absname} ${INSTALL_DIR}
             babsname=$(basename ${absname})
             program=$(echo ${babsname} | cut -d_ -f1)
             ln -sf ${babsname} ${INSTALL_DIR}/${program}.Abs
         done

         INCLUDE_DIR=${BH_INSTALL_DIR}/include
         LIB_DIR=${BH_INSTALL_DIR}/lib
         mkdir -p ${INCLUDE_DIR} ${LIB_DIR}
         for fmod in ${BH_MIDAS_TOP_LEVEL_DIR}/compiledir/compiledir-obsIO-${ORDENV_PLAT}_${MIDAS_VERSION}/*.mod; do
             cp ${fmod} ${INCLUDE_DIR}
         done
         [ -f ${LIB_DIR}/libmidas.a ] && rm ${LIB_DIR}/libmidas.a
         for obj in ${BH_MIDAS_TOP_LEVEL_DIR}/compiledir/compiledir-obsIO-${ORDENV_PLAT}_${MIDAS_VERSION}/*.o; do
             [ ${obj} = obsIO.o ] && continue
             if [ -f ${obj} ]; then
                 ar cru ${LIB_DIR}/libmidas.a ${obj}
             fi
         done
        )""")

if __name__ == "__main__":

    dr, b = bhlib.init(sys.argv, bhlib.PackageBuilder)
    b.actions.set("init", _init)
    b.actions.set("pull", _pull)
    b.actions.set("make", _make)
    b.actions.set("install", _install)
    b.actions.set("package", actions.package.to_ssm)
 
    b.supported_platforms = [
        "ubuntu-18.04-skylake-64",
        "sles-15-skylake-64-xc50",
        "all"
    ]
    dr.run(b)
