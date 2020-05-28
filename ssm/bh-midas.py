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

           )""")


def _install(b):
    if b.platform ! "all":
        b.shell("""
        (set -ex

         INSTALL_DIR=${BH_INSTALL_DIR}/bin
         mkdir -p ${INSTALL_DIR}
         for prog in ${BH_MIDAS_TOP_LEVEL_DIR}/src/programs/*.f90; do
             progname=$(basename ${prog} .f90)
             absname=${BH_MIDAS_ABS}/midas-${progname}_${ORDENV_PLAT}-${MIDAS_VERSION}.Abs
             cp ${absname} ${INSTALL_DIR}
             babsname=$(basename ${absname})
             program=$(echo ${babsname} | cut -d_ -f1)
             ln -sf ${babsname} ${INSTALL_DIR}/${program}.Abs
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
