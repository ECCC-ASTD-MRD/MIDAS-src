#! /usr/bin/env python
# -*- coding: utf-8 -*-

from os import environ
import sys
import os
from bh import bhlib
from bh.actions import package as actions


def _init(b):
    global compiler

    environ["BH_PROJECT_NAME"]    = "SplitObs"
    
    # Identify a default compiler, according to the platform
    if b.platform.startswith("aix"):
        environ["AFSISIO_MAKESHIFT"] = "/home/binops/afsi/sio/env_AIX-powerpc7"
    elif b.platform.startswith("ubuntu-14.04"):
        environ["AFSISIO_MAKESHIFT"] = "/home/smco502/env_ubuntu-14.04-amd64-64"
    elif b.platform.startswith("ubuntu"):
        environ["AFSISIO_MAKESHIFT"] = "/home/binops/afsi/sio/env_Linux_x86-64/afsisio"

    environ["CONTROL_DIR"]    = "%(BH_PACKAGE_CONTROL_DIR)s/%(BH_PACKAGE_NAMES)s/.ssm.d" % environ

    # platform-specific
    if b.platform.dist == "aix":
        environ["OBJECT_MODE"] = b.platform.obj_mode

def _pull(b):
    b.shell("""
           (
            set -x
            git clone --config advice.detachedHead=false --quiet --branch ${BH_PULL_SOURCE_GIT_BRANCH} --depth 1 ${BH_PULL_SOURCE} ${BH_TOP_BUILD_DIR}
            )""")

def _clean(b):
    pass

def _make(b):
    global compiler

    b.shell("""
           (
            mkdir -p ${CONTROL_DIR}
            CONTROL_FILE=${CONTROL_DIR}/control.template
            echo \"Package:  x\"                                > ${CONTROL_FILE}
            echo \"Version:  x\"                               >> ${CONTROL_FILE}
            echo \"Platform:  x\"                              >> ${CONTROL_FILE}
            echo \"Maintainer: RPN-AD\"                        >> ${CONTROL_FILE}
            echo \"BuildInfo: Compiled with 'make splitobs_${BASE_ARCH}'\" >> ${CONTROL_FILE}
            echo \"           with source ${BH_PULL_SOURCE} for version ${BH_PULL_SOURCE_GIT_BRANCH}\" >> ${CONTROL_FILE}
            echo \"Description:  Manipulate observations files according to geographical criteria\">> ${CONTROL_FILE}

            cat > ${CONTROL_DIR}/control.json <<EOF
{
    \"package\": \"x\",
    \"version\": \"x\",
    \"platform\": \"x\",
    \"maintainer\": \"RPN-AD\",
    \"summary\": \"Manipulate observations files according to geographical criteria\",
    \"build_info\": \"git clone --branch ${BH_PULL_SOURCE_GIT_BRANCH} ${BH_PULL_SOURCE}; cd splitobs; make splitobs_${ORDENV_PLAT}\"
}
EOF
            set -ex
            cd ${BH_TOP_BUILD_DIR}
            export VERSION=${BH_PULL_SOURCE_GIT_BRANCH}
            if [ \"${BASE_ARCH}\" = AIX-powerpc7 ]; then
                s.use gmake as make
            fi
            make splitobs_%s

           )""" % b.platform)


def _install(b):
    b.shell("""
        ( set -ex
         INSTALL_DIR=${BH_INSTALL_DIR}/bin
         mkdir -p ${INSTALL_DIR}
         cd ${BH_TOP_BUILD_DIR}
         mv splitobs_%s ${INSTALL_DIR}/splitobs.Abs
         mkdir -p ${BH_INSTALL_DIR}/src
         rm -r grids/
         tar cf ${BH_INSTALL_DIR}/src/${BH_PROJECT_NAME}_${BH_PULL_SOURCE_GIT_BRANCH}.tar *
         gzip   ${BH_INSTALL_DIR}/src/${BH_PROJECT_NAME}_${BH_PULL_SOURCE_GIT_BRANCH}.tar
        )""" % b.platform)

if __name__ == "__main__":
    dr, b = bhlib.init(sys.argv, bhlib.PackageBuilder)
    b.actions.set("init", _init)
    b.actions.set("pull",_pull)
    b.actions.set("make", _make)
    b.actions.set("install", _install)
    b.actions.set("package", actions.package.to_ssm)
 
    b.supported_platforms = [
        "ubuntu-18.04-skylake-64",
        "sles-15-skylake-64-xc50"
    ]
    dr.run(b)
