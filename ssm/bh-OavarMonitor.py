#! /usr/bin/env python
# -*- coding: utf-8 -*-

from os import environ
import sys
import os
from bh import bhlib
from bh.actions import package as actions


def _init(b):
    global compiler

    environ["BH_PROJECT_NAME"]    = "OavarMonitor"
    
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
            echo \"BuildInfo: Compiled with 'gcc oavar.monitor.c -o oavar.monitor.Abs'\" >> ${CONTROL_FILE}
            echo \"           with source ${BH_PULL_SOURCE} for version ${BH_PULL_SOURCE_GIT_BRANCH}\" >> ${CONTROL_FILE}
            echo \"Description:  Monitor the file 'VAR3D_STATUS.dot' to detect when 'rebm' files are ready to be processed.\">> ${CONTROL_FILE}

            cat > ${CONTROL_DIR}/control.json <<EOF
{
    \"package\": \"x\",
    \"version\": \"x\",
    \"platform\": \"x\",
    \"maintainer\": \"RPN-AD\",
    \"summary\": \"Monitor the file 'VAR3D_STATUS.dot' to detect when 'rebm' files are ready to be processed\",
    \"build_info\": \"Compiled with 'gcc oavar.monitor.c -o oavar.monitor.Abs' with source ${BH_PULL_SOURCE} for version ${BH_PULL_SOURCE_GIT_BRANCH}\"
}
EOF

            cd ${BH_TOP_BUILD_DIR}
            gcc oavar.monitor.c -o oavar.monitor.Abs
           )""")


def _install(b):
    b.shell("""
        (
         INSTALL_DIR=${BH_INSTALL_DIR}/bin
         mkdir -p ${INSTALL_DIR}
         cd ${BH_TOP_BUILD_DIR}
         mv oavar.monitor.Abs ${INSTALL_DIR}/oavar.monitor.Abs
         mkdir -p ${BH_INSTALL_DIR}/src
         tar cf ${BH_INSTALL_DIR}/src/${BH_PROJECT_NAME}_${BH_PULL_SOURCE_GIT_BRANCH}.tar *
         gzip ${BH_INSTALL_DIR}/src/${BH_PROJECT_NAME}_${BH_PULL_SOURCE_GIT_BRANCH}.tar
        )""")

if __name__ == "__main__":
    dr, b = bhlib.init(sys.argv, bhlib.PackageBuilder)
    b.actions.set("init", _init)
    b.actions.set("pull",_pull)
    b.actions.set("make", _make)
    b.actions.set("install", _install)
    b.actions.set("package", actions.package.to_ssm)
 
    b.supported_platforms = [
        "aix-7.1-ppc7-64",
        "ubuntu-10.04-amd64-64",
        "ubuntu-12.04-amd64-64",
        "ubuntu-14.04-amd64-64",
        "sles-11-broadwell-64-xc40"
    ]
    dr.run(b)
