#! /usr/bin/env python
# -*- coding: utf-8 -*-

from os import environ
import sys
import os
from bh import bhlib
from bh.actions import package as actions


def _init(b):
    global compiler

    environ["BH_PROJECT_NAME"]    = "FindTrials"
    
    environ["CONTROL_DIR"]    = "%(BH_PACKAGE_CONTROL_DIR)s/%(BH_PACKAGE_NAMES)s/.ssm.d" % environ

def _pull(b):
    b.shell("""
           (
            set -x
            git clone --config advice.detachedHead=false --quiet --branch ${BH_PULL_SOURCE_GIT_BRANCH} --depth 1 ${BH_PULL_SOURCE} ${BH_TOP_BUILD_DIR}
            )""")

def _clean(b):
    pass

def _make(b):
    b.shell("""
           (
            mkdir -p ${CONTROL_DIR}
            CONTROL_FILE=${CONTROL_DIR}/control.template
            echo \"Package:  x\"                                > ${CONTROL_FILE}
            echo \"Version:  x\"                               >> ${CONTROL_FILE}
            echo \"Platform:  x\"                              >> ${CONTROL_FILE}
            echo \"Maintainer: RPN-AD\"                        >> ${CONTROL_FILE}
            echo \"BuildInfo: envar.findTrials from ${BH_PULL_SOURCE} for version ${BH_PULL_SOURCE_GIT_BRANCH}\" >> ${CONTROL_FILE}
            echo \"Description: script to find trials in the assimilation window\"                               >> ${CONTROL_FILE}

            cat > ${CONTROL_DIR}/control.json <<EOF
{
    \"package\": \"x\",
    \"version\": \"x\",
    \"platform\": \"x\",
    \"maintainer\": \"RPN-AD\",
    \"summary\": \"Script to find trials in the assimilation window\",
    \"build_info\": \"'envar.findTrials' from ${BH_PULL_SOURCE} for version ${BH_PULL_SOURCE_GIT_BRANCH}\"
}
EOF
           )""")

def _install(b):
    b.shell("""
        (INSTALL_DIR=${BH_INSTALL_DIR}/bin
         mkdir -p ${INSTALL_DIR}
         cd ${BH_TOP_BUILD_DIR}
         cp envar.findTrials ${INSTALL_DIR}
        )""")

if __name__ == "__main__":
    dr, b = bhlib.init(sys.argv, bhlib.PackageBuilder)
    b.actions.set("init", _init)
    b.actions.set("pull", _pull)
    b.actions.set("clean", _clean)
    b.actions.set("make", _make)
    b.actions.set("install", _install)
    b.actions.set("package", actions.package.to_ssm)
 
    b.supported_platforms = [
        "all"
    ]
    dr.run(b)
