#! /usr/bin/env python
# -*- coding: utf-8 -*-

from os import environ
import sys
from bh import bhlib, actions
import os


def _init(b):
    environ["BH_PROJECT_NAME"]    = "VarScripts"

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
            echo \"Maintainer: arma (E. Lapalme, J. Blezius)\" >> ${CONTROL_FILE}
            echo \"BuildInfo: Copied from \${ARMNLIB}/modeles/ANAL_shared/scripts/env/r.tripotenml,\" >> ${CONTROL_FILE}
            echo \"           and ${BH_PULL_SOURCE} for version ${BH_PULL_SOURCE_GIT_BRANCH}\"        >> ${CONTROL_FILE}
            echo \"Description:  Scripts for Variational Assimilation Code\"         >> ${CONTROL_FILE}
           )""")

def _install(b):
    b.shell("""
        (INSTALL_DIR=${BH_INSTALL_DIR}/bin
         mkdir -p ${INSTALL_DIR}
         cd ${BH_TOP_BUILD_DIR}
         cp oavar.* r.tripotenml ${INSTALL_DIR}
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