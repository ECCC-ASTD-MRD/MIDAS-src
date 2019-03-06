#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# bh-AddAnalInc.py
 
from os import environ
import sys
from bh import bhlib, actions
import os

 
def _init(b):
    global compiler
    
    environ["BH_PROJECT_NAME"]                    = "AddAnalInc"
    environ["LIBRMN_REVISION"]                    = "rmn"
    
    environ["CONTROL_DIR"]    = "%(BH_PACKAGE_CONTROL_DIR)s/%(BH_PACKAGE_NAMES)s/.ssm.d" % environ

    # platform-specific
    if b.platform.dist == "aix":
        environ["OBJECT_MODE"] = b.platform.obj_mode

    if b.platform == "aix-7.1-ppc7-64":
        compiler="xlf13"
    elif b.platform == "ubuntu-14.04-amd64-64":
        compiler="intel-2016.1.156"
    elif b.platform in ["ubuntu-10.04-amd64-64","ubuntu-12.04-amd64-64"]:
        compiler="intel13sp1u2"
    elif b.platform in "sles-11-broadwell-64-xc40":
        compiler="PrgEnv-intel-5.2.82"
    else:
        print "ERROR: The platform '%s' is not supported for this application" % b.platform
        sys.exit(1)


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

    ## This string is only used when 'AAI_COMPILATION_METHOD=Compiling_AAI' which is the old method of compiling
    if compiler == "xlf13":
        ssmuse_st = """
            . ssmuse-sh -d hpcs/201402/00/base -d hpcs/ext/xlf_13.1.0.10
            s.use gmake as make
            s.use /home/ordenv/ssm-domains/ssm-setup-1.0/dot-profile-setup_1.0_multi/hostbin/AIX/gnu_find as find
            s.use /home/ordenv/ssm-domains/ssm-setup-1.0/dot-profile-setup_1.0_multi/hostbin/AIX/gnu_find as gnu_find
            """

    elif compiler == "pgi9xx":
        ssmuse_st = """
            . ssmuse-sh -d hpcs/13b/03/base
            . s.ssmuse.dot pgi9xx devtools
            """

    elif compiler == "pgi1301":
        ssmuse_st = """
           . ssmuse-sh -d hpcs/13b/04/base -d hpcs/13b/04/pgi1301
           """

    elif compiler == "pgi1401":
        ssmuse_st = """
            . ssmuse-sh -d /ssm/net/hpcs/201402/01/base -d /ssm/net/hpcs/201402/01/pgi1401
            """

    elif compiler == "intel13sp1u2":
        ssmuse_st = """
            . ssmuse-sh -d /ssm/net/hpcs/201402/01/base -d /ssm/net/hpcs/201402/01/intel13sp1u2
        """

    elif compiler == "intel-2016.1.156":
        ssmuse_st="This code is not necessary!"
    elif compiler == "PrgEnv-intel-5.2.82":
        ssmuse_st="This code is not necessary!"

    b.shell("""
           (set -e

           if [ \"${AAI_COMPILATION_METHOD}\" = Compiling_AAI ]; then
              export SOURCE_CODE=${BH_TOP_BUILD_DIR}/src
              export COMPILE_DIRECTORY=${BH_TOP_BUILD_DIR}/

              %s

              . s.ssmuse.dot /ssm/net/rpn/libs/15.1
              . s.ssmuse.dot SI/ezinterpv_15.0.3

              ${BH_TOP_BUILD_DIR}/src/Compiling_AAI
           else
            cd ${BH_TOP_BUILD_DIR}
            
            . ./aai.setup_compilation.dot

            cd src
            make AddAnalInc \"LIBRMN_REVISION=${AAI_LIBRMN}\" \"LIBAPPL=${AAI_LIBAPPL}\" OPTF=\"${OPTF}\" OPTMIZ=${OPTMIZ}

            # Add the revision number to the file name
            cd ${BH_TOP_BUILD_DIR}/Build/${EC_ARCH}
            absolute=`ls *.Abs`
            mv $absolute ${absolute%%.Abs}_${BH_PULL_SOURCE_GIT_BRANCH}.Abs
           fi

            mkdir -p ${CONTROL_DIR}
            CONTROL_FILE=${CONTROL_DIR}/control.template
            echo \"Package:  x\"                > ${CONTROL_FILE}
            echo \"Version:  x\"               >> ${CONTROL_FILE}
            echo \"Platform:  x\"              >> ${CONTROL_FILE}
            echo \"Maintainer: Jeff Blezius\"                               >> ${CONTROL_FILE}
            echo \"BuildInfo: Compiled from GIT revision '${BH_PULL_SOURCE_GIT_BRANCH}'\"  >> ${CONTROL_FILE}
            echo \"Description:  Add Analyis Increment utility\"            >> ${CONTROL_FILE}

           )""" % (ssmuse_st))


def _install(b):
    global compiler
    
    if compiler == "xlf13":
        ssmuse_st = """
            . ssmuse-sh -d hpcs/ext/xlf_13.1.0.10"""

    elif compiler == "pgi9xx":
        ssmuse_st = """
            . s.ssmuse.dot pgi9xx"""

    elif compiler == "pgi1301":
        ssmuse_st = """
            . ssmuse-sh -d hpcs/13b/04/pgi1301"""

    elif compiler == "pgi1401":
        ssmuse_st = """
            . ssmuse-sh -d /ssm/net/hpcs/201402/01/base -d /ssm/net/hpcs/201402/01/pgi1401
            """

    elif compiler == "intel13sp1u2":
        ssmuse_st = """
            . ssmuse-sh -d /ssm/net/hpcs/201402/01/base -d /ssm/net/hpcs/201402/01/intel13sp1u2
            """

    elif compiler == "intel13sp1u2":
        ssmuse_st = """
            . ssmuse-sh -d /ssm/net/hpcs/201402/01/base -d /ssm/net/hpcs/201402/01/intel13sp1u2
            """

    elif compiler == "intel-2016.1.156":
        ssmuse_st = """
            . ${BH_BUILD_DIR}/aai.setup_compilation.dot
        """
    elif compiler == "PrgEnv-intel-5.2.82":
        ssmuse_st = """
            . ${BH_BUILD_DIR}/aai.setup_compilation.dot
        """

    b.shell("""
        (set -e
         %s

         INSTALL_DIR_BIN=${BH_INSTALL_DIR}/bin
         mkdir -p ${INSTALL_DIR_BIN}
         mv ${BH_BUILD_DIR}/Build/${EC_ARCH}/*.Abs ${INSTALL_DIR_BIN}/

         # Link without the revision number
         cd ${INSTALL_DIR_BIN}/
         absolute=`ls *.Abs`
         ln -s ${absolute} AddAnalInc.${absolute##*.}
        )""" % ssmuse_st)


if __name__ == "__main__":
    dr, b = bhlib.init(sys.argv, bhlib.PackageBuilder)
    b.actions.set("init", _init)
    b.actions.set("pull", _pull)
    b.actions.set("clean", _clean)
    b.actions.set("make", _make)
    b.actions.set("install", _install)
    b.actions.set("package", actions.package.to_ssm)
 
    b.supported_platforms = [
        "ubuntu-10.04-amd64-64",
        "ubuntu-12.04-amd64-64",
        "ubuntu-14.04-amd64-64",
        "aix-7.1-ppc7-64",
        "sles-11-broadwell-64-xc40"
    ]
    dr.run(b)
