#! /bin/sh

if [ "${ORDENV_PLAT}" != rhel-8-icelake-64 -a "${ORDENV_PLAT}" != rhel-9-graniterapids-64 ];then
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported."
    return 1
fi

##
## DO NOT MODIFY THIS FILE
##
## as it is part of the versioned repository and contains the default configuration
## for the build environment.
##
## If you need to change the default configuration provided by this file, either
## export the desired variables in your environment before sourcing it or define
## them in your profile.
##

__run_cmake=true
__show_instructions=true
__cd_to_build_directory=true
__fresh_build_directory=false

function __default_compiler {
    if [ "${ORDENV_PLAT}" = rhel-8-icelake-64 ]; then
        echo inteloneapi-2022.1.2
    elif [ "${ORDENV_PLAT}" = rhel-9-graniterapids-64 ]; then
        echo inteloneapi-2025.1.0
    else
        echo "config.dot.sh: This platform '${ORDENV_PLAT}' is not supported.  Only 'rhel-8-icelake-64' and 'rhel-9-graniterapids-64' are    ." >&2
        exit 1
    fi
}
__compiler=$(__default_compiler)

while [[ $# > 0 ]]; do
    arg=${1}

    if [[ "${arg}" = --build ]]; then
        if [[ $# -ge 2 ]]; then
            __midas_compile_dir_main=${2}
            shift
        else
            echo "config.dot.sh: this option '${arg}' needs a build directory" >&2
        fi
    elif [[ "${arg}" = --build-id ]]; then
        if [[ $# -ge 2 ]]; then
            __midas_compile_build_id=${2}
            shift
        else
            echo "config.dot.sh: this option '${arg}' needs a build id" >&2
        fi
    elif [[ "${arg}" = --compiler ]]; then
        if [[ $# -lt 2 ]]; then
            echo "config.dot.sh: The option '--compiler' must be followed by a compiler identificator" >&2
            exit 1
        elif [[ ${2} = -* || -z "${2}" ]]; then
            echo "config.dot.sh: The option '--compiler' must be followed by a compiler identificator" >&2
            exit 1
        fi
        __compiler=${2}
        shift
    elif [[ "${arg}" = --no-cmake ]]; then
        __run_cmake=false
    elif [[ "${arg}" = --cmake ]]; then
        __run_cmake=true
    elif [[ "${arg}" = --show-instructions ]]; then
        __show_instructions=true
    elif [[ "${arg}" = --cd-build ]]; then
        __cd_to_build_directory=true
    elif [[ "${arg}" = --no-show-instructions ]]; then
        __show_instructions=false
    elif [[ "${arg}" = --no-cd-build ]]; then
        __cd_to_build_directory=false
    elif [[ "${arg}" = --fresh ]]; then
        __fresh_build_directory=true
    elif [[ "${arg}" = --no-fresh ]]; then
        __fresh_build_directory=false
    elif [[ "${arg}" = -h || "${arg}" = -help || "${arg}" = --help ]]; then
        echo "config.dot.sh: "
        echo "        --build: explicitly specify a build directory"
        echo "        --build-id: specify an extension the build directory name"
        echo "        --compiler: set the compiler to use (default is ${__compiler})"
        echo "        --no-cmake: avoid running cmake to create the build directory and leave it to the user"
        echo "        --cmake: do run cmake to prepare the build directory (default)"
        echo "        --no-show-instructions: do not print any instructions for the user"
        echo "        --no-cd-build: do not move to build directory "
        echo "        --show-instructions: do print any instructions for the user (default)"
        echo "        --cd-build: do move to build directory (default)"
        echo "        --fresh: do clean to build directory to start fresh"
        echo "        --no-fresh: do not clean to build directory to use previous builds (default)"
        echo "        -h|-help|--help: show this help"
        __run_cmake=stop
        break
    elif [[ "${arg}" = -* ]]; then
        echo "config.dot.sh: this option '${arg}' is not allowed (see '-h' for allowed options)" >&2
        __run_cmake=stop
        break
    fi

    shift
done ## End of 'while [[ $# > 0 ]]'

if [ "${MIDAS_CONFIG_ALREADY_SOURCED}" == "true" ]; then
    echo "+-----------------------------------------+"
    echo "|                                         |"
    echo "| WARNING: config.dot.sh already sourced! |"
    echo "|     possible environment conflict       |"
    echo "|   open a new shell before proceeding    |"
    echo "|                                         |"
    echo "+-----------------------------------------+"
    return 1
fi

if [ "${__show_instructions}" != true -a "${__show_instructions}" != false ]; then
    echo "config.dot.sh: The variable '__show_instructions' can only be 'true' or 'false' and not '${__show_instructions}'." >&2
    __run_cmake=stop
fi

if [ "${__cd_to_build_directory}" != true -a "${__cd_to_build_directory}" != false ]; then
    echo "config.dot.sh: The variable '__cd_to_build_directory' can only be 'true' or 'false' and not '${__cd_to_build_directory}'." >&2
    __run_cmake=stop
fi

if [ "${__fresh_build_directory}" != true -a "${__fresh_build_directory}" != false ]; then
    echo "config.dot.sh: The variable '__fresh_build_directory' can only be 'true' or 'false' and not '${__fresh_build_directory}'." >&2
    __run_cmake=stop
fi

__supported_compilers="'inteloneapi-2022.1.2', 'inteloneapi-2025.1.0', 'gnu-14.2.0', 'gnu-15.2.0'"
if [ "${__compiler}" = list ]; then
    echo "The compiler supported are 'inteloneapi-2022.1.2', 'inteloneapi-2025.1.0', 'gnu-14.2.0' and 'gnu-15.2.0'."
    echo "The supported compilers are ${__supported_compilers}."
    __run_cmake=stop
elif [ "${__compiler}" != inteloneapi-2022.1.2 -a "${__compiler}" != gnu-14.2.0 -a "${__compiler}" != gnu-15.2.0 -a "${__compiler}" != inteloneapi-2025.1.0 ]; then
    echo "config.dot.sh: The compiler '${__compiler}' is not supported.  Only ${__supported_compilers} are." >&2
    __status=false
    __run_cmake=stop
fi

if [ "${__run_cmake}" != stop -a "${__run_cmake}" != true -a "${__run_cmake}" != false ]; then
    echo "config.dot.sh: The variable '__run_cmake' can only be 'stop', 'true' or 'false' and not '${__run_cmake}'." >&2
    __run_cmake=stop
fi

cmake_options() {
    typeset __cmake_options_options__=

    if [[ "${MIDAS_COMPILE_ADD_DEBUG_OPTIONS}" = yes || -n "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}" ]]; then
        __cmake_options_options__="-DCMAKE_BUILD_TYPE=Debug -DWITH_WARNINGS=TRUE -DEXTRA_CHECKS=TRUE"
    fi

    if [[ "${MIDAS_COMPILE_VERBOSE}" != FALSE ]]; then
        typeset __cmake_options_verbose__=-DCMAKE_VERBOSE_MAKEFILE=${MIDAS_COMPILE_VERBOSE}
        if [[ -n "${__cmake_options_options__}" ]]; then
            __cmake_options_options__="${__cmake_options_options__} ${__cmake_options_verbose__}"
        else
            __cmake_options_options__=${__cmake_options_verbose__}
        fi
    fi

    echo "${__cmake_options_options__}"
} ## End of 'debug_options'

if [ "${__run_cmake}" != stop ]; then

    __toplevel=$(git rev-parse --show-toplevel)
    __revstring=$(${__toplevel}/midas.version)
    __revnum=$(echo ${__revstring} | sed -e 's/v_\([^-]*\)-.*/\1/')
    __status=true

    ###########################################################
    ##
    ##  USER CONFIGURATION
    ##
    ###########################################################

    ## The variable 'MIDAS_COMPILE_DIR_MAIN' can be defined as
    ## 'build_directory_local_to_the_repository' and so the build
    ## directory will be in '${__toplevel}/compiledir' (see variable '${__compiledir_link}' below).
    ## If it is not defined, we will use a default which appends the
    ## basename of the toplevel directory to '${HOME}/data_maestro/ords/midas-bld'.

    MIDAS_COMPILE_ADD_DEBUG_OPTIONS=${MIDAS_COMPILE_ADD_DEBUG_OPTIONS:-no}
    MIDAS_COMPILE_APPEND_VERSION_ID_BUILDDIR=${MIDAS_COMPILE_APPEND_VERSION_ID_BUILDDIR:-true}
    MIDAS_COMPILE_CODECOVERAGE_DATAPATH=${MIDAS_COMPILE_CODECOVERAGE_DATAPATH:-}
    if [ -z "${MIDAS_COMPILE_FRONTEND}" ]; then
        if [ "${ORDENV_PLAT}" = rhel-8-icelake-64 ]; then
            MIDAS_COMPILE_FRONTEND=ppp5
        elif [ "${ORDENV_PLAT}" = rhel-9-graniterapids-64 ]; then
            MIDAS_COMPILE_FRONTEND=ppp7
        fi
    fi
    MIDAS_COMPILE_CLEAN=${MIDAS_COMPILE_CLEAN:-true}
    MIDAS_COMPILE_COMPF_GLOBAL=${MIDAS_COMPILE_COMPF_GLOBAL:-}
    MIDAS_COMPILE_HEADNODE_FRONTEND=${MIDAS_COMPILE_HEADNODE_FRONTEND:-false}
    MIDAS_COMPILE_JOBNAME=${MIDAS_COMPILE_JOBNAME:-midasCompilation}
    MIDAS_COMPILE_KEEP_LISTING=${MIDAS_COMPILE_KEEP_LISTING:-false}
    MIDAS_COMPILE_NCORES=${MIDAS_COMPILE_NCORES:-8}
    MIDAS_COMPILE_VERBOSE=${MIDAS_COMPILE_VERBOSE:-FALSE}
    MIDAS_COMPILE_OPTIMIZE_REPORT=${MIDAS_COMPILE_OPTIMIZE_REPORT:-no}

    __compiledir_link=${__compiledir_link:-${__toplevel}/compiledir}

    export MIDAS_SOURCE_DIR=${__toplevel}
    if [[ -z "${__midas_compile_dir_main}" ]]; then
        ## If 'MIDAS_COMPILE_DIR_MAIN' is equal to the special value
        ## 'build_directory_local_to_the_repository', then we keep the build
        ## directory local to the Git repository.
        if [ "${MIDAS_COMPILE_DIR_MAIN}" = build_directory_local_to_the_repository ]; then
            echo "Creating '${__compiledir_link}' since 'MIDAS_COMPILE_DIR_MAIN' is equal to this special value"
            [ ! -d "${__compiledir_link}" ] && mkdir ${__compiledir_link}
            __midas_compile_dir_main=${__compiledir_link}
        else
            if [ -n "${MIDAS_COMPILE_DIR_MAIN}" ]; then
                if [[ "${MIDAS_COMPILE_DIR_MAIN}" != /* ]]; then
                    echo "Please provide of value of MIDAS_COMPILE_DIR_MAIN which is an absolute path"
                    echo "   MIDAS_COMPILE_DIR_MAIN=${MIDAS_COMPILE_DIR_MAIN}"
                    echo "was given"
                    return 1
                fi
                __midas_compile_dir_main=${MIDAS_COMPILE_DIR_MAIN}
            else
                ## If the variable 'MIDAS_COMPILE_DIR_MAIN' is not defined,
                ## we add the leaf part of the toplevel directory to
                ## '${HOME}/data_maestro/ords/midas-bld'.
                __toplevel_leaf=$(basename ${__toplevel})
                __midas_compile_dir_main=${HOME}/data_maestro/ords/midas-bld/${__toplevel_leaf}
            fi

            if [ ! -d "${__midas_compile_dir_main}" ]; then
                mkdir -p ${__midas_compile_dir_main}
            fi

            ##  linking the build directory where it used to be
            if [ -d "${__compiledir_link}" -o -L "${__compiledir_link}" ]; then
                echo "${__compiledir_link} already exists: not creating link."
            else
                ln -s ${__midas_compile_dir_main} ${__compiledir_link}
            fi
        fi ## End of 'else' associated to 'if [ "${MIDAS_COMPILE_DIR_MAIN}" = build_directory_local_to_the_repository ]'
    fi ## End of 'if [[ -z "${__midas_compile_dir_main}" ]]'

    export MIDAS_COMPILE_DIR_MAIN=${__midas_compile_dir_main}

    ###########################################################
    ##  LESS-USER-FRIENDLY CONFIGURATION
    ##
    ##  these should not be changed unless you know what you're doing
    ##  it can impact the maestro testing suite or the cleaning targets
    ##  in unwated ways
    MIDAS_ABS_LEAFDIR=${MIDAS_ABS_LEAFDIR:-midas_abs}
    __install_always_midas=true

    if [ -n "${__midas_compile_build_id}" ]; then
        __versionid_build=-${__midas_compile_build_id}
    elif [ "${MIDAS_COMPILE_APPEND_VERSION_ID_BUILDDIR}" = true ]; then
        __versionid_build=-${__revstring}
    else
        __versionid_build=
    fi
    export MIDAS_COMPILE_DIR_BUILD=${MIDAS_COMPILE_DIR_MAIN}/midas_bld${__versionid_build}/${__compiler}

    ###########################################################
    ##  compilation needed for compilation
    ##
    ## -- should not change that

    # User-specified compilation options
    #export MIDAS_COMPILE_COMPF_GLOBAL="-DCODEPRECISION_INCR_REAL_SINGLE"
    #export MIDAS_COMPILE_COMPF_GLOBAL="-DCODEPRECISION_SPECTRANS_REAL_SINGLE"
    if [ -n "${MIDAS_COMPILE_COMPF_GLOBAL}" ];then
        echo "..."
        echo "... Additional user-specified compilation options = ${MIDAS_COMPILE_COMPF_GLOBAL}"
        echo "..."
        export MIDAS_COMPILE_COMPF_GLOBAL
    fi

    if [ "${__compiler}" = gnu-14.2.0 ]; then
        __original_ordenv_plat__=${ORDENV_PLAT}

        echo "... loading gnu 14.2.0 compiler"
        export MODULEPATH=/home/sidr000/modules:${MODULEPATH}
        . r.load.dot rpn/code-tools/20260219/env/exp/gnu-14.2.0
    elif [ "${__compiler}" = gnu-15.2.0 ]; then
        __original_ordenv_plat__=${ORDENV_PLAT}

        echo "... loading gnu 15.2.0 compiler"
        export MODULEPATH=/home/sidr000/modules:${MODULEPATH}
        . r.load.dot rpn/code-tools/20260219/env/exp/gnu-15.2.0
    else
        echo "... loading rpn/code-tools/20251217/env/${__compiler}"
        . r.load.dot rpn/code-tools/20260219/env/${__compiler}
    fi

    echo "... loading eccc/mrd/rpn/libs/20260401-beta"
    . r.load.dot eccc/mrd/rpn/libs/20260401-beta

    echo "... loading eccc/mrd/rpn/utils/20260401-beta/burp-tools_20.0.16-${COMP_ARCH}_${ORDENV_PLAT}"
    . r.load.dot eccc/mrd/rpn/utils/20260401-beta/burp-tools_20.0.16-${COMP_ARCH}_${ORDENV_PLAT}

    if [ "${__compiler}" = gnu-14.2.0 ]; then
        ## This is needed for CMake to find the library
        export f90sqlite_INSTALLDIR=/home/erv000/SSM/midas-external-libs/1.0-rc2/${COMP_ARCH}/${ORDENV_PLAT}

        ## Loading the compiler affects 'ORDENV_PLAT' and for the
        ## rest, we need to use the original value.
        export ORDENV_PLAT=${__original_ordenv_plat__}
        unset __original_ordenv_plat__

        echo "... loading perftools, udfsqlite and f90sqlite for ${COMP_ARCH}"
        . r.load.dot /home/erv000/SSM/midas-external-libs/1.0-rc2/${COMP_ARCH}
    elif [ "${__compiler}" = gnu-15.2.0 ]; then
        ## This is needed for CMake to find the library
        echo "... providing path for f90sqlite, hdf5 and netcdf to CMake"
        export f90sqlite_INSTALLDIR=/home/mad001/ssm/f90sqlite_1.5-gnu-15.2.0_rhel-9-amd64-64
        export hdf5_INSTALLDIR=/home/mad001/ssm/hdf5_1.14.6-gnu-15.2.0_rhel-9-amd64-64
        export netcdf_INSTALLDIR=/home/mad001/ssm/netcdf-fortran_4.6.2-gnu-15.2.0_rhel-9-amd64-64

        ## Loading the compiler affects 'ORDENV_PLAT' and for the
        ## rest, we need to use the original value.
        export ORDENV_PLAT=${__original_ordenv_plat__}
        unset __original_ordenv_plat__

        echo "... loading perftools, and udfsqlite for ${COMP_ARCH}"
        . ssmuse-sh -x /fs/homeu3/eccc/mrd/rpnad/mad001/ssm/perftools_2.1-gnu-15.2.0_rhel-9-amd64-64
        . ssmuse-sh -x /fs/homeu3/eccc/mrd/rpnad/mad001/ssm/udfsqlite_1.19-gnu-15.2.0_rhel-9-amd64-64
    else
        echo "... loading eccc/cmd/cmda/libs/20260401-beta/${COMP_ARCH}"
        . ssmuse-sh -d eccc/cmd/cmda/libs/20260401-beta/${COMP_ARCH}

        echo "... loading hdf5-netcdf4"
        if [ "${ORDENV_PLAT}" = rhel-8-icelake-64 ]; then
            . ssmuse-sh -d main/opt/hdf5-netcdf4/serial/static/${COMP_ARCH}/01
        elif [ "${ORDENV_PLAT}" = rhel-9-graniterapids-64 ]; then
            . ssmuse-sh -d main/opt/hdf5-netcdf4/parallel/intelmpi-2025.1.0/alllib/${COMP_ARCH}/01
        fi

        if [ "${ORDENV_PLAT}" = rhel-8-icelake-64 ]; then
            perftools_version=2.0
        elif [ "${ORDENV_PLAT}" = rhel-9-graniterapids-64 ]; then
            perftools_version=2.1
        fi

        echo "... loading main/opt/perftools/perftools-${perftools_version}/${COMP_ARCH}"
        . ssmuse-sh -x main/opt/perftools/perftools-${perftools_version}/${COMP_ARCH}
    fi

    if [ "${MIDAS_COMPILE_ADD_DEBUG_OPTIONS:-no}" = yes ]; then
        __rttovdebug__=-debug
    else
        __rttovdebug__=
    fi
    if [ "${__compiler}" = gnu-15.2.0 ]; then
        __rttov_path__=/fs/site7/eccc/mrd/rpnad/erv000/midas/rttov_v13-17/tmp/rttov13/2.1.0-17-g93b9e11
        echo "... loading  ${__rttov_path__}${__rttovdebug__}"
        . ssmuse-sh -x ${__rttov_path__}${__rttovdebug__}
        export rttov_INSTALLDIR=${__rttov_path__}${__rttovdebug__}
    else
        export RTTOV_VERSION=2.1.0 ## This variable is used in '../CMakeLists.txt' for the script 'midas-config'
        __rttov_path__=/fs/ssm/eccc/mrd/rpn/anl/rttov13/${RTTOV_VERSION}/${COMP_ARCH}${__rttovdebug__}
        echo "... loading ${__rttov_path__}"
        . r.load.dot ${__rttov_path__}
    fi

    unset __rttov_path__ __rttovdebug__

    # add compiler option to produce reports on code optimization and deactivate cleaning
    if [ "${MIDAS_COMPILE_OPTIMIZE_REPORT:-no}" = yes ]; then
        echo "... > !WARNING! Compiler optimization reports will be produced in the compile directory."
        echo "... >           To be able to see them, we ensure cleaning is not activated."
        MIDAS_COMPILE_CLEAN=false
    fi

    if [ "${MIDAS_COMPILE_ADD_DEBUG_OPTIONS:-no}" = yes ]; then
        echo "... > !WARNING! You are compiling in DEBUG MODE"
    fi

    if [ -n "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}" ]; then
        echo "... > !WARNING! You are compiling in CODE COVERAGE MODE"
        [[ "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}" != /* ]] && {
            echo "Please provide an absolute path to variable 'MIDAS_COMPILE_CODECOVERAGE_DATAPATH'"
            echo "This value was given: ${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}"
            __status=false
        }
        [ ! -d "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}" ] && mkdir -p ${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}
        [ ! -d "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}" ] && {
            echo "Could not create the directory ${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}"
            __status=false
        }
    fi

    export MIDAS_ABS_LEAFDIR

    echo
    if [ -d "${MIDAS_COMPILE_DIR_BUILD}" ]; then
        if [ "${__fresh_build_directory}" = true ]; then
            echo "Erase the previous build directory: ${MIDAS_COMPILE_DIR_BUILD}"
            rm -rf ${MIDAS_COMPILE_DIR_BUILD}
            echo "Create the build directory: ${MIDAS_COMPILE_DIR_BUILD}"
            mkdir -p ${MIDAS_COMPILE_DIR_BUILD}
        else
            echo "The build directory already exist: ${MIDAS_COMPILE_DIR_BUILD}"
            echo "Leaving it there..."
            echo
        fi
    else
        echo "Create the build directory: ${MIDAS_COMPILE_DIR_BUILD}"
        mkdir -p ${MIDAS_COMPILE_DIR_BUILD}
    fi

    if [ -d ${MIDAS_COMPILE_DIR_MAIN}/${MIDAS_ABS_LEAFDIR} ]; then
        echo "The install directory already exist: ${MIDAS_COMPILE_DIR_MAIN}/${MIDAS_ABS_LEAFDIR}"
        echo "Leaving it there..."
        echo
    else
        echo "Create the install directory: ${MIDAS_COMPILE_DIR_MAIN}/${MIDAS_ABS_LEAFDIR}"
        mkdir ${MIDAS_COMPILE_DIR_MAIN}/${MIDAS_ABS_LEAFDIR}
    fi

    if [ -d "${MIDAS_COMPILE_DIR_BUILD}" ]; then
        __previous_directory=${PWD}

        echo "Moving to ${MIDAS_COMPILE_DIR_BUILD}"
        cd ${MIDAS_COMPILE_DIR_BUILD}

        if [ "${__run_cmake}" = true ]; then
            if [ -f ${MIDAS_COMPILE_DIR_BUILD}/CMakeCache.txt ]; then
                cat <<EOF
It seems 'cmake' has already been run in build directory ${MIDAS_COMPILE_DIR_BUILD}

You can rerun by yourself 'cmake' with the commands
   cd ${MIDAS_COMPILE_DIR_BUILD}
   cmake $(cmake_options) ${MIDAS_SOURCE_DIR}
EOF
            else
                echo "Running cmake $(cmake_options) ${MIDAS_SOURCE_DIR}"
                echo
                cmake $(cmake_options) ${MIDAS_SOURCE_DIR}
                if [ $? -ne 0 ]; then
                    echo
                    echo "CMake failed!"
                    echo
                    __status=false
                    __show_instructions=false
                fi
            fi
        fi ## End of 'if [ "${__run_cmake}" = true ]'

        if [ "${__show_instructions}" = true ]; then
            echo
            if [ "${__run_cmake}" = true ]; then
                cat <<EOF
The build directory has just been prepared by 'cmake'.
You have been moved into the build directory: ${PWD}

You can now compile all of the programs by simply using
EOF
            elif [ "${__run_cmake}" = false ]; then
                if [ ! -f ${MIDAS_COMPILE_DIR_BUILD}/CMakeCache.txt ]; then
                    cat <<EOF
You can run by yourself 'cmake' with the commands
   cd ${MIDAS_COMPILE_DIR_BUILD}
   cmake $(cmake_options) ${MIDAS_SOURCE_DIR}
and build the programs using
EOF
                else
                cat <<EOF
The build directory has already been prepared by 'cmake'.
You have been moved into the build directory: ${PWD}
You can now compile all of the programs by simply using
EOF
                fi
            fi

            cat <<EOF
   make
or just a single program
   make \${program} ## for example 'make var'

You can compile all modules by using
   make midas

Or you can compile a single module by going to 'src/modules' (but only
after you have first compiled all modules or a program at least once)
   cd src/modules
   make gridStateVector_mod.o

Then, to compile any remaining dependent modules and the program and
also create the executable, you need to first go back to the main
build directory
   cd ../..
   make var
EOF
        fi ## End of 'if [ "${__show_instructions}" = true ]'

        if [ "${__cd_to_build_directory}" = false ]; then
            echo "Moving back to ${__previous_directory}"
            cd ${__previous_directory}
        fi

        export MIDAS_COMPILE_FRONTEND
        export MIDAS_COMPILE_JOBNAME
        export MIDAS_COMPILE_VERBOSE
    else
        echo "The build directory ${MIDAS_COMPILE_DIR_BUILD} does not exist"
        __status=false
    fi

    ## this variable is used to prevent sourcing a second time the confiduration
    ## doing so is risky since some environment variables might already be set
    export MIDAS_CONFIG_ALREADY_SOURCED=true

    # config return status
    ${__status}

fi ## End of 'if [ "${__run_cmake}" != stop ]'

# config return status
${__status}
