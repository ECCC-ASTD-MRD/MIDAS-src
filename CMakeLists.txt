cmake_minimum_required(VERSION 3.16)

option(SQLITE_FORTRAN_SUPPORT "Enable SQLite support in fortran code (links with 'f90sqlite' and 'udfsqlite' internal ECCC librairies)" ON)
option(MKL_SUPPORT "Use Intel MKL mathematical library for blas and lapack routines" OFF)
option(DGEMMTR_SUPPORT "Use the BLAS routine for matrix multiplication if the resulting matrix is known to be symmetric (Intel MKL is implementing this routine as 'dgemmt'" OFF)

#----- Append EC specific module path
if(DEFINED ENV{EC_CMAKE_MODULE_PATH})
  list(APPEND CMAKE_MODULE_PATH $ENV{EC_CMAKE_MODULE_PATH})
endif()

include(ec_init)           # Include EC specific cmake utils
ec_git_version()           # Get the version from the git repository
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}" LANGUAGES C Fortran)

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  message(STATUS "Activating MKL support")
  set(MKL_SUPPORT ON)
  message(STATUS "Use 'dgemmt' as DGEMMTR")
  set(DGEMMTR_SUPPORT "dgemmt")
endif()

## Find the 'base' version of the project (the latest tag of the project)
execute_process(
  COMMAND sh "-c" "${CMAKE_SOURCE_DIR}/midas.version --abbrev=0 HEAD"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE MIDAS_BASE_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

## This variable must be set before loading 'ec_compiler_presets'
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
  set(STRICT TRUE)
endif()

ec_build_info()              # Generate build include file
include(ec_compiler_presets) # Include compiler presets

find_package(MPI REQUIRED)

# To build without OpenMP, you have to add -DWITH_OPENMP=FALSE to the cmake command line
set(WITH_OPENMP TRUE CACHE BOOL "Control whether to use OpenMP")
include(ec_openmp)

find_package(rmn ${rmn_REQ_VERSION} COMPONENTS static REQUIRED)
find_package(vgrid ${vgrid_REQ_VERSION} COMPONENTS static REQUIRED)
find_package(rpncomm ${rpncomm_REQ_VERSION} REQUIRED)
find_package(random ${random_REQ_VERSION} COMPONENTS REQUIRED)
find_package(burp ${burp_REQ_VERSION} COMPONENTS static REQUIRED)
find_package(SQLite3 ${sqlite3_REQ_VERSION} REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C Fortran HL)
if(NOT MKL_SUPPORT)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

# Define libraries bundle
if(DEFINED ENV{RTTOV_VERSION})
    set(rttov_VERSION $ENV{RTTOV_VERSION})
else()
    set(rttov_VERSION "unknown")
endif()

set(rttov_FIND_COMPONENTS rttov_coef_io rttov_hdf rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main)
set(rttov_LIBRARIES "") # Initialize list as empty string/list
foreach(component ${rttov_FIND_COMPONENTS})
  if(DEFINED ENV{rttov_INSTALLDIR})
    find_library(${component}_lib NAMES ${component} REQUIRED HINTS PATHS $ENV{rttov_INSTALLDIR}/lib)
  else()
    find_library(${component}_lib NAME ${component} HINTS ENV LIBRARY_PATH)
  endif()
  list(APPEND rttov_LIBRARIES ${${component}_lib})
endforeach()

find_library(hpcoperf NAMES hpcoperf REQUIRED HINTS PATHS ENV LIBRARY_PATH)

if (SQLITE_FORTRAN_SUPPORT)
  set(SQLITE_Libraries f90sqlite udfsqlite SQLite::SQLite3)
else()
  set(SQLITE_Libraries SQLite::SQLite3)
endif()

find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)

find_package(NetCDF OPTIONAL_COMPONENTS Fortran)
if (NETCDF_Fortran_FOUND)
  set(HDF5_Libraries NetCDF::NetCDF_Fortran)
else()
  set(HDF5_Libraries netcdff netcdf)
endif()

list(APPEND HDF5_Libraries hdf5::hdf5_fortran hdf5::hdf5_hl_fortran hdf5::hdf5_hl hdf5::hdf5 ZLIB::ZLIB CURL::libcurl)

if (NOT MKL_SUPPORT)
  list(APPEND HDF5_Libraries ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  ## If code codecoverage is activated, we use '-O0'
  if( DEFINED ENV{MIDAS_COMPILE_CODECOVERAGE_DATAPATH} )
    set(CMAKE_Fortran_FLAGS "-O0 -prof-gen=srcpos -prof-dir=$ENV{MIDAS_COMPILE_CODECOVERAGE_DATAPATH}")
  else()
    # Necessary, as compiler sets it to '-O2' by default
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast-transcendentals -no-prec-div -fpic -ip -no-prec-sqrt")
  endif()

  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all -fp-speculation=safe -init=snan,arrays")

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl -check noarg_temp_created -no-wrap-margin -warn all -warn noexternal -static-intel")

  # Since we have '-stand f08' in the default compiler options, we replace for '-stand f18'.
  string(REPLACE " -stand f08" " -stand f18" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

  if("$ENV{MIDAS_COMPILE_OPTIMIZE_REPORT}" STREQUAL "yes")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopt-report=5")
  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
  ## If code codecoverage is activated, we use '-O0' because the compiler sets it to '-O2' by default
  if( DEFINED ENV{MIDAS_COMPILE_CODECOVERAGE_DATAPATH} )
    set(CMAKE_Fortran_FLAGS "-O0 -prof-gen=srcpos -prof-dir=$ENV{MIDAS_COMPILE_CODECOVERAGE_DATAPATH}")
  endif()

  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all -fp-speculation=safe -init=snan,arrays")

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl -check noarg_temp_created -no-wrap-margin -warn all -warn noexternal")

  # Since we have '-stand f08' in the default compiler options, we replace for '-stand f18'.
  string(REPLACE " -stand f08" " -stand f18" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

  if("$ENV{MIDAS_COMPILE_OPTIMIZE_REPORT}" STREQUAL "yes")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopt-report=5")
  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  ## This variable is used to add some libraries at link time in 'src/programs/CMakeLists.txt'
  set(MKLROOT "/fs/ssm/main/opt/intelcomp/master/inteloneapi_2022.1.2_multi/oneapi/mkl/latest")

  if("$ENV{MIDAS_COMPILE_OPTIMIZE_REPORT}" STREQUAL "yes")
    message(WARNING "(EC) The option 'MIDAS_COMPILE_OPTIMIZE_REPORT' is not supported for this compiler")
  endif()

  if( DEFINED ENV{MIDAS_COMPILE_CODECOVERAGE_DATAPATH} )
    message(WARNING "(EC) The option 'MIDAS_COMPILE_CODECOVERAGE_DATAPATH' is not supported for this compiler")
  endif()

  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wno-conversion -Wno-compare-reals -Wno-unused-function -Wno-tabs")
else()
  if(DEFINED ENV{EC_ARCH})
    message(FATAL_ERROR "This compiler is not supported: ${CMAKE_Fortran_COMPILER_ID}")
  else()
    message(WARNING "This compiler is not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
  endif()
endif()

link_libraries("-Wl,--disable-new-dtags")

# 'MIDAS_COMPILE_COMPF_GLOBAL' is set in 'src/config.dot.sh'
if(DEFINED ENV{MIDAS_COMPILE_COMPF_GLOBAL})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} $ENV{MIDAS_COMPILE_COMPF_GLOBAL}")
  message(STATUS "Appending CMAKE_Fortran_FLAGS with options from MIDAS_COMPILE_COMPF = $ENV{MIDAS_COMPILE_COMPF}")
endif()

# Printing flags
message(STATUS "(EC) CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
message(STATUS "(EC) CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
  message(STATUS "(EC) CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}")
  message(STATUS "(EC) CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(DEFINED ENV{MIDAS_COMPILE_DIR_MAIN})
    set(CMAKE_INSTALL_PREFIX "$ENV{MIDAS_COMPILE_DIR_MAIN}/$ENV{MIDAS_ABS_LEAFDIR}" CACHE FILEPATH "Install directory prefix" FORCE)
  else()
    if(DEFINED ENV{EC_ARCH})
      message(STATUS "(EC) \${MIDAS_COMPILE_DIR_MAIN} is not set.  Set 'CMAKE_INSTALL_PREFIX' to '${CMAKE_SOURCE_DIR}/compiledir/midas_abs'")
      set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/compiledir/midas_abs" CACHE FILEPATH "Install directory prefix" FORCE)
    else()
      message(STATUS "Cannot guess and set 'CMAKE_INSTALL_PREFIX' from environment")
    endif()
  endif()
endif()

message(STATUS "(EC) Sources directory is: ${CMAKE_SOURCE_DIR}")
message(STATUS "(EC) Build   directory is: ${CMAKE_BINARY_DIR}")
message(STATUS "(EC) Install directory is: ${CMAKE_INSTALL_PREFIX}")

# for the special name used by MIDAS for programs and tools
if(DEFINED ENV{ORDENV_PLAT})
  set(AddToProgName _$ENV{ORDENV_PLAT})
endif()

add_subdirectory(src/modules)
add_subdirectory(src/programs)
add_subdirectory(tools)
#Don't add dumpBmatrix for now: if needed, it can be compiled with its own Makefile
#add_subdirectory(scripts/convenient_tools/dumpBmatrix)

## The 'cmake_rpn' macro 'ec_build_config' contains an install target
## which generates the 'midas-config' script.
ec_build_config()

## External libraries on demand compilation
## to compile them, define the variable MIDAS_COMPILE_EXTLIBS to include the libraries needed
## Exp: MIDAS_COMPILE_EXTLIBS=libudsqlite,libf90sqlite
include(ExternalProject)

if(DEFINED ENV{MIDAS_COMPILE_EXTLIBS})
  set(EXTLIBS_TO_BUILD $ENV{MIDAS_COMPILE_EXTLIBS})

  if(EXTLIBS_TO_BUILD MATCHES "f90sqlite")
    message(STATUS "(EC) Adding libf90sqlite external library target")
    ExternalProject_Add(libf90sqlite
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/extlibs/libf90sqlite/src
      BINARY_DIR ${CMAKE_BINARY_DIR}/extlibs/libf90sqlite-build
      DOWNLOAD_COMMAND ""
      UPDATE_COMMAND ""
      CONFIGURE_COMMAND ""
      BUILD_COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}/extlibs/libf90sqlite/src && make lib"
      INSTALL_COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}/extlibs/libf90sqlite/src && mv libf90sqlite.a ${CMAKE_BINARY_DIR}/extlibs/libf90sqlite-build"
    )
  endif()

  if(EXTLIBS_TO_BUILD MATCHES "udfsqlite")
    message(STATUS "(EC) Adding libudfsqlite external library target")
    ExternalProject_Add(libudfsqlite
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/extlibs/libudfsqlite/src
      BINARY_DIR ${CMAKE_BINARY_DIR}/extlibs/libudfsqlite-build
      DOWNLOAD_COMMAND ""
      UPDATE_COMMAND ""
      CONFIGURE_COMMAND ""
      BUILD_COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}/extlibs/libudfsqlite/src && make"
      INSTALL_COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}/extlibs/libudfsqlite/src && make install DESTDIR=${CMAKE_BINARY_DIR}/extlibs/libudfsqlite-build"
    )
  endif()

  if(EXTLIBS_TO_BUILD MATCHES "perftools")
    message(STATUS "(EC) Adding perftools external library target")
    ExternalProject_Add(perftools
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/extlibs/perftools/src
      BINARY_DIR ${CMAKE_BINARY_DIR}/extlibs/perftools-build
      DOWNLOAD_COMMAND ""
      UPDATE_COMMAND ""
      CONFIGURE_COMMAND ""
      BUILD_COMMAND make -C ${CMAKE_SOURCE_DIR}/extlibs/perftools/src
      INSTALL_COMMAND make -C ${CMAKE_SOURCE_DIR}/extlibs/perftools/src install DESTDIR=${CMAKE_BINARY_DIR}/extlibs/perftools-build
    )
  endif()
endif()

## Build the resources file 'maestro/suites/midas_system_tests/resources/resources.def'
if(DEFINED ENV{ORDENV_PLAT})
  if("$ENV{ORDENV_PLAT}" STREQUAL "rhel-8-icelake-64")
    set(RESOURCE_FILE_EXTENSION "hpcr-u2.in")
  elseif("$ENV{ORDENV_PLAT}" STREQUAL "rhel-9-graniterapids-64")
    set(RESOURCE_FILE_EXTENSION "hpcr-u3.in")
  endif()

  ## Build the file only for the two platforms supported
  if (DEFINED RESOURCE_FILE_EXTENSION)
    set(RESOURCES_FILE "${CMAKE_CURRENT_SOURCE_DIR}/maestro/suites/midas_system_tests/resources/resources.def")
    set(RESOURCES_FILE_FOR_PLATFORM "${RESOURCES_FILE}_${RESOURCE_FILE_EXTENSION}")

    if(NOT IS_READABLE ${RESOURCES_FILE_FOR_PLATFORM})
        message(FATAL_ERROR "${RESOURCES_FILE_FOR_PLATFORM} not found!")
    endif()

    if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND ${CMAKE_BUILD_TYPE} MATCHES Debug)
       set(CHECK_RESULTS_CATCHUP 9)
       set(CLEAN_UNITTEST_CATCHUP 9)
       set(LETKF_GLB15KM_NUM_THREADS 4)
       set(VAR_NITERMAX 10)
       set(VAR1D_IDEALIZED_NITERMAX 10)
       set(LETKF_GLB10KM_CATCHUP 9)
    elseif (${CMAKE_BUILD_TYPE} MATCHES Debug)
       set(CHECK_RESULTS_CATCHUP 9)
       set(CLEAN_UNITTEST_CATCHUP 9)
       set(LETKF_GLB15KM_NUM_THREADS 4)
       set(VAR_NITERMAX 10)
       set(LETKF_GLB10KM_CATCHUP 2)
    else()
       set(CHECK_RESULTS_CATCHUP 2)
       set(CLEAN_UNITTEST_CATCHUP 2)
       set(LETKF_GLB15KM_NUM_THREADS 2)
       set(VAR_NITERMAX 70)
       set(VAR1D_IDEALIZED_NITERMAX 80)
       set(LETKF_GLB10KM_CATCHUP 2)
    endif()

    set(RESOURCES_HEADER "#########################\n## DO NOT EDIT THIS FILE!\n## Instead, edit ${RESOURCES_FILE_FOR_PLATFORM}")
    set(RESOURCES_FOOTER "## DO NOT EDIT THIS FILE!\n########################")

    configure_file(${RESOURCES_FILE_FOR_PLATFORM} ${RESOURCES_FILE} @ONLY)

    set(MIDAS_TEST_DOT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/maestro/suites/midas_system_tests/abs.dot")
    set(MIDAS_TEST_DOT_FILE_IN "${MIDAS_TEST_DOT_FILE}.in")

    configure_file(${MIDAS_TEST_DOT_FILE_IN} ${MIDAS_TEST_DOT_FILE} @ONLY)
  endif()
endif()
