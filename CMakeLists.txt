cmake_minimum_required(VERSION 3.16)

#----- Append EC specific module path
list(APPEND CMAKE_MODULE_PATH $ENV{EC_CMAKE_MODULE_PATH})

include(ec_init)           # Include EC specific cmake utils
ec_git_version()           # Get the version from the git repository
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}" LANGUAGES C Fortran)

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
find_package(SQLite3 ${sqlite3_REQ_VERSION} REQUIRED)
find_package(HDF5 ${hdf5_REQ_VERSION} REQUIRED COMPONENTS C Fortran HL)

# Define libraries bundle
if(DEFINED ENV{RTTOV_VERSION})
    set(rttov_VERSION $ENV{RTTOV_VERSION})
else()
    set(rttov_VERSION "unknown")
endif()

set(rttov_FIND_COMPONENTS rttov_coef_io rttov_hdf rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main)
foreach(component ${rttov_FIND_COMPONENTS})
    find_library(${component}_lib NAME ${component} HINTS ENV LIBRARY_PATH)
    set(rttov_LIBRARIES ${rttov_LIBRARIES} ${${component}_lib})
endforeach()

set(SQLITE_Libraries f90sqlite udfsqlite SQLite::SQLite3)
set(HDF5_Libraries netcdff netcdf hdf5::hdf5_fortran hdf5::hdf5_hl_fortran hdf5::hdf5_hl hdf5::hdf5 z curl)

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
else()
  message(FATAL_ERROR "Compiler is not supported: ${CMAKE_Fortran_COMPILER_ID}")
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

if(DEFINED ENV{MIDAS_COMPILE_DIR_MAIN})
  set(CMAKE_INSTALL_PREFIX "$ENV{MIDAS_COMPILE_DIR_MAIN}/$ENV{MIDAS_ABS_LEAFDIR}" CACHE FILEPATH "Install directory prefix" FORCE)
else()
  message(STATUS "(EC) \${MIDAS_COMPILE_DIR_MAIN} is not set.  Set 'CMAKE_INSTALL_PREFIX' to '${CMAKE_SOURCE_DIR}/compiledir/midas_abs'")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/compiledir/midas_abs" CACHE FILEPATH "Install directory prefix" FORCE)
endif()

message(STATUS "(EC) Sources directory is: ${CMAKE_SOURCE_DIR}")
message(STATUS "(EC) Build directory is: ${CMAKE_BINARY_DIR}")
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

## Build the resources file 'maestro/suites/midas_system_tests/resources/resources.def'
if(DEFINED ENV{ORDENV_PLAT})
  if("$ENV{ORDENV_PLAT}" STREQUAL "rhel-8-icelake-64")
    set(RESOURCE_FILE_EXTENSION "hpcr-u2")
  elseif("$ENV{ORDENV_PLAT}" STREQUAL "rhel-9-graniterapids-64")
    set(RESOURCE_FILE_EXTENSION "hpcr-u3")
  endif()

  ## Build the file only for the two platforms supported
  if (DEFINED RESOURCE_FILE_EXTENSION)
    set(RESOURCES_FILE "${CMAKE_CURRENT_SOURCE_DIR}/maestro/suites/midas_system_tests/resources/resources.def")
    set(RESOURCES_FILE_FOR_PLATFORM "${RESOURCES_FILE}_${RESOURCE_FILE_EXTENSION}")

    if(NOT IS_READABLE ${RESOURCES_FILE_FOR_PLATFORM})
        message(FATAL_ERROR "${RESOURCES_FILE_FOR_PLATFORM} not found!")
    endif()

    if (${CMAKE_BUILD_TYPE} MATCHES Debug)
       set(CHECK_RESULTS_CATCHUP 9)
       set(CLEAN_UNITTEST_CATCHUP 9)
       set(LETKF_GLB15KM_NUM_THREADS 4)
       set(VAR_NITERMAX 10)
    else()
       set(CHECK_RESULTS_CATCHUP 2)
       set(CLEAN_UNITTEST_CATCHUP 2)
       set(LETKF_GLB15KM_NUM_THREADS 2)
       set(VAR_NITERMAX 70)
    endif()

    set(LETKF_GLB10KM_CATCHUP 2)

    set(RESOURCES_HEADER "#########################\n## DO NOT EDIT THIS FILE!\n## Instead, edit ${RESOURCES_FILE_FOR_PLATFORM}")
    set(RESOURCES_FOOTER "## DO NOT EDIT THIS FILE!\n########################")

    configure_file(${RESOURCES_FILE_FOR_PLATFORM} ${RESOURCES_FILE} @ONLY)

    set(MIDAS_TEST_DOT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/maestro/suites/midas_system_tests/abs.dot")
    set(MIDAS_TEST_DOT_FILE_IN "${MIDAS_TEST_DOT_FILE}.in")

    configure_file(${MIDAS_TEST_DOT_FILE_IN} ${MIDAS_TEST_DOT_FILE} @ONLY)
  endif()
endif()
