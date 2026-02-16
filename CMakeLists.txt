cmake_minimum_required(VERSION 3.16)

#----- Append EC specific module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../rpn/cmake_rpn/modules $ENV{EC_CMAKE_MODULE_PATH})
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

ec_build_info()              # Generate build include file
include(ec_compiler_presets) # Include compiler presets

find_package(MPI REQUIRED)

# To build without OpenMP, you have to add -DWITH_OPENMP=FALSE to the cmake command line
set(WITH_OPENMP TRUE CACHE BOOL "Control whether to use OpenMP")
include(ec_openmp)

find_package(rmn COMPONENTS static REQUIRED)
find_package(vgrid COMPONENTS static REQUIRED)
find_package(rpncomm REQUIRED)
find_package(random COMPONENTS REQUIRED)
find_package(burp COMPONENTS static REQUIRED)
find_package(HDF5 COMPONENTS C Fortran HL)

# Define libraries bundle
if(DEFINED ENV{RTTOV_VERSION})
    set(rttov_VERSION $ENV{RTTOV_VERSION})
else()
    set(rttov_VERSION "unknown")
endif()

set(rttov_LIBRARY_PATH $ENV{rttov_INSTALLDIR}/lib)
set(rttov_FIND_COMPONENTS rttov_coef_io rttov_hdf rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main)
set(rttov_LIBRARIES "") # Initialize list as empty string/list
foreach(component ${rttov_FIND_COMPONENTS})
  find_library(${component} NAMES ${component} REQUIRED
    HINTS PATHS ${rttov_LIBRARY_PATH})
  list(APPEND rttov_LIBRARIES ${${component}})
endforeach()

find_library(hpcoperf NAMES hpcoperf REQUIRED
  HINTS PATHS ENV perftools_LIBRARY_PATH)

set(HDF5_Libraries netcdff netcdf hdf5::hdf5_fortran hdf5::hdf5_hl_fortran hdf5::hdf5_hl hdf5::hdf5 z curl)

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  # Necessary, as compiler presets sets it to -O2
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast-transcendentals -no-prec-div -fpic -ip -no-prec-sqrt")

  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all -fp-speculation=safe -init=snan,arrays")

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl -check noarg_temp_created -no-wrap-margin -warn all -warn noexternal -static-intel")

  # Since we have '-stand f08' in the default compiler options, we replace for '-stand f18'.
  string(REPLACE " -stand f08" " -stand f18" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

  if("$ENV{MIDAS_COMPILE_OPTIMIZE_REPORT}" STREQUAL "yes")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopt-report=5")
  endif()

  if( DEFINED ENV{MIDAS_COMPILE_CODECOVERAGE_DATAPATH} )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -prof-gen=srcpos -prof-dir=$ENV{MIDAS_COMPILE_CODECOVERAGE_DATAPATH}")
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

if(DEFINED ENV{EC_ARCH})
  if(DEFINED ENV{MIDAS_COMPILE_DIR_MAIN})
    set(CMAKE_INSTALL_PREFIX "$ENV{MIDAS_COMPILE_DIR_MAIN}/$ENV{MIDAS_ABS_LEAFDIR}" CACHE FILEPATH "Install directory prefix" FORCE)
  else()
    message(status "(EC) \${MIDAS_COMPILE_DIR_MAIN} is not set.  Make sure you source the src/config.dot.sh file")
  endif()
endif()

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
