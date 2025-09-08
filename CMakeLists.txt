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
    find_library(${component}_lib NAME ${component}
       HINTS ENV LIBRARY_PATH)
    set(rttov_LIBRARIES ${rttov_LIBRARIES} ${${component}_lib})
endforeach()

set(SQLITE_Libraries SQLite::SQLite3 f90sqlite udfsqlite)
set(HDF5_Libraries netcdff netcdf hdf5::hdf5_fortran hdf5::hdf5_hl_fortran hdf5::hdf5_hl hdf5::hdf5 z curl)

# Necessary, as compiler presets sets it to -O2
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast-transcendentals -no-prec-div -fpic -ip -no-prec-sqrt")

# 'MIDAS_COMPILE_COMPF' is set in 'src/config.dot.sh'
if(DEFINED ENV{MIDAS_COMPILE_COMPF})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} $ENV{MIDAS_COMPILE_COMPF}")
  message(STATUS "Appending CMAKE_Fortran_FLAGS with options from MIDAS_COMPILE_COMPF = $ENV{MIDAS_COMPILE_COMPF}")
else()
  message(WARNING "The environment variable 'MIDAS_COMPILE_COMPF' does not exist.  It should be defined in 'src/config.dot.sh'.")
endif()

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -static-intel")
endif()

link_libraries("-Wl,--disable-new-dtags")

# Printing flags
message(STATUS "(EC) CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
message(STATUS "(EC) CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")

if (${CMAKE_BUILD_TYPE} MATCHES Debug )
  message(STATUS "(EC) CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}")
  message(STATUS "(EC) CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}")
endif()

if(DEFINED ENV{MIDAS_COMPILE_DIR_MAIN})
  set(CMAKE_INSTALL_PREFIX "$ENV{MIDAS_COMPILE_DIR_MAIN}/$ENV{MIDAS_ABS_LEAFDIR}" CACHE FILEPATH "Install directory prefix" FORCE)
else()
  message(FATAL_ERROR "(EC) \${MIDAS_COMPILE_DIR_MAIN} is not set.  Make sure you source the src/config.dot.sh file")
endif()

message(STATUS "(EC) Sources directory is: ${CMAKE_SOURCE_DIR}")
message(STATUS "(EC) Build directory is: $ENV{MIDAS_COMPILE_DIR_BUILD}")
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
