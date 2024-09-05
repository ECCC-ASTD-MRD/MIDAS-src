cmake_minimum_required(VERSION 3.16)

#----- Append EC specific module path
list(APPEND CMAKE_MODULE_PATH $ENV{EC_CMAKE_MODULE_PATH})

include(ec_init)           # Include EC specific cmake utils
ec_git_version()           # Get the version from the git repository
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}" LANGUAGES C Fortran)
# get MIDAS version from midas.version.sh, instead of using cmake_rpn ec_git_version() macro
execute_process(
  COMMAND sh "-c" "${CMAKE_SOURCE_DIR}/midas.version.sh" 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
set(PROJECT_VERSION ${VERSION})
message(STATUS "(EC) ${NAME} VERSION = ${PROJECT_VERSION}")

ec_build_info()              # Generate build include file
include(ec_compiler_presets) # Include compiler presets

include(CTest)
add_custom_target(check COMMAND CTEST_OUTPUT_ON_FAILURE=true ${CMAKE_CTEST_COMMAND})

enable_testing()

find_package(MPI REQUIRED)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

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
set(rttov_FIND_COMPONENTS rttov_parallel rttov_emis_atlas rttov_other rttov_mw_scatt rttov_main rttov_coef_io rttov_hdf)
foreach(component ${rttov_FIND_COMPONENTS})
    find_library(${component}_lib NAME ${component}
       HINTS ENV LIBRARY_PATH)
#       HINTS ${rttov_ROOT}/lib)
    set(rttov_LIBRARIES ${rttov_LIBRARIES} ${${component}_lib})
endforeach()

set(SQLITE_Libraries SQLite::SQLite3 f90sqlite udfsqlite)
set(HDF5_Libraries netcdff netcdf hdf5::hdf5_fortran hdf5::hdf5_hl_fortran hdf5::hdf5_hl hdf5::hdf5 z curl)

# Necessary, as compiler presets sets it to -O2
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  # FOPTMIZ and COMPF are set in src/config.dot.sh
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} $ENV{FOPTMIZ} $ENV{COMPF} -static-intel")
  message(STATUS "FOPTMIZ = $ENV{FOPTMIZ}")
  message(STATUS "COMPF = $ENV{COMPF}")
endif()

# For now, we delete -warn errors from Fortran flags set in the config file,
# as there are errors because of the new flag -stand f08, and the
# compilation fails
string(REPLACE " -warn errors" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

link_libraries("-Wl,--disable-new-dtags")

# Former flags
#set(CMAKE_Fortran_FLAGS "-DAMD64 -DLINUX_X86_64 -DLittle_Endian -DWITH_OpenMP -DWITH_intel -D__INTEL_COMPILER=2021 -O3 -align array32byte -assume byterecl -check noarg_temp_created -convert big_endian -diag-disable=cpu-dispatch -fast-transcendentals -fp-model source -fpe0 -fpic -ip -no-prec-div -no-prec-sqrt -no-wrap-margin -qmkl -qopenmp -static-intel -traceback -warn all -warn errors -warn noexternal -xCORE-AVX512")

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
message(STATUS "(EC) Build directory is: $ENV{MIDAS_COMPILE_DIR_MAIN}")
message(STATUS "(EC) Install directory is: ${CMAKE_INSTALL_PREFIX}")

# for the special name used by MIDAS for programs and tools
if(DEFINED ENV{ORDENV_PLAT})
  set(AddToProgName $ENV{ORDENV_PLAT}-${PROJECT_VERSION})
endif()

add_subdirectory(src/modules)
add_subdirectory(src/programs)
add_subdirectory(tools)
#Don't add dumpBmatrix for now: if needed, it can be compiled with its own Makefile
#add_subdirectory(scripts/convenient_tools/dumpBmatrix)

ec_build_config()      # Create build configuration script
configure_file(config.in ${CMAKE_BINARY_DIR}/${NAME}-config @ONLY)
