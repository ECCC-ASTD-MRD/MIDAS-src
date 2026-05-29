file(GLOB PROJECT_MOD_FILES CONFIGURE_DEPENDS *.F90 *.f90)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

add_library(midas ${PROJECT_MOD_FILES})
add_dependencies(midas ${PROJECT_NAME}_build_info)

target_include_directories(midas PUBLIC $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
if(DEFINED ENV{rttov_INSTALLDIR})
  target_include_directories(midas PUBLIC $ENV{rttov_INSTALLDIR}/include)
  target_include_directories(midas PUBLIC $ENV{rttov_INSTALLDIR}/mod)
endif()

target_link_libraries(midas PUBLIC
  rmn::rmn-ompi
  rpncomm::rpncomm
  vgrid::vgrid
  ${hpcoperf}
  random::random
  ${HDF5_Libraries}
  burp::burp
  ${SQLITE_Libraries}
)

if(MKL_SUPPORT)
  target_compile_definitions(midas PRIVATE MKL_SUPPORT)
endif()

if (DGEMMTR_SUPPORT)
  target_compile_definitions(midas PRIVATE DGEMMTR_SUPPORT=${DGEMMTR_SUPPORT})
endif()

# Adding include path for mod files (especially for netcdf and hdf5)
set(FFLAGS "")
if(DEFINED ENV{CPATH})
  string(REPLACE ":" ";" CPATH_LIST $ENV{CPATH})
  foreach(CPATH_DIR ${CPATH_LIST})
    # Only add paths that actually contain .mod files
    # Typically from netcdf-fortran and hdf5
    list(APPEND FFLAGS "-I${CPATH_DIR}")
  endforeach()
endif()

target_compile_options(midas PRIVATE ${FFLAGS})
