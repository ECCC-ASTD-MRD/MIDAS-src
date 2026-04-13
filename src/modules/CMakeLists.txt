file(GLOB PROJECT_MOD_FILES CONFIGURE_DEPENDS *.F90 *.f90)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

add_library(midas ${PROJECT_MOD_FILES})
add_dependencies(midas ${PROJECT_NAME}_build_info)

target_include_directories(midas PUBLIC $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
if(DEFINED ENV{rttov_INSTALLDIR})
  target_include_directories(midas PUBLIC $ENV{rttov_INSTALLDIR}/include)
  target_include_directories(midas PUBLIC $ENV{rttov_INSTALLDIR}/mod)
endif()
target_include_directories(midas PUBLIC $ENV{f90sqlite_INSTALLDIR}/include)
if(DEFINED ENV{hdf5_INSTALLDIR})
  target_include_directories(midas PUBLIC $ENV{hdf5_INSTALLDIR}/include)
endif()
if(DEFINED ENV{netcdf_INSTALLDIR})
  target_include_directories(midas PUBLIC $ENV{netcdf_INSTALLDIR}/include)
endif()

target_link_libraries(midas PUBLIC
  rmn::rmn-ompi
  rpncomm::rpncomm
  vgrid::vgrid
  ${hpcoperf}
  random::random
  ${HDF5_Libraries}
  burp::burp
)

if(MKL_SUPPORT)
  target_compile_definitions(midas PRIVATE MKL_SUPPORT)
endif()

if (DGEMMTR_SUPPORT)
  target_compile_definitions(midas PRIVATE DGEMMTR_SUPPORT=${DGEMMTR_SUPPORT})
endif()
