file(GLOB PROJECT_MOD_FILES CONFIGURE_DEPENDS *.F90 *.f90)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

add_library(midas ${PROJECT_MOD_FILES})
add_dependencies(midas ${PROJECT_NAME}_build_info)

target_include_directories(midas PUBLIC $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
#target_include_directories(midas PUBLIC $ENV{rttov_INSTALLDIR}/include)
#target_include_directories(midas PUBLIC $ENV{rttov_INSTALLDIR}/mod)

target_link_libraries(midas PUBLIC
  rmn::rmn-ompi
  rpncomm::rpncomm
  vgrid::vgrid
  hpcoperf
)
