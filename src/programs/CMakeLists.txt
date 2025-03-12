add_executable(energyNorm energyNorm.f90)
target_link_libraries(energyNorm PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
random::random ${HDF5_Libraries})

add_executable(letkf letkf.f90)
target_link_libraries(letkf PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(ensPostProcess ensPostProcess.f90)
target_link_libraries(ensPostProcess PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} random::random)

set(PROJECT_PROGRAMS energyNorm letkf ensPostProcess)

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND MKL_SUPPORT)
   ## The variable 'MKLROOT' defined in '$(git rev-parse --show-toplevel)/CMakeLists.txt'
   ## From 'https://hpc.ncsu.edu/Software/Libraries.php?app=MKL'
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread")
endif()

foreach(prog ${PROJECT_PROGRAMS})
  target_link_directories(${prog} PUBLIC $ENV{rttov_INSTALLDIR}/lib)
  if(USE_LOCAL_EXTLIBS)
    target_link_libraries(${prog} PUBLIC perftools f90sqlite)
  endif()
endforeach()

list(APPEND FFLAGS "${MPI_COMPILE_FLAGS}")

foreach(prog ${PROJECT_PROGRAMS})

  set_target_properties(${prog} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    COMPILE_FLAGS "${FFLAGS}"
    LINK_FLAGS "${MPI_LINK_FLAGS}"
  )

  add_custom_command(
      TARGET ${prog}
      POST_BUILD

      COMMAND ${CMAKE_COMMAND} -E echo Move '${prog}' to ${CMAKE_INSTALL_PREFIX}/midas-${prog}${AddToProgName}-`${CMAKE_SOURCE_DIR}/midas.version`.Abs
      COMMAND mv ${prog} ${CMAKE_INSTALL_PREFIX}/midas-${prog}${AddToProgName}-`${CMAKE_SOURCE_DIR}/midas.version`.Abs

      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endforeach()

install(TARGETS energyNorm letkf ensPostProcess RUNTIME DESTINATION bin)
