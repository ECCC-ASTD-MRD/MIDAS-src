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

foreach(prog ${PROJECT_PROGRAMS})
  set_target_properties(${prog} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
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
