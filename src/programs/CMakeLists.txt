add_executable(energyNorm energyNorm.f90)
target_link_libraries(energyNorm PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
random::random ${HDF5_Libraries})

add_executable(letkf letkf.f90)
target_link_libraries(letkf PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

set(PROJECT_PROGRAMS energyNorm letkf)

foreach(prog ${PROJECT_PROGRAMS})
  set_target_properties(${prog} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
    LINK_FLAGS "${MPI_LINK_FLAGS}"
    OUTPUT_NAME "midas-${prog}.Abs"
  )
endforeach()

install(TARGETS energyNorm letkf RUNTIME DESTINATION bin)

