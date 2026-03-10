add_executable(adjointTest adjointTest.f90)
target_link_libraries(adjointTest PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(analysisErrorOI analysisErrorOI.f90)
target_link_libraries(analysisErrorOI PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(calcStats calcStats.f90)
target_link_libraries(calcStats PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas random hpcoperf)

add_executable(dfs dfs.f90)
target_link_libraries(dfs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas random hpcoperf
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule)

add_executable(diagBmatrix diagBmatrix.f90)
target_link_libraries(diagBmatrix PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas random hpcoperf)

add_executable(diagHBHt diagHBHt.f90)
target_link_libraries(diagHBHt PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(energyNorm energyNorm.f90)
target_link_libraries(energyNorm PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
random ${HDF5_Libraries} hpcoperf)

add_executable(ensembleH ensembleH.f90)
target_link_libraries(ensembleH PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(ensPostProcess ensPostProcess.f90)
target_link_libraries(ensPostProcess PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(extractBmatrixFor1Dvar extractBmatrixFor1Dvar.f90)
target_link_libraries(extractBmatrixFor1Dvar PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(genCoeff genCoeff.f90)
target_link_libraries(genCoeff PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(letkf letkf.f90)
target_link_libraries(letkf PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(obsImpact obsImpact.f90)
target_link_libraries(obsImpact PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(obsSelection obsSelection.f90)
target_link_libraries(obsSelection PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(oMinusF oMinusF.f90)
target_link_libraries(oMinusF PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(prepcma prepcma.f90)
target_link_libraries(prepcma PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(pseudoSSTobs pseudoSSTobs.f90)
target_link_libraries(pseudoSSTobs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(randomPert randomPert.f90)
target_link_libraries(randomPert PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas random hpcoperf)

add_executable(sstBias sstBias.f90)
target_link_libraries(sstBias PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(sstTrial sstTrial.f90)
target_link_libraries(sstTrial PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(var var.f90)
target_link_libraries(var PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid  midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(var1D var1D.f90)
target_link_libraries(var1D PRIVATE MPI::MPI_Fortran  rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

set(PROJECT_PROGRAMS adjointTest analysisErrorOI calcStats dfs diagBmatrix diagHBHt energyNorm ensembleH
ensPostProcess extractBmatrixFor1Dvar genCoeff letkf obsImpact obsSelection oMinusF prepcma pseudoSSTobs
randomPert sstBias sstTrial var var1D)

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

add_custom_command(
    OUTPUT update_test_environment.cmake

    COMMAND ${CMAKE_COMMAND} -E echo Executing ${CMAKE_SOURCE_DIR}/maestro/suites/midas_system_tests/update_abs.dot.sh ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_SOURCE_DIR}/maestro/suites/midas_system_tests/update_abs.dot.sh ${CMAKE_SOURCE_DIR}

    COMMENT "Preparing the testing environment"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
add_custom_target(prepare_test DEPENDS update_test_environment.cmake)

set_property(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_CURRENT_BINARY_DIR}/*.Abs"
)
