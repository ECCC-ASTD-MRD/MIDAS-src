add_executable(adjointTest adjointTest.f90)
target_link_libraries(adjointTest PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(analysisErrorOI analysisErrorOI.f90)
target_link_libraries(analysisErrorOI PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module random hpcoperf)

add_executable(calcStats calcStats.f90)
target_link_libraries(calcStats PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas irc random hpcoperf)

add_executable(diagBmatrix diagBmatrix.f90)
target_link_libraries(diagBmatrix PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas irc random hpcoperf)

add_executable(diagHBHt diagHBHt.f90)
target_link_libraries(diagHBHt PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(ensembleH ensembleH.f90)
target_link_libraries(ensembleH PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(ensPostProcess ensPostProcess.f90)
target_link_libraries(ensPostProcess PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(extractBmatrixFor1Dvar extractBmatrixFor1Dvar.f90)
target_link_libraries(extractBmatrixFor1Dvar PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(genCoeff genCoeff.f90)
target_link_libraries(genCoeff PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(letkf letkf.f90)
target_link_libraries(letkf PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(obsImpact obsImpact.f90)
target_link_libraries(obsImpact PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(obsSelection obsSelection.f90)
target_link_libraries(obsSelection PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(oMinusF oMinusF.f90)
target_link_libraries(oMinusF PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(prepcma prepcma.f90)
target_link_libraries(prepcma PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(pseudoSSTobs pseudoSSTobs.f90)
target_link_libraries(pseudoSSTobs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(randomPert randomPert.f90)
target_link_libraries(randomPert PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas irc random hpcoperf)

add_executable(sstBias sstBias.f90)
target_link_libraries(sstBias PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(sstTrial sstTrial.f90)
target_link_libraries(sstTrial PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(thinning thinning.f90)
target_link_libraries(thinning PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(var var.f90)
target_link_libraries(var PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid  midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(var1D var1D.f90)
target_link_libraries(var1D PRIVATE MPI::MPI_Fortran  rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

set(PROJECT_PROGRAMS adjointTest analysisErrorOI calcStats diagBmatrix diagHBHt ensembleH ensPostProcess extractBmatrixFor1Dvar genCoeff letkf obsImpact obsSelection oMinusF prepcma pseudoSSTobs randomPert sstBias sstTrial thinning var var1D)

foreach(prog ${PROJECT_PROGRAMS})
  set_target_properties(${prog} PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME midas-${prog}_${AddToProgName}.Abs
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endforeach()

install(TARGETS ${PROJECT_PROGRAMS} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
