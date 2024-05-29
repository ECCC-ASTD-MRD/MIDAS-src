add_executable(adjointTest.Abs adjointTest.f90)
target_link_libraries(adjointTest.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(analysisErrorOI.Abs analysisErrorOI.f90)
target_link_libraries(analysisErrorOI.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module random hpcoperf)

add_executable(calcStats.Abs calcStats.f90)
target_link_libraries(calcStats.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas irc random hpcoperf)

add_executable(diagBmatrix.Abs diagBmatrix.f90)
target_link_libraries(diagBmatrix.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas irc random hpcoperf)

add_executable(diagHBHt.Abs diagHBHt.f90)
target_link_libraries(diagHBHt.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(ensembleH.Abs ensembleH.f90)
target_link_libraries(ensembleH.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(ensPostProcess.Abs ensPostProcess.f90)
target_link_libraries(ensPostProcess.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(extractBmatrixFor1Dvar.Abs extractBmatrixFor1Dvar.f90)
target_link_libraries(extractBmatrixFor1Dvar.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(genCoeff.Abs genCoeff.f90)
target_link_libraries(genCoeff.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(letkf.Abs letkf.f90)
target_link_libraries(letkf.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(obsImpact.Abs obsImpact.f90)
target_link_libraries(obsImpact.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(obsSelection.Abs obsSelection.f90)
target_link_libraries(obsSelection.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(oMinusF.Abs oMinusF.f90)
target_link_libraries(oMinusF.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(prepcma.Abs prepcma.f90)
target_link_libraries(prepcma.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(pseudoSSTobs.Abs pseudoSSTobs.f90)
target_link_libraries(pseudoSSTobs.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(randomPert.Abs randomPert.f90)
target_link_libraries(randomPert.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas irc random hpcoperf)

add_executable(sstBias.Abs sstBias.f90)
target_link_libraries(sstBias.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(sstTrial.Abs sstTrial.f90)
target_link_libraries(sstTrial.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(thinning.Abs thinning.f90)
target_link_libraries(thinning.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(var.Abs var.f90)
target_link_libraries(var.Abs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid  midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

add_executable(var1D.Abs var1D.f90)
target_link_libraries(var1D.Abs PRIVATE MPI::MPI_Fortran  rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas 
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp_module irc random hpcoperf)

set(PROJECT_PROGRAMS adjointTest.Abs analysisErrorOI.Abs calcStats.Abs diagBmatrix.Abs diagHBHt.Abs ensembleH.Abs ensPostProcess.Abs extractBmatrixFor1Dvar.Abs genCoeff.Abs letkf.Abs obsImpact.Abs obsSelection.Abs oMinusF.Abs prepcma.Abs pseudoSSTobs.Abs randomPert.Abs sstBias.Abs sstTrial.Abs thinning.Abs var.Abs var1D.Abs)

foreach(prog ${PROJECT_PROGRAMS})
  get_filename_component(prog_we ${prog} NAME_WE)
  set_target_properties(${prog} PROPERTIES 
    VERSION ${PROJECT_VERSION}
    POSITION_INDEPENDENT_CODE  ON
    OUTPUT_NAME midas-${prog}
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
    LINK_FLAGS "${MPI_LINK_FLAGS}")
  add_custom_command(TARGET ${prog} POST_BUILD
    # Symbolic link to complete name with ORDENV_PLAT variable
    COMMAND ln -sfv midas-${prog} midas-${prog_we}${AddToProgName}.Abs)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/midas-${prog_we}${AddToProgName}.Abs DESTINATION ${CMAKE_INSTALL_PREFIX})
endforeach()

install(TARGETS ${PROJECT_PROGRAMS} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
