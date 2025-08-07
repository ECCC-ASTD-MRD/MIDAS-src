add_executable(adjointTest adjointTest.f90)
target_link_libraries(adjointTest PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(analysisErrorOI analysisErrorOI.f90)
target_link_libraries(analysisErrorOI PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

add_executable(calcStats calcStats.f90)
target_link_libraries(calcStats PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas irc random hpcoperf)

add_executable(dfs dfs.f90)
target_link_libraries(dfs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas irc random hpcoperf
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule)

add_executable(diagBmatrix diagBmatrix.f90)
target_link_libraries(diagBmatrix PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas irc random hpcoperf)

add_executable(diagHBHt diagHBHt.f90)
target_link_libraries(diagHBHt PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(energyNorm energyNorm.f90)
target_link_libraries(energyNorm PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
irc random ${HDF5_Libraries} hpcoperf)

add_executable(ensembleH ensembleH.f90)
target_link_libraries(ensembleH PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(ensPostProcess ensPostProcess.f90)
target_link_libraries(ensPostProcess PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(extractBmatrixFor1Dvar extractBmatrixFor1Dvar.f90)
target_link_libraries(extractBmatrixFor1Dvar PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(genCoeff genCoeff.f90)
target_link_libraries(genCoeff PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(letkf letkf.f90)
target_link_libraries(letkf PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(obsImpact obsImpact.f90)
target_link_libraries(obsImpact PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(obsSelection obsSelection.f90)
target_link_libraries(obsSelection PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(oMinusF oMinusF.f90)
target_link_libraries(oMinusF PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(prepcma prepcma.f90)
target_link_libraries(prepcma PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(pseudoSSTobs pseudoSSTobs.f90)
target_link_libraries(pseudoSSTobs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(randomPert randomPert.f90)
target_link_libraries(randomPert PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas irc random hpcoperf)

add_executable(sstBias sstBias.f90)
target_link_libraries(sstBias PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(sstTrial sstTrial.f90)
target_link_libraries(sstTrial PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(var var.f90)
target_link_libraries(var PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid  midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

add_executable(var1D var1D.f90)
target_link_libraries(var1D PRIVATE MPI::MPI_Fortran  rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule irc random hpcoperf)

set(PROJECT_PROGRAMS adjointTest analysisErrorOI calcStats dfs diagBmatrix diagHBHt energyNorm ensembleH
ensPostProcess extractBmatrixFor1Dvar genCoeff letkf obsImpact obsSelection oMinusF prepcma pseudoSSTobs
randomPert sstBias sstTrial var var1D)

foreach(prog ${PROJECT_PROGRAMS})
  set_target_properties(${prog} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME midas-${prog}_${AddToProgName}.Abs
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
    LINK_FLAGS "${MPI_LINK_FLAGS}")

  add_custom_command(
      TARGET ${prog}
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E echo "Copy midas-${prog}_${AddToProgName}.Abs to ${CMAKE_INSTALL_PREFIX}"
      COMMAND cp -f midas-${prog}_${AddToProgName}.Abs ${CMAKE_INSTALL_PREFIX}/midas-${prog}_${AddToProgName}-`${CMAKE_SOURCE_DIR}/midas.version`.Abs
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endforeach()

install(TARGETS ${PROJECT_PROGRAMS} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})

file(WRITE ${CMAKE_BINARY_DIR}/install_programs_script.cmake "
    execute_process(
        COMMAND ${CMAKE_SOURCE_DIR}/midas.version
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE TRUE_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND ${CMAKE_SOURCE_DIR}/src/update_resources.sh
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src
    )
    foreach(prog ${PROJECT_PROGRAMS})
        set(old_name \"${CMAKE_INSTALL_PREFIX}/midas-\${prog}_${AddToProgName}.Abs\")
        set(new_name \"${CMAKE_INSTALL_PREFIX}/midas-\${prog}_${AddToProgName}-\${TRUE_VERSION}.Abs\")
        file(RENAME  \"\${old_name}\" \"\${new_name}\")
    endforeach()
")

# Install the custom script to run after the targets are installed
install(SCRIPT ${CMAKE_BINARY_DIR}/install_programs_script.cmake)

# install only the compiled programs
add_custom_target(installabs
  COMMAND ${CMAKE_SOURCE_DIR}/midas.version > ${CMAKE_INSTALL_PREFIX}/midas.version.temp

  COMMAND /bin/mkdir -p ${CMAKE_INSTALL_PREFIX}

  COMMAND for program in ${PROJECT_PROGRAMS}\; do progin=${CMAKE_CURRENT_BINARY_DIR}/midas-\$\${program}_${AddToProgName}.Abs\; [ -f \$\${progin} ] && cp \$\${progin} ${CMAKE_INSTALL_PREFIX}/midas-\$\${program}_${AddToProgName}-`cat ${CMAKE_INSTALL_PREFIX}/midas.version.temp`.Abs\; done\; rm -f ${CMAKE_INSTALL_PREFIX}/midas.version.temp
)

set_property(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_CURRENT_BINARY_DIR}/*.Abs")
