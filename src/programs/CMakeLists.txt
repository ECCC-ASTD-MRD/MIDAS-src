add_executable(adjointTest adjointTest.f90)
target_link_libraries(adjointTest PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(analysisErrorOI analysisErrorOI.f90)
target_link_libraries(analysisErrorOI PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(calcStats calcStats.f90)
target_link_libraries(calcStats PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas random::random)

add_executable(createLamTemplateGrids createLamTemplateGrids.f90)
target_link_libraries(createLamTemplateGrids PRIVATE rmn::rmn-ompi vgrid::vgrid midas random::random)

add_executable(dfs dfs.f90)
target_link_libraries(dfs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas random::random
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp)

add_executable(diagBmatrix diagBmatrix.f90)
target_link_libraries(diagBmatrix PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas random::random)

add_executable(diagHBHt diagHBHt.f90)
target_link_libraries(diagHBHt PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(energyNorm energyNorm.f90)
target_link_libraries(energyNorm PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
random::random ${HDF5_Libraries})

add_executable(ensembleH ensembleH.f90)
target_link_libraries(ensembleH PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(ensPostProcess ensPostProcess.f90)
target_link_libraries(ensPostProcess PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(extractBmatrixFor1Dvar extractBmatrixFor1Dvar.f90)
target_link_libraries(extractBmatrixFor1Dvar PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(genCoeff genCoeff.f90)
target_link_libraries(genCoeff PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(letkf letkf.f90)
target_link_libraries(letkf PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(obsImpact obsImpact.f90)
target_link_libraries(obsImpact PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(obsSelection obsSelection.f90)
target_link_libraries(obsSelection PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(oMinusF oMinusF.f90)
target_link_libraries(oMinusF PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(prepcma prepcma.f90)
target_link_libraries(prepcma PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(pseudoOceanIceObs pseudoOceanIceObs.f90)
target_link_libraries(pseudoOceanIceObs PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(randomPert randomPert.f90)
target_link_libraries(randomPert PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas random::random)

add_executable(sstBias sstBias.f90)
target_link_libraries(sstBias PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(sstTrial sstTrial.f90)
target_link_libraries(sstTrial PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(var var.f90)
target_link_libraries(var PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid  midas
  ${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

add_executable(var1D var1D.f90)
target_link_libraries(var1D PRIVATE MPI::MPI_Fortran  rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
  ${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burp::burp random::random)

set(PROJECT_PROGRAMS adjointTest analysisErrorOI calcStats createLamTemplateGrids dfs diagBmatrix diagHBHt
  energyNorm ensembleH ensPostProcess extractBmatrixFor1Dvar genCoeff letkf obsImpact obsSelection oMinusF
  prepcma pseudoOceanIceObs randomPert sstBias sstTrial var var1D)

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
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

set_property(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_CURRENT_BINARY_DIR}/*.Abs"
)
