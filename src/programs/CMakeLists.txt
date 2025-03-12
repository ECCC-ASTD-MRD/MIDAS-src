add_executable(energyNorm energyNorm.f90)
target_link_libraries(energyNorm PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
random ${HDF5_Libraries} hpcoperf)

add_executable(letkf letkf.f90)
target_link_libraries(letkf PRIVATE MPI::MPI_Fortran rmn::rmn-ompi rpncomm::rpncomm vgrid::vgrid midas
${SQLITE_Libraries} ${rttov_LIBRARIES} ${HDF5_Libraries} burpmodule random hpcoperf)

set(PROJECT_PROGRAMS energyNorm letkf)

foreach(prog ${PROJECT_PROGRAMS})
  set_target_properties(${prog} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
    LINK_FLAGS "${MPI_LINK_FLAGS}"
  )

  add_custom_command(
      TARGET ${prog}
      POST_BUILD

      COMMAND ${CMAKE_COMMAND} -E echo Copy '${prog}' to ${CMAKE_INSTALL_PREFIX}/midas-${prog}${AddToProgName}-`${CMAKE_SOURCE_DIR}/midas.version`.Abs
      COMMAND cp ${prog} ${CMAKE_INSTALL_PREFIX}/midas-${prog}${AddToProgName}-`${CMAKE_SOURCE_DIR}/midas.version`.Abs

      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endforeach()

add_custom_command(
    OUTPUT update_test_environment.cmake

    COMMAND ${CMAKE_COMMAND} -E echo Executing ${CMAKE_SOURCE_DIR}/src/update_resources.sh
    COMMAND ${CMAKE_SOURCE_DIR}/src/update_resources.sh

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
