# monitor.Abs
add_executable(midas.monitor.Abs monitor/monitor.c)
set_target_properties(midas.monitor.Abs PROPERTIES 
   VERSION ${PROJECT_VERSION}
   OUTPUT_NAME midas.monitor.Abs)
add_dependencies(midas.monitor.Abs ${PROJECT_NAME}_build_info)
install(TARGETS midas.monitor.Abs DESTINATION ${CMAKE_INSTALL_PREFIX})
get_filename_component(monitor_we midas.monitor.Abs NAME_WLE)
add_custom_command(TARGET midas.monitor.Abs POST_BUILD
  # Symbolic link to complete name with ORDENV_PLAT variable
  COMMAND ln -sfv midas.monitor.Abs ${monitor_we}${AddToProgName}.Abs)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${monitor_we}${AddToProgName}.Abs DESTINATION ${CMAKE_INSTALL_PREFIX})

#splitobs.Abs
#DEBUG_OPTIONS = -debug -optc =-w3 =-fp-speculation=safe
#OPT_LEVEL = 2
add_executable(midas.splitobs.Abs splitobs/splitobs.c splitobs/fstdlib.c)
target_compile_options(midas.splitobs.Abs PRIVATE "-O2")
target_link_libraries(midas.splitobs.Abs PRIVATE
  MPI::MPI_C
  SQLite::SQLite3
  burp_c
  rmn::rmn-ompi
)
add_dependencies(midas.splitobs.Abs ${PROJECT_NAME}_build_info)
set_target_properties(midas.splitobs.Abs PROPERTIES 
   VERSION ${PROJECT_VERSION}
   OUTPUT_NAME midas.splitobs.Abs)
install(TARGETS midas.splitobs.Abs DESTINATION ${CMAKE_INSTALL_PREFIX})
get_filename_component(splitobs_we midas.splitobs.Abs NAME_WLE)
add_custom_command(TARGET midas.splitobs.Abs POST_BUILD
  # Symbolic link to complete name with ORDENV_PLAT variable
  COMMAND ln -sfv midas.splitobs.Abs ${splitobs_we}${AddToProgName}.Abs)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${splitobs_we}${AddToProgName}.Abs DESTINATION ${CMAKE_INSTALL_PREFIX})

# scripts
install(PROGRAMS 
    findTrials/midas.findTrials
    DESTINATION scripts)
install(PROGRAMS 
    midas_scripts/midas.check_ensemble 
    midas_scripts/midas.launch 
    midas_scripts/midas.mpi 
    midas_scripts/midas.mpi_barrier 
    midas_scripts/midas.mpirun 
    midas_scripts/midas.reunir_obs 
    midas_scripts/midas.reunir_obs_mpi 
    midas_scripts/midas.splitobs 
    midas_scripts/midas.tripotenml 
    DESTINATION scripts)
 
if(BUILD_TESTING)
   find_program(BASH_COMMAND bash)
   add_test(NAME test_monitor COMMAND ${BASH_COMMAND} ${CMAKE_SOURCE_DIR}/tools/monitor/test.sh ${CMAKE_BINARY_DIR}/tools/midas.monitor.Abs)
endif()
