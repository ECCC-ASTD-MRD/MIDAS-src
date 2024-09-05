# monitor.Abs
#add_executable(midas.monitor.Abs monitor/monitor.c)

#splitobs.Abs
#DEBUG_OPTIONS = -debug -optc =-w3 =-fp-speculation=safe
#OPT_LEVEL = 2
add_executable(splitobs splitobs/splitobs.c splitobs/fstdlib.c)
target_compile_options(splitobs PRIVATE "-O2")
target_link_libraries(splitobs PRIVATE
  MPI::MPI_C
  SQLite::SQLite3
  burp_c
  rmn::rmn-ompi
)
set_target_properties(splitobs PROPERTIES 
   OUTPUT_NAME midas.splitobs_${AddToProgName}.Abs)
install(TARGETS splitobs DESTINATION ${CMAKE_INSTALL_PREFIX})

if(BUILD_TESTING)
   find_program(BASH_COMMAND bash)
   add_test(NAME test_monitor COMMAND ${BASH_COMMAND} ${CMAKE_SOURCE_DIR}/tools/monitor/test.sh ${CMAKE_BINARY_DIR}/tools/midas.monitor.Abs)
endif()
