add_executable(splitobs splitobs/splitobs.c splitobs/fstdlib.c)
add_dependencies(splitobs ${PROJECT_NAME}_build_info)
target_compile_options(splitobs PRIVATE "-O2")
target_link_libraries(splitobs PRIVATE
  SQLite::SQLite3
  burpc.a
  rmn::rmn
)
set_target_properties(splitobs PROPERTIES
   OUTPUT_NAME midas.splitobs_${AddToProgName}.Abs)

add_custom_command(
    TARGET splitobs
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Copy midas.splitobs_${AddToProgName}.Abs to ${CMAKE_INSTALL_PREFIX}"
    COMMAND cp -f midas.splitobs_${AddToProgName}.Abs ${CMAKE_INSTALL_PREFIX}/midas.splitobs_${AddToProgName}-`${CMAKE_SOURCE_DIR}/midas.version`.Abs
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(monitor monitor/monitor.c)
# monitor is not associated to a specific platform, thus the shorter suffix
set_target_properties(monitor PROPERTIES
   OUTPUT_NAME midas.monitor.Abs)

add_custom_command(
    TARGET monitor
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Copy midas.monitor.Abs to ${CMAKE_INSTALL_PREFIX}"
    COMMAND cp -f midas.monitor.Abs ${CMAKE_INSTALL_PREFIX}/midas.monitor-`${CMAKE_SOURCE_DIR}/midas.version`.Abs
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

if(BUILD_TESTING)
   find_program(BASH_COMMAND bash)
   add_test(NAME test_monitor COMMAND ${BASH_COMMAND} ${CMAKE_SOURCE_DIR}/tools/monitor/test.sh ${CMAKE_BINARY_DIR}/tools/midas.monitor.Abs)
endif()
