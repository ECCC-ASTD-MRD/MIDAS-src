# monitor.Abs
#add_executable(midas.monitor.Abs monitor/monitor.c)

#splitobs.Abs
#DEBUG_OPTIONS = -debug -optc =-w3 =-fp-speculation=safe
#OPT_LEVEL = 2
add_executable(splitobs splitobs/splitobs.c splitobs/fstdlib.c)
target_compile_options(splitobs PRIVATE "-O2")
target_link_libraries(splitobs PRIVATE
  SQLite::SQLite3
  burpc.a
  rmn::rmn
)
set_target_properties(splitobs PROPERTIES
   OUTPUT_NAME midas.splitobs_${AddToProgName}.Abs)
install(TARGETS splitobs DESTINATION ${CMAKE_INSTALL_PREFIX})

add_executable(monitor monitor/monitor.c)
# monitor is not associated to a specific platform, thus the shorter suffix
set_target_properties(monitor PROPERTIES
   OUTPUT_NAME midas.monitor.Abs)
install(TARGETS monitor DESTINATION ${CMAKE_INSTALL_PREFIX})

file(WRITE ${CMAKE_BINARY_DIR}/install_splitobs_script.cmake "
    execute_process(
        COMMAND ${CMAKE_SOURCE_DIR}/midas.version.sh
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE TRUE_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
        set(old_name \"${CMAKE_INSTALL_PREFIX}/midas.splitobs_${AddToProgName}.Abs\")
        set(new_name \"${CMAKE_INSTALL_PREFIX}/midas.splitobs_${AddToProgName}-\${TRUE_VERSION}.Abs\")
        file(RENAME  \"\${old_name}\" \"\${new_name}\")
")

file(WRITE ${CMAKE_BINARY_DIR}/install_monitor_script.cmake "
    execute_process(
        COMMAND ${CMAKE_SOURCE_DIR}/midas.version.sh
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE TRUE_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
        set(old_name \"${CMAKE_INSTALL_PREFIX}/midas.monitor.Abs\")
        set(new_name \"${CMAKE_INSTALL_PREFIX}/midas.monitor-\${TRUE_VERSION}.Abs\")
        file(RENAME  \"\${old_name}\" \"\${new_name}\")
")

# Install the custom scripts to run after the targets are installed
install(SCRIPT ${CMAKE_BINARY_DIR}/install_splitobs_script.cmake)
install(SCRIPT ${CMAKE_BINARY_DIR}/install_monitor_script.cmake)

if(BUILD_TESTING)
   find_program(BASH_COMMAND bash)
   add_test(NAME test_monitor COMMAND ${BASH_COMMAND} ${CMAKE_SOURCE_DIR}/tools/monitor/test.sh ${CMAKE_BINARY_DIR}/tools/midas.monitor.Abs)
endif()
