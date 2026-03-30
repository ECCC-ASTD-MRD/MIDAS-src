add_executable(splitobs splitobs/splitobs.c splitobs/fstdlib.c)
add_dependencies(splitobs ${PROJECT_NAME}_build_info)

target_link_libraries(splitobs PRIVATE
  SQLite::SQLite3
  burp::burpc-static
  rmn::rmn
)
if (DEFINED CMAKE_INSTALL_RPATH)
  set_target_properties(splitobs PROPERTIES
     BUILD_RPATH ${CMAKE_INSTALL_RPATH}
  )
endif()

add_custom_command(
    TARGET splitobs
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Copy 'splitobs' to ${CMAKE_INSTALL_PREFIX}/midas.splitobs${AddToProgName}-`${CMAKE_SOURCE_DIR}/midas.version`.Abs"
    COMMAND cp splitobs ${CMAKE_INSTALL_PREFIX}/midas.splitobs${AddToProgName}-`${CMAKE_SOURCE_DIR}/midas.version`.Abs
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(monitor monitor/monitor.c)

add_custom_command(
    TARGET monitor
    POST_BUILD
    # monitor is not associated to a specific platform, thus the shorter suffix
    COMMAND ${CMAKE_COMMAND} -E echo "Copy 'monitor' to ${CMAKE_INSTALL_PREFIX}/midas.monitor-`${CMAKE_SOURCE_DIR}/midas.version`.Abs"
    COMMAND cp monitor ${CMAKE_INSTALL_PREFIX}/midas.monitor-`${CMAKE_SOURCE_DIR}/midas.version`.Abs
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
