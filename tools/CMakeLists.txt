add_executable(splitobs splitobs/splitobs.c splitobs/fstdlib.c splitobs/options.c splitobs/splitobs_sql.c splitobs/splitobs_burp.c splitobs/splitobs_ascii.c splitobs/checkgrid.c splitobs/gzUtils.c)
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
