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

set_target_properties(splitobs PROPERTIES OUTPUT_NAME "midas.splitobs.Abs")
install(TARGETS splitobs RUNTIME DESTINATION bin)
