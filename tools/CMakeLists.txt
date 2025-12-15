add_executable(splitobs splitobs/splitobs.c splitobs/fstdlib.c)
add_dependencies(splitobs ${PROJECT_NAME}_build_info)

find_package(SQLite3 ${sqlite3_REQ_VERSION} REQUIRED)

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
