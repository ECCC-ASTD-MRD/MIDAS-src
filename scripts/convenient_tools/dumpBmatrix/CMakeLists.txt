add_executable(dumpBmatrix dumpBmatrix.f90)
target_link_libraries(dumpBmatrix PRIVATE rmn::rmn)

install(TARGETS dumpBmatrix DESTINATION scripts)
