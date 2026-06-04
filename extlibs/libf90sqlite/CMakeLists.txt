cmake_minimum_required(VERSION 3.16)
project(f90sqlite Fortran C)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Library output directory")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Archive output directory")


# Source files
set(SOURCES
  src/fsqlite.f90
  src/csqlite.c
)

# Create the static library
add_library(f90sqlite STATIC ${SOURCES})

# Link against udfsqlite target (this propagates its include directories)
target_link_libraries(f90sqlite PRIVATE udfsqlite SQLite::SQLite3)

target_include_directories(f90sqlite
  PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
)

set_target_properties(f90sqlite PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

# Set compiler flags
target_compile_options(f90sqlite PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:-Wall -Wextra>
  $<$<COMPILE_LANGUAGE:C>:-Wall -Wextra>
)

# Install library and modules
install(TARGETS f90sqlite
  LIBRARY DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  ARCHIVE DESTINATION ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
)

install(DIRECTORY ${CMAKE_BINARY_DIR}/include/
  DESTINATION ${CMAKE_BINARY_DIR}/include
  FILES_MATCHING PATTERN "*.mod"
)
