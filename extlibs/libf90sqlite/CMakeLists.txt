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

# Set include directories for csqlite.c
set(LIBUDFSQLITE_INCLUDE_DIR "" CACHE PATH "Path to libudfsqlite include directory")
if(NOT LIBUDFSQLITE_INCLUDE_DIR)
    message(FATAL_ERROR "LIBUDFSQLITE_INCLUDE_DIR must be set")
endif()
target_include_directories(f90sqlite PRIVATE ${LIBUDFSQLITE_INCLUDE_DIR})

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
