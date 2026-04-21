cmake_minimum_required(VERSION 3.16)
project(udfsqlite C)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


# Find required libraries
find_package(MPI COMPONENTS C)
if(NOT MPI_C_FOUND)
  message(FATAL_ERROR "MPI C compiler not found. You should `. r.load.dot rpn/code-tools/...`")
endif()

find_package(rmn REQUIRED)
find_library(LIBZ NAMES z REQUIRED)
find_library(LIBM NAMES m REQUIRED)

add_library(udfsqlite STATIC src/udfsqlite.c)

target_compile_definitions(udfsqlite PRIVATE SQLITE_CORE)

target_link_libraries(udfsqlite PRIVATE
  rmn::rmn
  MPI::MPI_C
  ${LIBZ}
  ${LIBM}
)

install(TARGETS udfsqlite
  ARCHIVE DESTINATION ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
)
