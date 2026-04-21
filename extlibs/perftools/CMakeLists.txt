cmake_minimum_required(VERSION 3.16)
project(perftools C)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

find_package(OpenMP REQUIRED)

add_subdirectory(src/jio)
add_subdirectory(src/timing)
