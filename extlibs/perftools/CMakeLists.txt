cmake_minimum_required(VERSION 3.16)
project(perftools C)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

find_package(OpenMP REQUIRED)

add_subdirectory(src/jio)
add_subdirectory(src/timing)

# Create interface library that combines jio and timing
add_library(perftools INTERFACE)
target_link_libraries(perftools INTERFACE jio timing)

install(TARGETS perftools)
