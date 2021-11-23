# Based on tutorial at https://raymii.org/s/tutorials/Cpp_project_setup_with_cmake_and_unit_tests.html

cmake_minimum_required(VERSION 3.10)
project(FourthOrderMicromagneticSolver)

set(CMAKE_CXX_STANDARD 14)

# Required for tests to be able to find header files
include_directories(src)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(lib/googletest)
