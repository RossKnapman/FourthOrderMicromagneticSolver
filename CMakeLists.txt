# Based on tutorial at https://raymii.org/s/tutorials/Cpp_project_setup_with_cmake_and_unit_tests.html

cmake_minimum_required(VERSION 3.10)
project(FourthOrderMicromagneticSolver)

set(CMAKE_CXX_STANDARD 14)

# Required for tests to be able to find header files
include_directories(src)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(lib/googletest)

#cmake_minimum_required(VERSION 3.14)
#project(FourthOrderMicromagneticSolver)
#
## GoogleTest requires at least C++11
#set(CMAKE_CXX_STANDARD 11)
#
#include(FetchContent)
#FetchContent_Declare(
#  googletest
#  URL https://github.com/google/googletest/archive/3e0e32ba300ce8afe695ad3ba7e81b21b7cf237.zip
#)
## For Windows: Prevent overriding the parent project's compiler/linker settings
#set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#FetchContent_MakeAvailable(googletest)
#
#enable_testing()
#
#add_executable(
#  VectorField_test
#  VectorField_test.cpp
#)
#target_link_libraries(
#  VectorField_test
#  gtest_main
#)
#
#include(GoogleTest)
#gtest_discover_tests(VectorField_test)