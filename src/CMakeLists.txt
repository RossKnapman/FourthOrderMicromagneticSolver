set (BINARY ${CMAKE_PROJECT_NAME})

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.hpp *.cpp)

add_executable(${BINARY} ${SOURCES})

# For including the code in unit tests
add_library(${BINARY}_lib STATIC ${SOURCES})