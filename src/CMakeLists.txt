set (BINARY ${CMAKE_PROJECT_NAME})

set(SOURCES main.cpp)

add_executable(${BINARY} ${SOURCES})

# For including the code in unit tests
add_library(${BINARY}_lib STATIC ${SOURCES})