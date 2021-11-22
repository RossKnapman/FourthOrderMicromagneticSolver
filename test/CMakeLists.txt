set(BINARY ${CMAKE_PROJECT_NAME}_test)

set(TEST_SOURCES VectorField_test.cpp)

add_executable(${BINARY} ${TEST_SOURCES})

add_test(NAME ${BINARY} COMMAND ${BINARY})

# Use binary from src and add tests as target
target_link_libraries(${BINARY} PUBLIC {CMAKE_PROJECT_NAME}_lib gtest)