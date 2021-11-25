#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include "Initialise.hpp"
#include "Constants.hpp"

TEST(InitialiseTest, InitialisePosition)
{
    // Test initialising the position array
    VectorField *pos = new VectorField(N, N, 2);
    initialisePosition(pos);
    int i = 97;
    int j = 128;
    int idx = pos->getIndex(i, j, 1);
    float expected = j * DELTA;
    EXPECT_EQ(pos->data[idx], expected);
}
