#include <gtest/gtest.h>
#include "Maths.hpp"

TEST(MathsTest, CrossProduct)
{
    float a0 = 1;
    float a1 = 4;
    float a2 = -3;
    float b0 = 7;
    float b1 = 1;
    float b2 = 4;
    float out0;
    float out1;
    float out2;

    crossProduct(&a0, &a1, &a2, &b0, &b1, &b2, &out0, &out1, &out2);

    EXPECT_EQ(out0, 18.);
    EXPECT_EQ(out1, -25.);
    EXPECT_EQ(out2, -27.);
}