#include <gtest/gtest.h>
#include "Maths.hpp"

#include <cmath>

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

    EXPECT_FLOAT_EQ(out0, 19.);
    EXPECT_FLOAT_EQ(out1, -25.);
    EXPECT_FLOAT_EQ(out2, -27.);
}

TEST(MathsTest, NormaliseTest)
{
    // Test normalisation of the vector field
    VectorField *vectorField = new VectorField(N, N, 3);
    int i = 50;
    int j = 71;
    int xIndex = vectorField->getIndex(i, j, 0);
    int yIndex = vectorField->getIndex(i, j, 1);
    int zIndex = vectorField->getIndex(i, j, 2);
    vectorField->data[xIndex] = 1;
    vectorField->data[yIndex] = 2;
    vectorField->data[zIndex] = -3;
    normalise(vectorField);
    EXPECT_FLOAT_EQ(vectorField->data[xIndex], 1/sqrt(14));
    EXPECT_FLOAT_EQ(vectorField->data[yIndex], 2/sqrt(14));
    EXPECT_FLOAT_EQ(vectorField->data[zIndex], -3/sqrt(14));
}