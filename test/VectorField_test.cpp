#include <gtest/gtest.h>
#include "VectorField.hpp"

TEST(VectorFieldTest, AssignAndRecall)
{
    VectorField *vectorField = new VectorField(256, 256, 3);
    vectorField[getIndex(3, 4, 5)] = 5.;
    float outValue = vectorField[getIndex(3, 4, 5)];
    EXPECT_EQ(outValue, 5.);
}