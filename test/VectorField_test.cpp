#include <gtest/gtest.h>
#include "VectorField.hpp"

TEST(VectorFieldTest, AssignAndRecall)
{
    VectorField *vectorField = new VectorField(256, 256, 3);
    vectorField->setValue(3, 4, 5, 5.);
    float outValue = vectorField->getValue(3, 4, 5);
    EXPECT_EQ(outValue, 5.);
}