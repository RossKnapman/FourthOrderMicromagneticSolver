#include <gtest/gtest.h>
#include "VectorField.hpp"

TEST(VectorFieldTest, AssignAndRecall)
{
    VectorField *vectorField = new VectorField(256, 256, 3);
    vectorField->data[vectorField->getIndex(3, 4, 5)] = 5.;
    float outValue = vectorField->data[vectorField->getIndex(3, 4, 5)];
    EXPECT_EQ(outValue, 5.);
}

TEST(VectorFieldTest, multiplyByScalar)
{
    VectorField *vectorField = new VectorField(256, 256, 3);
    vectorField->data[400] = 5.;
    *vectorField = vectorField->multiplyByScalar(5.);
    EXPECT_EQ(vectorField->data[400], 25.);
}