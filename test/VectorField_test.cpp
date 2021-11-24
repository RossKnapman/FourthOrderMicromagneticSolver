#include <gtest/gtest.h>
#include "VectorField.hpp"

TEST(VectorFieldTest, AssignAndRecall)
{
    VectorField *vectorField = new VectorField(256, 256, 3);
    vectorField->data[vectorField->getIndex(3, 4, 5)] = 5.;
    float outValue = vectorField->data[vectorField->getIndex(3, 4, 5)];
    EXPECT_EQ(outValue, 5.);
}

TEST(VectorFieldTest, ElementWiseAddition)
{
    VectorField *a = new VectorField(2, 2, 2);
    VectorField *b = new VectorField(2, 2, 2);
    a->data[4] = 3;
    b->data[4] = 4;
    VectorField *result = new VectorField(2, 2, 2);
    *result = a->add(*b);
    EXPECT_EQ(result->data[4], 7);
}

TEST(VectorFieldTest, MultiplyByScalar)
{
    VectorField *vectorField = new VectorField(256, 256, 3);
    vectorField->data[400] = 5.;
    *vectorField = vectorField->multiplyByScalar(5.);
    EXPECT_EQ(vectorField->data[400], 25.);
}