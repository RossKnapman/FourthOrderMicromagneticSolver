#include <gtest/gtest.h>
#include "VectorField.hpp"
#include "Constants.hpp"

TEST(VectorFieldTest, ConstructVector)
{
    // Check that the vector field initialises with the correct diemsions
    VectorField *vectorField = new VectorField(N, N, 3);
    EXPECT_EQ(vectorField->Nx, N);
    EXPECT_EQ(vectorField->Ny, N);
    EXPECT_EQ(vectorField->vecSize, 3);
    EXPECT_EQ(vectorField->data.size(), N*N*3);
}

TEST(VectorFieldTest, AssignAndRecall)
{
    // Test assigning a value to a given element of the vector field, then recalling it
    VectorField *vectorField = new VectorField(N, N, 3);
    vectorField->data[vectorField->getIndex(3, 4, 5)] = 5.;
    float outValue = vectorField->data[vectorField->getIndex(3, 4, 5)];
    EXPECT_FLOAT_EQ(outValue, 5.);
}

TEST(VectorFieldTest, CalculateIndex)
{
    // Test calculation of 1D index from 3D indices
    VectorField *vectorField = new VectorField(N, N, 2);
    int i = 10;
    int j = 20;
    int component = 1;
    int index = vectorField->getIndex(i, j, component);
    int expected = ((i*N + j) * 2) + component;
    EXPECT_EQ(index, expected);
}

TEST(VectorFieldTest, ElementWiseAddition)
{
    // Test elementwise addition of two vector fields
    VectorField *a = new VectorField(2, 2, 2);
    VectorField *b = new VectorField(2, 2, 2);
    a->data[4] = 3;
    b->data[4] = 4;
    VectorField *result = new VectorField(2, 2, 2);
    *result = *a + *b;
    EXPECT_FLOAT_EQ(result->data[4], 7);
}

TEST(VectorFieldTest, ElementWiseSubtraction)
{
    // Test elementwise subtraction of two vector fields
    VectorField *a = new VectorField(2, 2, 2);
    VectorField *b = new VectorField(2, 2, 2);
    a->data[4] = 3;
    b->data[4] = 4;
    VectorField *result = new VectorField(2, 2, 2);
    *result = *a - *b;
    EXPECT_FLOAT_EQ(result->data[4], -1);
}

TEST(VectorFieldTest, MultipleElementWiseAddition)
{
    // Test elementwise addition of three vector fields
    VectorField *a = new VectorField(2, 2, 2);
    VectorField *b = new VectorField(2, 2, 2);
    VectorField *c = new VectorField(2, 2, 2);
    a->data[4] = 3;
    b->data[4] = 4;
    c->data[4] = 10;
    VectorField *result = new VectorField(2, 2, 2);
    *result = *a + *b + *c;
    EXPECT_FLOAT_EQ(result->data[4], 17);
}

TEST(VectorFieldTest, MultipleElementwiseSubtraction)
{
    // Test elementwise addition of two vector fields, then further subtraction of another
    VectorField *a = new VectorField(2, 2, 2);
    VectorField *b = new VectorField(2, 2, 2);
    VectorField *c = new VectorField(2, 2, 2);
    a->data[4] = 3;
    b->data[4] = 4;
    c->data[4] = 10;
    VectorField *result = new VectorField(2, 2, 2);
    *result = *a + *b - *c;
    EXPECT_FLOAT_EQ(result->data[4], -3);
}

TEST(VectorFieldTest, MultiplyByScalar)
{
    // Test elementwise multiplication of a vector field by a scalar
    VectorField *vectorField = new VectorField(256, 256, 3);
    vectorField->data[400] = 5.;
    *vectorField = *vectorField * 5.;
    EXPECT_FLOAT_EQ(vectorField->data[400], 25.);
}

TEST(VectorFieldTest, AddMultipliedByScalar)
{
    // Test adding a vector field to another vector field that is multiplied by a scalar
    VectorField a(2, 2, 2);
    VectorField b(2, 2, 2);
    VectorField c(2, 2, 2);
    a.data[4] = 5;
    b.data[4] = 10;
    c.data[4] = 25;
    VectorField result(2, 2, 2);
    result = a + b*0.5 + c*5;
    EXPECT_FLOAT_EQ(result.data[4], 135);
}
