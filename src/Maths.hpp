#pragma once

#include "Constants.hpp"
#include "VectorField.hpp"

#include <cmath>

void crossProduct(float* a0, float* a1, float* a2, float* b0, float* b1, float* b2, float* out0, float* out1, float* out2)
{
    (*out0) = (*a1)*(*b2) - (*b1)*(*a2);
    (*out1) = (*a2)*(*b0) - (*a0)*(*b2);
    (*out2) = (*a0)*(*b1) - (*a1)*(*b0);
}


void normalise(VectorField *vec)
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            float magnitude = sqrt(vec->data[vec->getIndex(i, j, 0)] * vec->data[vec->getIndex(i, j, 0)] + vec->data[vec->getIndex(i, j, 1)] * vec->data[vec->getIndex(i, j, 1)] + vec->data[vec->getIndex(i, j, 2)] * vec->data[vec->getIndex(i, j, 2)]);
            vec->data[vec->getIndex(i, j, 0)] = vec->data[vec->getIndex(i, j, 0)] * (1/magnitude);
            vec->data[vec->getIndex(i, j, 1)] = vec->data[vec->getIndex(i, j, 1)] * (1/magnitude);
            vec->data[vec->getIndex(i, j, 2)] = vec->data[vec->getIndex(i, j, 2)] * (1/magnitude);
        }
    }
}
