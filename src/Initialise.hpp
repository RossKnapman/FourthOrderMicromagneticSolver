#include "VectorField.hpp"
#include "Constants.hpp"

#include <cmath>
#include <iostream>
#include <string>

void initialisePosition(VectorField *pos)
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            pos->data[pos->getIndex(i, j, 0)] = i * DELTA;
            pos->data[pos->getIndex(i, j, 1)] = j * DELTA;
        }
    }
}


void initialiseMagnetizationSkyrmion(VectorField *m, VectorField *pos, float pol, float charge, float w)
{
    // Initialise with a skyrmion in the middle of the grid
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            float xMax = pos->data[pos->getIndex(N-1, N-1, 0)];

            float x = pos->data[pos->getIndex(i, j, 0)] - xMax/2;
            float y = pos->data[pos->getIndex(i, j, 1)] - xMax/2;
            float r2 = x*x + y*y;
            float r = sqrt(r2);
            float w2 = w*w;

            float mz = 2 * pol * (exp(-r2 / w2) - 0.5);
            float mx = (x * charge / r) * (1 - abs(mz));
            float my = (y * charge / r) * (1 - abs(mz));
            
            m->data[m->getIndex(i, j, 0)] = mx;
            m->data[m->getIndex(i, j, 1)] = my;
            m->data[m->getIndex(i, j, 2)] = mz;
        }
    }
}
