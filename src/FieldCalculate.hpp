#pragma once

#include "VectorField.hpp"
#include "Constants.hpp"
#include "Maths.hpp"


void calculateMagneticField(VectorField *m, VectorField *B, const float (&Ba)[3])
{
    for (int i=2; i<N-2; i++)
    {
        for (int j=2; j<N-2; j++)
        {
            // Reset to zero (adding the applied magnetic field first)
            B->data[B->getIndex(i, j, 0)] = Ba[0];
            B->data[B->getIndex(i, j, 1)] = Ba[1];
            B->data[B->getIndex(i, j, 2)] = Ba[2];

            // Central spin
            B->data[B->getIndex(i, j, 0)] -= (4 / (DELTA*DELTA)) * m->data[m->getIndex(i, j, 0)];
            B->data[B->getIndex(i, j, 1)] -= (4 / (DELTA*DELTA)) * m->data[m->getIndex(i, j, 1)];
            B->data[B->getIndex(i, j, 2)] -= (4 / (DELTA*DELTA)) * m->data[m->getIndex(i, j, 2)];
            B->data[B->getIndex(i, j, 0)] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i, j, 0)];
            B->data[B->getIndex(i, j, 1)] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i, j, 1)];
            B->data[B->getIndex(i, j, 2)] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i, j, 2)];

            // Cell directly to left
            B->data[B->getIndex(i, j, 0)] += m->data[m->getIndex(i-1, j, 0)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 1)] += m->data[m->getIndex(i-1, j, 1)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 2)] += m->data[m->getIndex(i-1, j, 2)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 0)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i-1, j, 0)];
            B->data[B->getIndex(i, j, 1)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i-1, j, 1)];
            B->data[B->getIndex(i, j, 2)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i-1, j, 2)];

            // Cell directly to right
            B->data[B->getIndex(i, j, 0)] += m->data[m->getIndex(i+1, j, 0)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 1)] += m->data[m->getIndex(i+1, j, 1)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 2)] += m->data[m->getIndex(i+1, j, 2)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 0)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i+1, j, 0)];
            B->data[B->getIndex(i, j, 1)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i+1, j, 1)];
            B->data[B->getIndex(i, j, 2)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i+1, j, 2)];

            // Cell directly below
            B->data[B->getIndex(i, j, 0)] += m->data[m->getIndex(i, j-1, 0)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 1)] += m->data[m->getIndex(i, j-1, 1)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 2)] += m->data[m->getIndex(i, j-1, 2)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 0)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i, j-1, 0)];
            B->data[B->getIndex(i, j, 1)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i, j-1, 1)];
            B->data[B->getIndex(i, j, 2)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i, j-1, 2)];

            // Cell directly above
            B->data[B->getIndex(i, j, 0)] += m->data[m->getIndex(i, j+1, 0)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 1)] += m->data[m->getIndex(i, j+1, 1)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 2)] += m->data[m->getIndex(i, j+1, 2)] / (DELTA*DELTA);
            B->data[B->getIndex(i, j, 0)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i, j+1, 0)];
            B->data[B->getIndex(i, j, 1)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i, j+1, 1)];
            B->data[B->getIndex(i, j, 2)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i, j+1, 2)];

            // Cell to bottom-left
            B->data[B->getIndex(i, j, 0)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i-1, j-1, 0)];
            B->data[B->getIndex(i, j, 1)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i-1, j-1, 1)];
            B->data[B->getIndex(i, j, 2)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i-1, j-1, 2)];

            // Cell to bottom-right
            B->data[B->getIndex(i, j, 0)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i+1, j-1, 0)];
            B->data[B->getIndex(i, j, 1)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i+1, j-1, 1)];
            B->data[B->getIndex(i, j, 2)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i+1, j-1, 2)];

            // Cell to top-left
            B->data[B->getIndex(i, j, 0)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i-1, j+1, 0)];
            B->data[B->getIndex(i, j, 1)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i-1, j+1, 1)];
            B->data[B->getIndex(i, j, 2)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i-1, j+1, 2)];

            // Cell to top-right
            B->data[B->getIndex(i, j, 0)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i+1, j+1, 0)];
            B->data[B->getIndex(i, j, 1)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i+1, j+1, 1)];
            B->data[B->getIndex(i, j, 2)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m->data[m->getIndex(i+1, j+1, 2)];

            // Cell two to left
            B->data[B->getIndex(i, j, 0)] -= m->data[m->getIndex(i-2, j, 0)] / (DELTA*DELTA*DELTA*DELTA);
            B->data[B->getIndex(i, j, 1)] -= m->data[m->getIndex(i-2, j, 1)] / (DELTA*DELTA*DELTA*DELTA);
            B->data[B->getIndex(i, j, 2)] -= m->data[m->getIndex(i-2, j, 2)] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two to right
            B->data[B->getIndex(i, j, 0)] -= m->data[m->getIndex(i+2, j, 0)] / (DELTA*DELTA*DELTA*DELTA);
            B->data[B->getIndex(i, j, 1)] -= m->data[m->getIndex(i+2, j, 1)] / (DELTA*DELTA*DELTA*DELTA);
            B->data[B->getIndex(i, j, 2)] -= m->data[m->getIndex(i+2, j, 2)] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two below
            B->data[B->getIndex(i, j, 0)] -= m->data[m->getIndex(i, j-2, 0)] / (DELTA*DELTA*DELTA*DELTA);
            B->data[B->getIndex(i, j, 1)] -= m->data[m->getIndex(i, j-2, 1)] / (DELTA*DELTA*DELTA*DELTA);
            B->data[B->getIndex(i, j, 2)] -= m->data[m->getIndex(i, j-2, 2)] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two above
            B->data[B->getIndex(i, j, 0)] -= m->data[m->getIndex(i, j+2, 0)] / (DELTA*DELTA*DELTA*DELTA);
            B->data[B->getIndex(i, j, 1)] -= m->data[m->getIndex(i, j+2, 1)] / (DELTA*DELTA*DELTA*DELTA);
            B->data[B->getIndex(i, j, 2)] -= m->data[m->getIndex(i, j+2, 2)] / (DELTA*DELTA*DELTA*DELTA);

        }
    }
}


void calculateTimeDerivative(VectorField *m, VectorField *B, VectorField *dmdt)
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            float mCrossB0;
            float mCrossB1;
            float mCrossB2;
            float mCrossmCrossB0;
            float mCrossmCrossB1;
            float mCrossmCrossB2;

            float m0 = m->data[m->getIndex(i, j, 0)];
            float m1 = m->data[m->getIndex(i, j, 1)];
            float m2 = m->data[m->getIndex(i, j, 2)];

            float B0 = B->data[B->getIndex(i, j, 0)];
            float B1 = B->data[B->getIndex(i, j, 1)];
            float B2 = B->data[B->getIndex(i, j, 2)];

            crossProduct(&m0, &m1, &m2, &B0, &B1, &B2, &mCrossB0, &mCrossB1, &mCrossB2);
            crossProduct(&m0, &m1, &m2, &mCrossB0, &mCrossB1, &mCrossB2, &mCrossmCrossB0, &mCrossmCrossB1, &mCrossmCrossB2);
            dmdt->data[dmdt->getIndex(i, j, 0)] = mCrossB0 + ALPHA * mCrossmCrossB0;
            dmdt->data[dmdt->getIndex(i, j, 1)] = mCrossB1 + ALPHA * mCrossmCrossB1;
            dmdt->data[dmdt->getIndex(i, j, 2)] = mCrossB2 + ALPHA * mCrossmCrossB2;
        }
    }
}
