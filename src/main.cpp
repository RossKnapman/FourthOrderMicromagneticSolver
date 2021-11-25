#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <cmath>

#include "VectorField.hpp"
#include "Maths.hpp"
#include "Constants.hpp"

using namespace std;

float t = 0;  // Total simulation time


void normalise(VectorField &vec)
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            float magnitude = sqrt(vec.data[vec.getIndex(i, j, 0)] * vec.data[vec.getIndex(i, j, 0)] + vec.data[vec.getIndex(i, j, 1)] * vec.data[vec.getIndex(i, j, 1)] + vec.data[vec.getIndex(i, j, 2)] * vec.data[vec.getIndex(i, j, 2)]);
            vec.data[vec.getIndex(i, j, 0)] = vec.data[vec.getIndex(i, j, 0)] * (1/magnitude);
            vec.data[vec.getIndex(i, j, 1)] = vec.data[vec.getIndex(i, j, 1)] * (1/magnitude);
            vec.data[vec.getIndex(i, j, 2)] = vec.data[vec.getIndex(i, j, 2)] * (1/magnitude);
        }
    }
}

void calculateMagneticField(VectorField &m, VectorField &B, const float (&Ba)[3])
{
    for (int i=2; i<N-2; i++)
    {
        for (int j=2; j<N-2; j++)
        {
            // Reset to zero (adding the applied magnetic field first)
            B.data[B.getIndex(i, j, 0)] = Ba[0];
            B.data[B.getIndex(i, j, 1)] = Ba[1];
            B.data[B.getIndex(i, j, 2)] = Ba[2];

            // Central spin
            B.data[B.getIndex(i, j, 0)] -= (4 / (DELTA*DELTA)) * m.data[m.getIndex(i, j, 0)];
            B.data[B.getIndex(i, j, 1)] -= (4 / (DELTA*DELTA)) * m.data[m.getIndex(i, j, 1)];
            B.data[B.getIndex(i, j, 2)] -= (4 / (DELTA*DELTA)) * m.data[m.getIndex(i, j, 2)];
            B.data[B.getIndex(i, j, 0)] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i, j, 0)];
            B.data[B.getIndex(i, j, 1)] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i, j, 1)];
            B.data[B.getIndex(i, j, 2)] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i, j, 2)];

            // Cell directly to left
            B.data[B.getIndex(i, j, 0)] += m.data[m.getIndex(i-1, j, 0)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 1)] += m.data[m.getIndex(i-1, j, 1)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 2)] += m.data[m.getIndex(i-1, j, 2)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 0)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i-1, j, 0)];
            B.data[B.getIndex(i, j, 1)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i-1, j, 1)];
            B.data[B.getIndex(i, j, 2)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i-1, j, 2)];

            // Cell directly to right
            B.data[B.getIndex(i, j, 0)] += m.data[m.getIndex(i+1, j, 0)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 1)] += m.data[m.getIndex(i+1, j, 1)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 2)] += m.data[m.getIndex(i+1, j, 2)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 0)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i+1, j, 0)];
            B.data[B.getIndex(i, j, 1)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i+1, j, 1)];
            B.data[B.getIndex(i, j, 2)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i+1, j, 2)];

            // Cell directly below
            B.data[B.getIndex(i, j, 0)] += m.data[m.getIndex(i, j-1, 0)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 1)] += m.data[m.getIndex(i, j-1, 1)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 2)] += m.data[m.getIndex(i, j-1, 2)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 0)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i, j-1, 0)];
            B.data[B.getIndex(i, j, 1)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i, j-1, 1)];
            B.data[B.getIndex(i, j, 2)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i, j-1, 2)];

            // Cell directly above
            B.data[B.getIndex(i, j, 0)] += m.data[m.getIndex(i, j+1, 0)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 1)] += m.data[m.getIndex(i, j+1, 1)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 2)] += m.data[m.getIndex(i, j+1, 2)] / (DELTA*DELTA);
            B.data[B.getIndex(i, j, 0)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i, j+1, 0)];
            B.data[B.getIndex(i, j, 1)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i, j+1, 1)];
            B.data[B.getIndex(i, j, 2)] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i, j+1, 2)];

            // Cell to bottom-left
            B.data[B.getIndex(i, j, 0)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i-1, j-1, 0)];
            B.data[B.getIndex(i, j, 1)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i-1, j-1, 1)];
            B.data[B.getIndex(i, j, 2)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i-1, j-1, 2)];

            // Cell to bottom-right
            B.data[B.getIndex(i, j, 0)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i+1, j-1, 0)];
            B.data[B.getIndex(i, j, 1)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i+1, j-1, 1)];
            B.data[B.getIndex(i, j, 2)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i+1, j-1, 2)];

            // Cell to top-left
            B.data[B.getIndex(i, j, 0)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i-1, j+1, 0)];
            B.data[B.getIndex(i, j, 1)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i-1, j+1, 1)];
            B.data[B.getIndex(i, j, 2)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i-1, j+1, 2)];

            // Cell to top-right
            B.data[B.getIndex(i, j, 0)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i+1, j+1, 0)];
            B.data[B.getIndex(i, j, 1)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i+1, j+1, 1)];
            B.data[B.getIndex(i, j, 2)] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[m.getIndex(i+1, j+1, 2)];

            // Cell two to left
            B.data[B.getIndex(i, j, 0)] -= m.data[m.getIndex(i-2, j, 0)] / (DELTA*DELTA*DELTA*DELTA);
            B.data[B.getIndex(i, j, 1)] -= m.data[m.getIndex(i-2, j, 1)] / (DELTA*DELTA*DELTA*DELTA);
            B.data[B.getIndex(i, j, 2)] -= m.data[m.getIndex(i-2, j, 2)] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two to right
            B.data[B.getIndex(i, j, 0)] -= m.data[m.getIndex(i+2, j, 0)] / (DELTA*DELTA*DELTA*DELTA);
            B.data[B.getIndex(i, j, 1)] -= m.data[m.getIndex(i+2, j, 1)] / (DELTA*DELTA*DELTA*DELTA);
            B.data[B.getIndex(i, j, 2)] -= m.data[m.getIndex(i+2, j, 2)] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two below
            B.data[B.getIndex(i, j, 0)] -= m.data[m.getIndex(i, j-2, 0)] / (DELTA*DELTA*DELTA*DELTA);
            B.data[B.getIndex(i, j, 1)] -= m.data[m.getIndex(i, j-2, 1)] / (DELTA*DELTA*DELTA*DELTA);
            B.data[B.getIndex(i, j, 2)] -= m.data[m.getIndex(i, j-2, 2)] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two above
            B.data[B.getIndex(i, j, 0)] -= m.data[m.getIndex(i, j+2, 0)] / (DELTA*DELTA*DELTA*DELTA);
            B.data[B.getIndex(i, j, 1)] -= m.data[m.getIndex(i, j+2, 1)] / (DELTA*DELTA*DELTA*DELTA);
            B.data[B.getIndex(i, j, 2)] -= m.data[m.getIndex(i, j+2, 2)] / (DELTA*DELTA*DELTA*DELTA);

        }
    }
}

void calculateTimeDerivative(const VectorField &m, const VectorField &B, VectorField &dmdt)
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

            float m0 = m.data[m.getIndex(i, j, 0)];
            float m1 = m.data[m.getIndex(i, j, 1)];
            float m2 = m.data[m.getIndex(i, j, 2)];

            float B0 = B.data[B.getIndex(i, j, 0)];
            float B1 = B.data[B.getIndex(i, j, 1)];
            float B2 = B.data[B.getIndex(i, j, 2)];

            crossProduct(&m0, &m1, &m2, &B0, &B1, &B2, &mCrossB0, &mCrossB1, &mCrossB2);
            crossProduct(&m0, &m1, &m2, &mCrossB0, &mCrossB1, &mCrossB2, &mCrossmCrossB0, &mCrossmCrossB1, &mCrossmCrossB2);
            dmdt.data[dmdt.getIndex(i, j, 0)] = mCrossB0 + ALPHA * mCrossmCrossB0;
            dmdt.data[dmdt.getIndex(i, j, 1)] = mCrossB1 + ALPHA * mCrossmCrossB1;
            dmdt.data[dmdt.getIndex(i, j, 2)] = mCrossB2 + ALPHA * mCrossmCrossB2;
        }
    }
}

void step(VectorField &m, VectorField &B)
{

    // Dorman-Prince (RK45) integration
    VectorField k1(N, N, 3);
    VectorField k2(N, N, 3);
    VectorField k3(N, N, 3);
    VectorField k4(N, N, 3);
    VectorField k5(N, N, 3);
    VectorField k6(N, N, 3);


    // Calculate k1
    calculateTimeDerivative(m, B, k1);
    k1 = k1.multiplyByScalar(h);

    // Calculate k2
    calculateTimeDerivative(m.add(k1.multiplyByScalar(1./5.)), B, k2);
    k2 = k2.multiplyByScalar(h);

    // Calculate k3
    calculateTimeDerivative(m.add(k1.multiplyByScalar(3./40.).add(k2.multiplyByScalar(9./40.))), B, k3);
    k3 = k3.multiplyByScalar(h);

    // Calculate k4
    calculateTimeDerivative(m.add(k1.multiplyByScalar(44./45.).subtract(k2.multiplyByScalar(56./15.).add(k3.multiplyByScalar(32./9)))), B, k4);
    k4 = k4.multiplyByScalar(h);

    // Calculate k5
    calculateTimeDerivative(m.add(k1.multiplyByScalar(19372./6561.).subtract(k2.multiplyByScalar(25360./2187.).add(k3.multiplyByScalar(64448./6561.).subtract(k4.multiplyByScalar(212./729.))))), B, k5);
    k5 = k5.multiplyByScalar(h);

    // Calculate k6
    calculateTimeDerivative(m.add(k1.multiplyByScalar(9017./3168.).subtract(k2.multiplyByScalar(355./33.).subtract(k3.multiplyByScalar(46732./5247.).add(k4.multiplyByScalar(49./176.).subtract(k5.multiplyByScalar(5103./18656.)))))), B, k6);
    k6 = k6.multiplyByScalar(h);

    m = m.add(k1.multiplyByScalar(35./384.).add(k3.multiplyByScalar(500./1113.).add(k4.multiplyByScalar(125./192.).subtract(k5.multiplyByScalar(2187./6784.).add(k6.multiplyByScalar(11./64.))))));
    normalise(m);
}

void writeFile(string name, VectorField &m)
{
    ofstream file;
    file.open(name);

    // Write out the header
    file << "# OOMMF OVF 2.0" << endl;
    file << "# Segment count: 1" << endl;
    file << "# Begin: Segment" << endl;
    file << "# Begin: Header" << endl;
    file << "# Title: m" << endl;
    file << "# meshtype: rectangular" << endl;
    file << "# meshunit: dimensionless" << endl;
    file << "# valueunit: dimensionless" << endl;
    file << "# xmin: 0" << endl;
    file << "# ymin: 0" << endl;
    file << "# zmin: 0" << endl;
    file << "# xmax: " << to_string(DELTA * N) << endl;
    file << "# ymax: " << to_string(DELTA * N) << endl;
    file << "# zmax: " << to_string(DELTA) << endl;
    file << "# valuedim: 3" << endl;
    file << "# valuelabels: m_x m_y m_z" << endl;
    file << "# Desc: Total simulation time: " << to_string(t) << endl;
    file << "# xbase: " << to_string(DELTA/2) << endl;
    file << "# ybase: " << to_string(DELTA/2) << endl;
    file << "# zbase: " << to_string(DELTA/2) << endl;
    file << "# xnodes: " << to_string(N) << endl;
    file << "# ynodes: " << to_string(N) << endl;
    file << "# znodes: 1" << endl;
    file << "# xstepsize: " << to_string(DELTA) << endl;
    file << "# ystepsize: " << to_string(DELTA) << endl;
    file << "# zstepsize: " << to_string(DELTA) << endl;
    file << "# End: Header" << endl;
    file << "# Begin: Data Text" << endl;

    // OVF uses the Fortran convention, i.e. x changes fastest
    for (int j=0; j<N; j++)
    {
        for (int i=0; i<N; i++)
        {
            file << m.data[m.getIndex(i, j, 0)] << " " << m.data[m.getIndex(i, j, 1)] << " " << m.data[m.getIndex(i, j, 2)] << endl;
        }
    }

    file << "# End: Data Text" << endl;
    file << "# End: Segment" << endl;

    file.close();
}

int main()
{
    VectorField m(N, N, 3);  // Declare magnetization array
    VectorField B(N, N, 3);  // Declare emergent field array
    VectorField pos(N, N, 2);  // Declare position array

    const float Ba[3] = {0., 0., 1.};  // Applied magnetic field

    // initialisePosition(pos);
    // initialiseMagnetization(m, pos, -1, 1, 3);

    // calculateMagneticField(m, B, Ba);

    // VectorField dmdt(N, N, 3);
    // calculateTimeDerivative(m, B, dmdt);

    // int counter = 0;
    // while (counter < 100)
    // {
    //     cout << counter << endl;
    //     calculateMagneticField(m, B, Ba);
    //     step(m, B);
    //     char outName[16];
    //     sprintf(outName, "data/m%06d.ovf", counter);
    //     writeFile(outName, m);
    //     t += h;
    //     counter++;
    // }

    return 0;
}
