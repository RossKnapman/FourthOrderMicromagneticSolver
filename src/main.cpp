#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <cmath>

#include "VectorField.hpp"

#define N 512
#define DELTA 0.2
#define ALPHA 1e-3
#define h 1e-2
#define T_FIN 0.1

using namespace std;

float t = 0;  // Total simulation time

void crossProduct(float* a, float* b, float* out)
{
    out[0] = a[1] * b[2] - b[1] * a[2];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

void normalise(float (&a)[3], float (&out)[3])
{
    float magnitude = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    out[0] = a[0] / magnitude;
    out[1] = a[1] / magnitude;
    out[2] = a[2] / magnitude;
}

void normalise(VectorField<float, N, N, 2> &vec)
{
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            float magnitude = sqrt(vec.data[i][j][0] * vec.data[i][j][0] + vec.data[i][j][1] * vec.data[i][j][1] + vec.data[i][j][2] * vec.data[i][j][2];
            vec[i][j][0] = vec[i][j][0] * (1/magnitude);
            vec[i][j][1] = vec[i][j][1] * (1/magnitude);
            vec[i][j][2] = vec[i][j][2] * (1/magnitude);
}

void initialisePosition(VectorField<float, N, N, 2> &pos)
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            pos.data[i][j][0] = i * DELTA;
            pos.data[i][j][1] = j * DELTA;
        }
    }
}

void initialiseMagnetization(VectorField<float, N, N, 3> &m, VectorField<float, N, N, 2> &pos, float pol, float charge, float w)
{
    // Initialise with a skyrmion in the middle of the grid
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            float xMax = pos.data[N-1][N-1][0];

            float x = pos.data[i][j][0] - xMax/2;
            float y = pos.data[i][j][1] - xMax/2;
            float r2 = x*x + y*y;
            float r = sqrt(r2);
            float w2 = w*w;

            float mz = 2 * pol * (exp(-r2 / w2) - 0.5);
            float mx = (x * charge / r) * (1 - abs(mz));
            float my = (y * charge / r) * (1 - abs(mz));
            
            m.data[i][j][0] = mx;
            m.data[i][j][1] = my;
            m.data[i][j][2] = mz;
        }
    }
}

void calculateMagneticField(VectorField<float, N, N, 3> &m, VectorField<float, N, N, 3> &B, const float (&Ba)[3])
{
    for (int i=2; i<N-2; i++)
    {
        for (int j=2; j<N-2; j++)
        {
            // Reset to zero (adding the applied magnetic field first)
            B.data[i][j][0] = Ba[0];
            B.data[i][j][1] = Ba[1];
            B.data[i][j][2] = Ba[2];

            // Central spin
            B.data[i][j][0] -= (4 / (DELTA*DELTA)) * m.data[i][j][0];
            B.data[i][j][1] -= (4 / (DELTA*DELTA)) * m.data[i][j][1];
            B.data[i][j][2] -= (4 / (DELTA*DELTA)) * m.data[i][j][2];
            B.data[i][j][0] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i][j][0];
            B.data[i][j][1] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i][j][1];
            B.data[i][j][2] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i][j][2];

            // Cell directly to left
            B.data[i][j][0] += m.data[i-1][j][0] / (DELTA*DELTA);
            B.data[i][j][1] += m.data[i-1][j][1] / (DELTA*DELTA);
            B.data[i][j][2] += m.data[i-1][j][2] / (DELTA*DELTA);
            B.data[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i-1][j][0];
            B.data[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i-1][j][1];
            B.data[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i-1][j][2];

            // Cell directly to right
            B.data[i][j][0] += m.data[i+1][j][0] / (DELTA*DELTA);
            B.data[i][j][1] += m.data[i+1][j][1] / (DELTA*DELTA);
            B.data[i][j][2] += m.data[i+1][j][2] / (DELTA*DELTA);
            B.data[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i+1][j][0];
            B.data[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i+1][j][1];
            B.data[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i+1][j][2];

            // Cell directly below
            B.data[i][j][0] += m.data[i][j-1][0] / (DELTA*DELTA);
            B.data[i][j][1] += m.data[i][j-1][1] / (DELTA*DELTA);
            B.data[i][j][2] += m.data[i][j-1][2] / (DELTA*DELTA);
            B.data[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i][j-1][0];
            B.data[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i][j-1][1];
            B.data[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i][j-1][2];

            // Cell directly above
            B.data[i][j][0] += m.data[i][j+1][0] / (DELTA*DELTA);
            B.data[i][j][1] += m.data[i][j+1][1] / (DELTA*DELTA);
            B.data[i][j][2] += m.data[i][j+1][2] / (DELTA*DELTA);
            B.data[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i][j+1][0];
            B.data[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i][j+1][1];
            B.data[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i][j+1][2];

            // Cell to bottom-left
            B.data[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i-1][j-1][0];
            B.data[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i-1][j-1][1];
            B.data[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i-1][j-1][2];

            // Cell to bottom-right
            B.data[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i+1][j-1][0];
            B.data[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i+1][j-1][1];
            B.data[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i+1][j-1][2];

            // Cell to top-left
            B.data[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i-1][j+1][0];
            B.data[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i-1][j+1][1];
            B.data[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i-1][j+1][2];

            // Cell to top-right
            B.data[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i+1][j+1][0];
            B.data[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i+1][j+1][1];
            B.data[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m.data[i+1][j+1][2];

            // Cell two to left
            B.data[i][j][0] -= m.data[i-2][j][0] / (DELTA*DELTA*DELTA*DELTA);
            B.data[i][j][1] -= m.data[i-2][j][1] / (DELTA*DELTA*DELTA*DELTA);
            B.data[i][j][2] -= m.data[i-2][j][2] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two to right
            B.data[i][j][0] -= m.data[i+2][j][0] / (DELTA*DELTA*DELTA*DELTA);
            B.data[i][j][1] -= m.data[i+2][j][1] / (DELTA*DELTA*DELTA*DELTA);
            B.data[i][j][2] -= m.data[i+2][j][2] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two below
            B.data[i][j][0] -= m.data[i][j-2][0] / (DELTA*DELTA*DELTA*DELTA);
            B.data[i][j][1] -= m.data[i][j-2][1] / (DELTA*DELTA*DELTA*DELTA);
            B.data[i][j][2] -= m.data[i][j-2][2] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two above
            B.data[i][j][0] -= m.data[i][j+2][0] / (DELTA*DELTA*DELTA*DELTA);
            B.data[i][j][1] -= m.data[i][j+2][1] / (DELTA*DELTA*DELTA*DELTA);
            B.data[i][j][2] -= m.data[i][j+2][2] / (DELTA*DELTA*DELTA*DELTA);

        }
    }
}

void calculateTimeDerivative(VectorField<float, N, N, 3> &m, VectorField<float, N, N, 3> &B, VectorField<float, N, N, 3> &dmdt)
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            float mCrossB[3];
            float mCrossmCrossB[3];
            crossProduct(m.data[i][j], B.data[i][j], mCrossB);
            crossProduct(m.data[i][j], mCrossB, mCrossmCrossB);
            dmdt.data[i][j][0] = mCrossB[0] + ALPHA * mCrossmCrossB[0];
            dmdt.data[i][j][1] = mCrossB[1] + ALPHA * mCrossmCrossB[1];
            dmdt.data[i][j][2] = mCrossB[2] + ALPHA * mCrossmCrossB[2];
        }
    }
}

void step(VectorField<float, N, N, 3> &m, VectorField<float, N, N, 3> &B)
{

    // Dorman-Prince (RK45) integration
    VectorField<float, N, N, 3> k1;
    VectorField<float, N, N, 3> k2;
    VectorField<float, N, N, 3> k3;
    VectorField<float, N, N, 3> k4;
    VectorField<float, N, N, 3> k5;
    VectorField<float, N, N, 3> k6;


    // Calculate k1
    calculateTimeDerivative(m, B, k1);
    k1 = k1 * h;

    // Calculate k2
    calculateTimeDerivative(m + k1*(1./5.), B, k2);
    k2 = k2 * h;

    // Calculate k3
    calculateTimeDerivative(m + k1*(3./40.) + k2*(9./40.), B, k3);
    k3 = k3 * h;

    // Calculate k4
    calculateTimeDerivative(m + k1*(44./45.) - k2*(56./15.) + k3*(32./9/), B, k4)
    k4 = k4 * h;

    // Calculate k5
    calculateTimeDerivative(m + k1*(19372./6561.) - k2*(25360./2187.) + k3*(64448./6561.) - k4*(212./729.), B, k5);
    k5 = k5 * h;

    // Calculate k6
    calculateTimeDerivative(m + k1*(9017./3168.) - k2*(355./33.) - k3*(46732./5247.) + k4*(49./176.) - k5*(5103./18656.), B, k6);
    k6 = k6 * h;

    m = m + k1*(35./384.) + k3*(500./1113.) + k4*(125./192.) - k5*(2187./6784.) + k6*(11./64.);
    normalise(m);
    }
}

void writeFile(string name, VectorField<float, N, N, 3> &m)
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
            file << m.data[i][j][0] << " " << m.data[i][j][1] << " " << m.data[i][j][2] << endl;
        }
    }

    file << "# End: Data Text" << endl;
    file << "# End: Segment" << endl;

    file.close();
}

int main()
{
    VectorField<float, N, N, 3> m;  // Declare magnetization array
    VectorField<float, N, N, 3> B;  // Declare emergent field array
    VectorField<float, N, N, 2> pos;  // Declare position array

    const float Ba[3] = {0., 0., 1.};  // Applied magnetic field

    initialisePosition(pos);
    initialiseMagnetization(m, pos, -1, 1, 3);

    calculateMagneticField(m, B, Ba);

    VectorField<float, N, N, 3> dmdt;
    calculateTimeDerivative(m, B, dmdt);
    writeFile("data/m000000.ovf", dmdt);

    VectorField<float, 1, 1, 3> test;
    test.data[0][0][0] = 10;
    test.data[0][0][1] = 21;
    test.data[0][0][2] = -100;
    normalise(test);

    std::cout << to_string(test.data[0][0][0]) << " " << to_string(test.data[0][0][1]) << " " << to_string(test.data[0][0][2]) << endl;


    // int counter = 0;
    // while (counter < 100)
    // {
    //     cout << counter << endl;
    //     calculateMagneticField(m, B, Ba);
    //     step(m, B);
    //     char outName[16];
    //     sprintf(outName, "data/m%06d.ovf", counter);
    //     writeFile(outName, r, m);
    //     t += h;
    //     counter++;
    // }

    return 0;
}