#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <cmath>

#include "vectorfield.hpp"

#define N 512
#define DELTA 0.2
#define ALPHA 1e-3
#define h 1e-2
#define T_FIN 0.1

using namespace std;

float t = 0;  // Total simulation time

void crossProduct(float (&a)[3], float (&b)[3], float (&out)[3])
{
    out[0] = a[1] * b[2] - b[1] * a[2];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

void arrayAdd2(float (&a)[N][N][3], float(&b)[N][N][3], float (&out)[N][N][3])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            out[i][j][0] = a[i][j][0] + b[i][j][0];
            out[i][j][1] = a[i][j][1] + b[i][j][1];
            out[i][j][2] = a[i][j][2] + b[i][j][2];
        }
    }
}

void arrayAdd3(float (&a)[N][N][3], float(&b)[N][N][3], float (&c)[N][N][3], float (&out)[N][N][3])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            out[i][j][0] = a[i][j][0] + b[i][j][0] + c[i][j][0];
            out[i][j][1] = a[i][j][1] + b[i][j][1] + c[i][j][1];
            out[i][j][2] = a[i][j][2] + b[i][j][2] + c[i][j][2];
        }
    }
}

void arrayAdd4(float (&a)[N][N][3], float(&b)[N][N][3], float (&c)[N][N][3], float (&d)[N][N][3], float (&out)[N][N][3])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            out[i][j][0] = a[i][j][0] + b[i][j][0] + c[i][j][0] + d[i][j][0];
            out[i][j][1] = a[i][j][1] + b[i][j][1] + c[i][j][1] + d[i][j][1];
            out[i][j][2] = a[i][j][2] + b[i][j][2] + c[i][j][2] + d[i][j][2];
        }
    }
}

void arrayAdd5(float (&a)[N][N][3], float(&b)[N][N][3], float (&c)[N][N][3], float (&d)[N][N][3], float (&e)[N][N][3], float (&out)[N][N][3])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            out[i][j][0] = a[i][j][0] + b[i][j][0] + c[i][j][0] + d[i][j][0] + e[i][j][0];
            out[i][j][1] = a[i][j][1] + b[i][j][1] + c[i][j][1] + d[i][j][1] + e[i][j][1];
            out[i][j][2] = a[i][j][2] + b[i][j][2] + c[i][j][2] + d[i][j][2] + e[i][j][2];
        }
    }
}

void arrayAdd6(float (&a)[N][N][3], float(&b)[N][N][3], float (&c)[N][N][3], float (&d)[N][N][3], float (&e)[N][N][3], float (&f)[N][N][3], float (&out)[N][N][3])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            out[i][j][0] = a[i][j][0] + b[i][j][0] + c[i][j][0] + d[i][j][0] + e[i][j][0] + f[i][j][0];
            out[i][j][1] = a[i][j][1] + b[i][j][1] + c[i][j][1] + d[i][j][1] + e[i][j][1] + f[i][j][1];
            out[i][j][2] = a[i][j][2] + b[i][j][2] + c[i][j][2] + d[i][j][2] + e[i][j][2] + f[i][j][2];
        }
    }
}

void normalise(float (&a)[3], float (&out)[3])
{
    float magnitude = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    out[0] = a[0] / magnitude;
    out[1] = a[1] / magnitude;
    out[2] = a[2] / magnitude;
}

void initialisePosition(float (&r)[N][N][2])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            r[i][j][0] = i * DELTA;
            r[i][j][1] = j * DELTA;
        }
    }
}

void initialiseMagnetization(float (&m)[N][N][3], float (&r)[N][N][2], float pol, float charge, float w)
{
    // Initialise with a skyrmion in the middle of the grid
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            float xMax = r[N-1][N-1][0];

            float x = r[i][j][0] - xMax/2;
            float y = r[i][j][1] - xMax/2;
            float r2 = x*x + y*y;
            float r = sqrt(r2);
            float w2 = w*w;

            float mz = 2 * pol * (exp(-r2 / w2) - 0.5);
            float mx = (x * charge / r) * (1 - abs(mz));
            float my = (y * charge / r) * (1 - abs(mz));
            
            m[i][j][0] = mx;
            m[i][j][1] = my;
            m[i][j][2] = mz;
        }
    }
}

void calculateMagneticField(float *m, float *B, const float (&Ba)[3])
{
    for (int i=2; i<N-2; i++)
    {
        for (int j=2; j<N-2; j++)
        {
            // Reset to zero (adding the applied magnetic field first)
            B[i][j][0] = Ba[0];
            B[i][j][1] = Ba[1];
            B[i][j][2] = Ba[2];

            // Central spin
            B[i][j][0] -= (4 / (DELTA*DELTA)) * m[i][j][0];
            B[i][j][1] -= (4 / (DELTA*DELTA)) * m[i][j][1];
            B[i][j][2] -= (4 / (DELTA*DELTA)) * m[i][j][2];
            // B[i][j][0] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j][0];
            // B[i][j][1] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j][1];
            // B[i][j][2] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j][2];

            // Cell directly to left
            B[i][j][0] += m[i-1][j][0] / (DELTA*DELTA);
            B[i][j][1] += m[i-1][j][1] / (DELTA*DELTA);
            B[i][j][2] += m[i-1][j][2] / (DELTA*DELTA);
            // B[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j][0];
            // B[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j][1];
            // B[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j][2];

            // Cell directly to right
            B[i][j][0] += m[i+1][j][0] / (DELTA*DELTA);
            B[i][j][1] += m[i+1][j][1] / (DELTA*DELTA);
            B[i][j][2] += m[i+1][j][2] / (DELTA*DELTA);
            // B[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j][0];
            // B[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j][1];
            // B[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j][2];

            // Cell directly below
            B[i][j][0] += m[i][j-1][0] / (DELTA*DELTA);
            B[i][j][1] += m[i][j-1][1] / (DELTA*DELTA);
            B[i][j][2] += m[i][j-1][2] / (DELTA*DELTA);
            // B[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j-1][0];
            // B[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j-1][1];
            // B[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j-1][2];

            // Cell directly above
            B[i][j][0] += m[i][j+1][0] / (DELTA*DELTA);
            B[i][j][1] += m[i][j+1][1] / (DELTA*DELTA);
            B[i][j][2] += m[i][j+1][2] / (DELTA*DELTA);
            // B[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j+1][0];
            // B[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j+1][1];
            // B[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j+1][2];

            // // Cell to bottom-left
            // B[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j-1][0];
            // B[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j-1][1];
            // B[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j-1][2];

            // // Cell to bottom-right
            // B[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j-1][0];
            // B[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j-1][1];
            // B[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j-1][2];

            // // Cell to top-left
            // B[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j+1][0];
            // B[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j+1][1];
            // B[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j+1][2];

            // // Cell to top-right
            // B[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j+1][0];
            // B[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j+1][1];
            // B[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j+1][2];

            // // Cell two to left
            // B[i][j][0] -= m[i-2][j][0] / (DELTA*DELTA*DELTA*DELTA);
            // B[i][j][1] -= m[i-2][j][1] / (DELTA*DELTA*DELTA*DELTA);
            // B[i][j][2] -= m[i-2][j][2] / (DELTA*DELTA*DELTA*DELTA);

            // // Cell two to right
            // B[i][j][0] -= m[i+2][j][0] / (DELTA*DELTA*DELTA*DELTA);
            // B[i][j][1] -= m[i+2][j][1] / (DELTA*DELTA*DELTA*DELTA);
            // B[i][j][2] -= m[i+2][j][2] / (DELTA*DELTA*DELTA*DELTA);

            // // Cell two below
            // B[i][j][0] -= m[i][j-2][0] / (DELTA*DELTA*DELTA*DELTA);
            // B[i][j][1] -= m[i][j-2][1] / (DELTA*DELTA*DELTA*DELTA);
            // B[i][j][2] -= m[i][j-2][2] / (DELTA*DELTA*DELTA*DELTA);

            // // Cell two above
            // B[i][j][0] -= m[i][j+2][0] / (DELTA*DELTA*DELTA*DELTA);
            // B[i][j][1] -= m[i][j+2][1] / (DELTA*DELTA*DELTA*DELTA);
            // B[i][j][2] -= m[i][j+2][2] / (DELTA*DELTA*DELTA*DELTA);

        }
    }
}

void calculateTimeDerivative(float (&m)[N][N][3], float (&B)[N][N][3], float(&dmdt)[N][N][3])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            float mCrossB[3];
            float mCrossmCrossB[3];
            crossProduct(m[i][j], B[i][j], mCrossB);
            crossProduct(m[i][j], mCrossB, mCrossmCrossB);
            dmdt[i][j][0] = mCrossB[0] + ALPHA * mCrossmCrossB[0];
            dmdt[i][j][1] = mCrossB[1] + ALPHA * mCrossmCrossB[1];
            dmdt[i][j][2] = mCrossB[2] + ALPHA * mCrossmCrossB[2];
        }
    }
}

void step(float (&m)[N][N][3], float (&B)[N][N][3])
{

    // Dorman-Prince (RK45) integration
    float k1[N][N][3];
    float k2[N][N][3];
    float k3[N][N][3];
    float k4[N][N][3];
    float k5[N][N][3];
    float k6[N][N][3];

    // Calculate k1
    calculateTimeDerivative(m, B, k1);
    scalarMultiply(h, k1, k1);

    // Calculate k2
    float k1fork2[N][N][3];
    float k2Field[N][N][3];
    scalarMultiply((1./5.), k1, k1fork2);
    arrayAdd2(m, k1fork2, k2Field);
    calculateTimeDerivative(k2Field, B, k2);
    scalarMultiply(h, k2, k2);

    // Calculate k3
    float k1fork3[N][N][3];
    float k2fork3[N][N][3];
    float k3Field[N][N][3];
    scalarMultiply((3./40.), k1, k1fork3);
    scalarMultiply((9./40.), k2, k2fork3);
    arrayAdd3(m, k1fork3, k2fork3, k3Field);
    calculateTimeDerivative(k3Field, B, k3);
    scalarMultiply(h, k3, k3);

    // Calculate k4
    float k1fork4[N][N][3];
    float k2fork4[N][N][3];
    float k3fork4[N][N][3];
    float k4Field[N][N][3];
    scalarMultiply((44./45.), k1, k1fork4);
    scalarMultiply(-(56./15.), k2, k2fork4);
    scalarMultiply((32./9.), k3, k3fork4);
    arrayAdd4(m, k1fork4, k2fork4, k3fork4, k4Field);
    calculateTimeDerivative(k4Field, B, k4);
    scalarMultiply(h, k4, k4);

    // Calculate k5
    float k1fork5[N][N][3];
    float k2fork5[N][N][3];
    float k3fork5[N][N][3];
    float k4fork5[N][N][3];
    float k5Field[N][N][3];
    scalarMultiply((19372./6561.), k1, k1fork5);
    scalarMultiply(-(25360./2187.), k2, k2fork5);
    scalarMultiply((64448./6561), k3, k3fork5);
    scalarMultiply(-(212./729.), k4, k4fork5);
    arrayAdd5(m, k1fork5, k1fork5, k3fork5, k4fork5, k5Field);
    calculateTimeDerivative(k5Field, B, k5);
    scalarMultiply(h, k5, k5);

    // Calculate k6
    float k1fork6[N][N][3];
    float k2fork6[N][N][3];
    float k3fork6[N][N][3];
    float k4fork6[N][N][3];
    float k5fork6[N][N][3];
    float k6Field[N][N][3];
    scalarMultiply((9017./3168.), k1, k1fork6);
    scalarMultiply(-(355./33.), k2, k2fork6);
    scalarMultiply(-(46732./5247.), k3, k3fork6);
    scalarMultiply((49./176.), k4, k4fork6);
    scalarMultiply(-(5103./18656.), k5, k5fork6);
    arrayAdd6(m, k1fork6, k2fork6, k3fork6, k4fork6, k5fork6, k6Field);
    calculateTimeDerivative(k6Field, B, k6);
    scalarMultiply(h, k6, k6);

    float k1Final[N][N][3];
    float k3Final[N][N][3];
    float k4Final[N][N][3];
    float k5Final[N][N][3];
    float k6Final[N][N][3];
    scalarMultiply((35./384.), k1, k1Final);
    scalarMultiply((500./1113.), k3, k3Final);
    scalarMultiply((125./192.), k4, k4Final);
    scalarMultiply(-(2187./6784.), k5, k5Final);
    scalarMultiply((11./84.), k6, k6Final);

    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            m[i][j][0] += k1Final[i][j][0];
            m[i][j][1] += k1Final[i][j][1];
            m[i][j][2] += k1Final[i][j][2];

            m[i][j][0] += k3Final[i][j][0];
            m[i][j][1] += k3Final[i][j][1];
            m[i][j][2] += k3Final[i][j][2];

            m[i][j][0] += k4Final[i][j][0];
            m[i][j][1] += k4Final[i][j][1];
            m[i][j][2] += k4Final[i][j][2];

            m[i][j][0] += k5Final[i][j][0];
            m[i][j][1] += k5Final[i][j][1];
            m[i][j][2] += k5Final[i][j][2];

            m[i][j][0] += k6Final[i][j][0];
            m[i][j][1] += k6Final[i][j][1];
            m[i][j][2] += k6Final[i][j][2];

            normalise(m[i][j], m[i][j]);
        }
    }
}

void writeFile(string name, float (&r)[N][N][2], float (&m)[N][N][3])
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
            file << m[i][j][0] << " " << m[i][j][1] << " " << m[i][j][2] << endl;
        }
    }

    file << "# End: Data Text" << endl;
    file << "# End: Segment" << endl;

    file.close();
}

int main()
{

    auto B = new float [N][N][3];  // Declare emergent field array
    auto r = new float [N][N][3];  // Declare position array
    auto m = new float [N][N][3];  // Declare magnetization array

    const float Ba[3] = {0., 0., 1.};  // Applied magnetic field

    initialisePosition(r);
    initialiseMagnetization(m, r, -1, 1, 3);
    calculateMagneticField(m, B, Ba);
    writeFile("Bfield.txt", r, B);

    int counter = 0;
    while (counter < 100)
    {
        cout << counter << endl;
        calculateMagneticField(m, B, Ba);
        step(m, B);
        char outName[16];
        sprintf(outName, "data/m%06d.ovf", counter);
        writeFile(outName, r, m);
        t += h;
        counter++;
    }

    return 0;
}