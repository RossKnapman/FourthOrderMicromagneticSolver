#include <iostream>
#include <array>
#include <cmath>

#define N 128
#define DELTA 0.2
#define h 1e-3
#define T_FIN 0.1

using namespace std;

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

void initialiseMagnetization(float (&m)[N][N][3], float (&r)[N][N][2], float X, float Y, float pol, float charge, float w)
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

void updateMagnetization(float (&m)[N][N][3])
{
    for (int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            m[i][j][0]++;
        }
    }
}

void updateMagneticField(float (&m)[N][N][3])
{
}

void writeFile(string name, float (&r)[N][N][2], float (&m)[N][N][3])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
        }
    }
}

int main()
{

    float t = 0;  // Time

    float r[N][N][2];  // Declare position array
    float m[N][N][3];  // Declare magnetization array
    float B[N][N][3];  // Declare magnetic field array

    initialisePosition(r);
    initialiseMagnetization(m, r, 10, 10, 1, 1, 3);

    // while (t < T_FIN)
    // {
    //     cout << t / T_FIN << endl;
    //     updateMagnetization(m);
    //     // updateMagneticField(B);
    //     t += h;
    // }

    writeFile("Test", r, m);

    return 0;
}