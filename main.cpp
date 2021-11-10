#include <iostream>
#include <fstream>
#include <string>
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

void updateMagnetization(float (&m)[N][N][3], float (&B)[N][N][3])
{
    for (int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            m[i][j][0]++;
        }
    }
}

void calculateMagneticField(float (&m)[N][N][3], float (&B)[N][N][3], const float (&Ba)[3])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            // Reset to zero (adding the applied magnetic field first)
            B[i][j][0] = Ba[0];
            B[i][j][1] = Ba[1];
            B[i][j][2] = Ba[2];

            // Central spin
            B[i][j][0] -= (4 / (DELTA*DELTA)) * m[i][j][0];
            B[i][j][1] -= (4 / (DELTA*DELTA)) * m[i][j][1];
            B[i][j][2] -= (4 / (DELTA*DELTA)) * m[i][j][2];
            B[i][j][0] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j][0];
            B[i][j][1] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j][1];
            B[i][j][2] -= (20 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j][2];

            // Cell directly to left
            B[i][j][0] += m[i-1][j][0] / (DELTA*DELTA);
            B[i][j][1] += m[i-1][j][1] / (DELTA*DELTA);
            B[i][j][2] += m[i-1][j][2] / (DELTA*DELTA);
            B[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j][0];
            B[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j][1];
            B[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j][2];

            // Cell directly to right
            B[i][j][0] += m[i+1][j][0] / (DELTA*DELTA);
            B[i][j][1] += m[i+1][j][1] / (DELTA*DELTA);
            B[i][j][2] += m[i+1][j][2] / (DELTA*DELTA);
            B[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j][0];
            B[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j][1];
            B[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j][2];

            // Cell directly below
            B[i][j][0] += m[i][j-1][0] / (DELTA*DELTA);
            B[i][j][1] += m[i][j-1][1] / (DELTA*DELTA);
            B[i][j][2] += m[i][j-1][2] / (DELTA*DELTA);
            B[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j-1][0];
            B[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j-1][1];
            B[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j-1][2];

            // Cell directly above
            B[i][j][0] += m[i][j+1][0] / (DELTA*DELTA);
            B[i][j][1] += m[i][j+1][1] / (DELTA*DELTA);
            B[i][j][2] += m[i][j+1][2] / (DELTA*DELTA);
            B[i][j][0] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j+1][0];
            B[i][j][1] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j+1][1];
            B[i][j][2] += (8 / (DELTA*DELTA*DELTA*DELTA)) * m[i][j+1][2];

            // Cell to bottom-left
            B[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j-1][0];
            B[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j-1][1];
            B[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j-1][2];

            // Cell to bottom-right
            B[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j-1][0];
            B[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j-1][1];
            B[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j-1][2];

            // Cell to top-left
            B[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j+1][0];
            B[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j+1][1];
            B[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i-1][j+1][2];

            // Cell to top-right
            B[i][j][0] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j+1][0];
            B[i][j][1] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j+1][1];
            B[i][j][2] -= (2 / (DELTA*DELTA*DELTA*DELTA)) * m[i+1][j+1][2];

            // Cell two to left
            B[i][j][0] -= m[i-2][j][0] / (DELTA*DELTA*DELTA*DELTA);
            B[i][j][1] -= m[i-2][j][1] / (DELTA*DELTA*DELTA*DELTA);
            B[i][j][2] -= m[i-2][j][2] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two to right
            B[i][j][0] -= m[i+2][j][0] / (DELTA*DELTA*DELTA*DELTA);
            B[i][j][1] -= m[i+2][j][1] / (DELTA*DELTA*DELTA*DELTA);
            B[i][j][2] -= m[i+2][j][2] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two below
            B[i][j][0] -= m[i][j-2][0] / (DELTA*DELTA*DELTA*DELTA);
            B[i][j][1] -= m[i][j-2][1] / (DELTA*DELTA*DELTA*DELTA);
            B[i][j][2] -= m[i][j-2][2] / (DELTA*DELTA*DELTA*DELTA);

            // Cell two above
            B[i][j][0] -= m[i][j+2][0] / (DELTA*DELTA*DELTA*DELTA);
            B[i][j][1] -= m[i][j+2][1] / (DELTA*DELTA*DELTA*DELTA);
            B[i][j][2] -= m[i][j+2][2] / (DELTA*DELTA*DELTA*DELTA);

        }
    }
}

void writeFile(string name, float (&r)[N][N][2], float (&m)[N][N][3])
{
    ofstream file;
    file.open(name);
    file << "x,y,mx,my,mz" << endl;

    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            file << r[i][j][0] << "," << r[i][j][1] << "," << m[i][j][0] << "," << m[i][j][1] << "," << m[i][j][2] << endl;
        }
    }
}

int main()
{

    float t = 0;  // Time

    float r[N][N][2];  // Declare position array
    float m[N][N][3];  // Declare magnetization array
    float B[N][N][3];  // Declare magnetic field array

    const float Ba[3] = {0., 0., 0.};  // Applied magnetic field

    initialisePosition(r);
    cout << "Setting initial magnetization" << endl;
    initialiseMagnetization(m, r, 1, 1, 3);
    cout << "Calculating effective field" << endl;
    calculateMagneticField(m, B, Ba);
    cout << "Done" << endl;

    // int counter = 0;
    // while (t < T_FIN)
    // {
    //     cout << t / T_FIN << endl;
    //     updateMagnetization(m, B);
    //     updateMagneticField(B);
    //     writeFile("data/" + to_string(counter) + ".txt", r, m);
    //     t += h;
    //     counter++;
    // }

    writeFile("Test", r, m);

    return 0;
}