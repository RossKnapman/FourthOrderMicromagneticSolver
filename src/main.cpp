#include <iostream>
#include <string>
#include <array>
#include <cmath>

#include "VectorField.hpp"
#include "Maths.hpp"
#include "Constants.hpp"
#include "Initialise.hpp"
#include "IO.hpp"
#include "FieldCalculate.hpp"
#include "Evolve.hpp"

using namespace std;

float t = 0;  // Total simulation time

int main()
{
    VectorField m(N, N, 3);  // Declare magnetization array
    VectorField B(N, N, 3);  // Declare emergent field array
    VectorField pos(N, N, 2);  // Declare position array

    const float Ba[3] = {0., 0., 1.};  // Applied magnetic field

    initialisePosition(&pos);
    initialiseMagnetizationSkyrmion(&m, &pos, -1, 1, 3);

    calculateMagneticField(&m, &B, Ba);

    VectorField dmdt(N, N, 3);
    calculateTimeDerivative(&m, &B, &dmdt);

    int counter = 0;
    while (counter < STEPS)
    {
        cout << counter << endl;
        calculateMagneticField(&m, &B, Ba);
        step(&m, &B);
        char outName[16];
        sprintf(outName, "data/m%06d.ovf", counter);
        writeFile(outName, &m, t);
        t += h;
        counter++;
    }

    return 0;
}
