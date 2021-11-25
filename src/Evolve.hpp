#pragma once

#include "VectorField.hpp"
#include "Constants.hpp"
#include "FieldCalculate.hpp"

void step(VectorField *m, const float (&Ba)[3])
{

    // Dorman-Prince (RK45) integration

    // Emergent field array
    VectorField B(N, N, 3);

    // Runge-Kutta slopes
    VectorField k1(N, N, 3);
    VectorField k2(N, N, 3);
    VectorField k3(N, N, 3);
    VectorField k4(N, N, 3);
    VectorField k5(N, N, 3);
    VectorField k6(N, N, 3);

    // Calculate k1
    calculateMagneticField(m, &B, Ba);
    calculateTimeDerivative(m, &B, &k1);
    k1 = k1*h;

    // Calculate k2
    VectorField newMagField = *m + k1*(1/5);
    calculateMagneticField(&newMagField, &B, Ba);
    calculateTimeDerivative(&newMagField, &B, &k2);
    k2 = k2*h;

    // Calculate k3
    newMagField = *m + k1*(3/40) + k2*(9/40);
    calculateMagneticField(&newMagField, &B, Ba);
    calculateTimeDerivative(&newMagField, &B, &k3);
    k3 = k3*h;

    // Calculate k4
    newMagField = *m + k1*(44/45) - k2*(56/15) + k3*(32/9);
    calculateMagneticField(&newMagField, &B, Ba);
    calculateTimeDerivative(&newMagField, &B, &k4);
    k4 = k4*h;

    // Calculate k5
    newMagField = *m + k1*(19372/6561) - k2*(25360/2187) + k3*(64448/6561) - k4*(212/729);
    calculateMagneticField(&newMagField, &B, Ba);
    calculateTimeDerivative(&newMagField, &B, &k5);
    k5 = k5*h;

    // Calculate k6
    newMagField = *m + k1*(9017/3168) - k2*(355/33) - k3*(46732/5247) + k4*(49/176) - k5*(5103/18656);
    calculateMagneticField(&newMagField, &B, Ba);
    calculateTimeDerivative(&newMagField, &B, &k6);
    k6 = k6*h;

    *m = *m + k1*(35/384) + k3*(500/1113) + k4*(125/192) - k5*(2187/6784) + k6*(11/84);

    normalise(m);
}
