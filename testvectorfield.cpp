#include <iostream>
#include "vectorfield.hpp"

#define N 16

using namespace std;

int main()
{
    VectorField<float, N, N, 3> vecField1;
    vecField1(4, 5, 1) = 5.;
    VectorField<float, N, N, 3> vecField2;
    vecField2(4, 5, 1) = 10.;
    VectorField<float, N, N, 3> vecField3;
    vecField3(4, 5, 1) = 100.;

    VectorField<float, N, N, 3> result = vecField1 + vecField2 + vecField3;

    // cout << result(511, 38, 1) << endl;
}