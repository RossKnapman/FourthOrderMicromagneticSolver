#include <iostream>
#include "vectorfield.hpp"

#define N 16

using namespace std;

int main()
{
    VectorField<float, N, N, 3> vecField1;
    vecField1.data[4][5][1]= 5.;
    VectorField<float, N, N, 3> vecField2;
    vecField2.data[4][5][1] = 10.;
    VectorField<float, N, N, 3> vecField3;
    vecField3.data[4][5][1] = 100.;

    VectorField<float, N, N, 3> vecField4;

    // vecField1 += vecField2;
    // cout << vecField1.data[4][5][1];

    vecField4 = vecField1 * 5;
    cout << vecField4.data[4][5][1] << endl;
    // delete &vecField1;

    // cout << result(511, 38, 1) << endl;
}