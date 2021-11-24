#pragma once

#include <iostream>
#include <vector>

class VectorField
{
    public:
        int Nx;
        int Ny;
        int vecSize;
        std::vector<float> data;

        VectorField(int Nx, int Ny, int vecSize) 
        {
            Nx = Nx;
            Ny = Ny;
            vecSize = vecSize;
            for (int idx=0; idx<Nx*Ny*vecSize; idx++)
                data.push_back(0.);
        }

        ~VectorField() {}

        int getIndex(int i, int j, int k)
        {
            int index = ((i*Ny + j) * vecSize) + k;
            return index;
        }

        VectorField add(VectorField vec)
        {
            VectorField returnVectorField = *this;

            for (int i=0; i<data.size(); i++)
                returnVectorField.data[i] = returnVectorField.data[i] + vec.data[i];

            return returnVectorField;
        }

        VectorField multiplyByScalar(float scalar)
        {
            VectorField returnVectorField = *this;

            for (int i=0; i<data.size(); i++)
                returnVectorField.data[i] = returnVectorField.data[i] * scalar;

            return returnVectorField;
        }
};
