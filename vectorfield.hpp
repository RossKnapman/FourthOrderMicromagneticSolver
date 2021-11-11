#include <iostream>

#ifndef VECTOR_FIELD
#define VECTOR_FIELD

template <typename T, int Nx, int Ny, int vecSize>
class VectorField
{
    public:
        T*** data;

        VectorField()
        {
            data = new T**[Nx];
            for (int i=0; i<Nx; i++)
            {
                data[i] = new T*[Ny];
        
                for (int j=0; j<Ny; j++)
                {
                    data[i][j] = new T[vecSize];
                }
            }
        }

        VectorField(VectorField const& other)  // Copy constructor
        {
            data = new T**[Nx];
            for (int i=0; i<Nx; i++)
            {
                data[i] = new T*[Ny];
        
                for (int j=0; j<Ny; j++)
                {
                    data[i][j] = new T[vecSize];
                    
                    for (int component=0; component<vecSize; component++)
                    {
                        data[i][j][component] = other.data[i][j][component];
                    }
                }
            }
        }

        ~VectorField() { }

        VectorField& operator+=(VectorField const& rhs)
        {
            for (int i=0; i<Nx; i++)
            {
                for (int j=0; j<Ny; j++)
                {
                    for (int component=0; component<vecSize; component++)
                    {
                        data[i][j][component] += rhs.data[i][j][component];
                    }
                }
            }
            return *this;
        }

        VectorField& operator+(VectorField const& rhs)
        {
            VectorField<T, Nx, Ny, vecSize> result(*this);
            return result += rhs;
        }

        VectorField& operator*=(float scalar)
        {
            for (int i=0; i<Nx; i++)
            {
                for (int j=0; j<Ny; j++)
                {
                    for (int component=0; component<vecSize; component++)
                    {
                        data[i][j][component] *= scalar;
                    }
                }
            }
            return *this;
        }

        VectorField& operator*(float scalar)
        {
            VectorField<T, Nx, Ny, vecSize> result(*this);
            return result *= scalar;
        }
};

#endif