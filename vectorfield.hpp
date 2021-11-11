#ifndef VECTOR_FIELD
#define VECTOR_FIELD

template <typename T, int Nx, int Ny, int vecSize>
class VectorField
{
    private:
        T*** data = new T**[Nx];

    public:
        T& operator() (int i, int j, int k) { return data[i][j][k]; }
        VectorField()
        {
            for (int i=0; i<Nx; i++)
            {
                data[i] = new T*[Ny];
        
                for (int j=0; j<Ny; j++)
                {
                    data[i][j] = new T[vecSize];
                }
            }
        }
        ~VectorField() { delete [] data; }
};

#endif