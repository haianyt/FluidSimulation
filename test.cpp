#include <iostream>
#include <Eigen/Core>
#include <Eigen/StdVector>
using namespace std;

typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ArrayXs;

struct Array3D
{
public:
    ArrayXs *data;
    int size;
    inline Array3D(int N)
    {
        data = new ArrayXs[N];
        for(int i = 0; i < N; i++){
            data[i] = ArrayXs(N,N);
        }
        size = N;
    }

    inline double &operator()(int i, int j, int k)
    {
        return data[i](j, k);
    }

    friend ostream &operator<<(ostream &os, const Array3D &A)
    {
        for (int i = 0; i < A.size; i++)
        {
            os << A.data[i] << endl
               << endl;
        }
        return os;
    }
};

int main(int argc, char const *argv[])
{
    cout << "hello" << endl;
    Array3D A = Array3D(2);
    // A.data[0] << 1,2,3,4;
    // A.data[1] << 5,6,7,8;
    cout << A;
    return 0;
}
