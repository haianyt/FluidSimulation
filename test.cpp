#include <iostream>
#include <Eigen/Core>
#include <Eigen/StdVector>
using namespace std;

template <typename T>
struct Array3D
{
public:
    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> *data;
    int x;
    int y;
    int z;
    inline Array3D(int x,int y, int z)
    {
        data = new Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>[z];
        for(int i = 0; i < z; i++){
            data[i] = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(x,y);
        }
        this->x = x;
        this->y = y;
        this->z = z;
    }

    inline double &operator()(int i, int j, int k)
    {
        return data[k](i, j);
    }

    friend std::ostream &operator<<(std::ostream &os, const Array3D &A)
    {
        for (int i = 0; i < A.z; i++)
        {
            os << A.data[i] << std::endl
               << std::endl;
        }
        return os;
    }

    // inline void operator=(Array3D &other){
    //     x = other.x;
    //     y = other.y;
    //     z = other.z;
    //     data = new Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>[z];
    //     for(int i = 0; i<z; i++){
    //         data[i] = other.data[i];
    //     }
    // }

    inline Array3D<T>& operator=(Array3D<T> &other){
        // std::cout << "hellllll";
        x = other.x;
        y = other.y;
        z = other.z;
        data = new Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>[z];
        for(int i = 0; i<z; i++){
            data[i] = other.data[i];
        }
        return *this;
    }

};

int main(int argc, char const *argv[])
{
    cout << "hello" << endl;
    Array3D<double> A = Array3D<double>(2,3,4);
    // A.data[0] << 1,2,3,4;
    // A.data[1] << 5,6,7,8;
    cout << A;
    // cout << A.x;
    Array3D<double> B = Array3D<double>(2,3,4);

    B = A;

    B(1,1,1) += 20;
    B(1,1,1) += 20;
    cout << B;

    cout << A;
    return 0;
}
