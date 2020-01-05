#ifndef __MATH_DEFS_H__
#define __MATH_DEFS_H__

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <limits>

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989

typedef double scalar;

// Define scalar-valued versions of the Eigen vector types
typedef Eigen::Matrix<scalar, 2, 1> Vector2s;
typedef Eigen::Matrix<scalar, 3, 1> Vector3s;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> VectorXs;

// Define scalar-valued versions of the Eigen matrix types
typedef Eigen::Matrix<scalar, 2, 2> Matrix2s;
typedef Eigen::Matrix<scalar, 3, 3> Matrix3s;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

// // Define scalar-valued
// typedef Eigen::Array<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ArrayXs;

// // Define bool-valued 2D array
// typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ArrayXb;

// Special scalar values (infinity, not-a-number, etc)
// TODO: Actually check that infinity and not-a-number are supported by selected type
#define SCALAR_INFINITY std::numeric_limits<scalar>::infinity()
#define SCALAR_NAN std::numeric_limits<scalar>::signaling_NaN()

template <typename T>
struct Array2D
{
public:
    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> D;
    int x;
    int y;
    Array2D(int x,int y)
    {
        D = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(x,y);

        this->x = x;
        this->y = y;
    }

    const T &operator()(int i, int j) const
    {
        return D(i, j);
    }

    T &operator()(int i, int j)
    {
        return D(i, j);
    }

    void setZero(){
        D.setZero();
    }
    void setOnes(){
        D.setOnes();
    }

    int size() const{
        return D.size();
    }

    const T * data() const {
        return D.data();
    }

    T * data() {
        return D.data();
    }

    static Array2D<T> Zero(int rows, int cols) {
        Array2D<T> D(rows,cols);
        D.D = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(rows,cols);
        return D;
    }

    friend std::ostream &operator<<(std::ostream &os, const Array2D &A)
    {
        os << A.D << std::endl;
        return os;
    }
};


template <typename T>
struct Array3D
{
public:
    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> *data;
    int size;
    int x;
    int y;
    int z;
    inline Array3D(int x,int y, int z)
    {
        data = new Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>[z];
        for(int i = 0; i < z; i++){
            data[i] = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(x,y);
        }
        size = z;
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
        for (int i = 0; i < A.size; i++)
        {
            os << A.data[i] << std::endl
               << std::endl;
        }
        return os;
    }
};


// Define scalar-valued
typedef Array2D<scalar> ArrayXs;

// Define bool-valued 2D array
typedef Array2D<bool> ArrayXb;



#endif
