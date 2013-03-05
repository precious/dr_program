#include <iostream>

typedef double real;

struct Vector3D {
    real x,y,z;
    Vector3D(real _x,real _y,real _z): x(_x), y(_y), z(_z) {};
};

extern "C" int laplace_(int*,Vector3D*,Vector3D*,Vector3D*);

//extern "C" int testfun_(int*);

using namespace std;

int main() {
    Vector3D v1(1.,2.,3.), v2(4.,5.,6.), v3(7.,8.,9.);

    int p = 1;
//    cout << testfun_(&p) << "!" << endl;

    cout << laplace_(&p,&v1,&v2,&v3) << "!" << endl;
    //cout << f << "  " << d << endl;

    return 0;
}
