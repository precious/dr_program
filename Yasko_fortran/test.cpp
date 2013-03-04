#include <iostream>

struct Vector3D {
    double x,y,z;
    Vector3D(double _x,double _y,double _z): x(_x), y(_y), z(_z) {};
};

extern "C" int laplace_(int,Vector3D*,Vector3D*,Vector3D*);

extern "C" int testfun_(int);

using namespace std;

int main() {
    Vector3D v1(1.,0.,0.), v2(0.,1.,0.), v3(0.,0.,1.);
    
    cout << testfun_(15) << "!" << endl;
    
    //cout << laplace_(10,&v1,&v2,&v3) << "!" << endl;
    //cout << f << "  " << d << endl;
        
    return 0;
}
