#ifndef TYPES_H
#define TYPES_H

#include <iostream>
#include <string>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <cstring>
#include <assert.h>

#include "types.h"
#include "constants.h"
#include "geometry_utils.h"
#include "time_utils.h"

using namespace std;

struct Point;
struct Plane;
struct Vector;
struct OrientedPlane;
struct Object3D;
struct Particle;

typedef OrientedPlane PlaneType;
typedef pair<Particle,PlaneType*> ParticlePolygon;

extern Point POINT_OF_ORIGIN; // (0,0,0)

// system of coordinates orientation
#define ORIENT_RIGHT_HANDED 1
#define ORIENT_LEFT_HANDED 0

// order - clockwise or counterclockwise
#define ORDER_CW 1
#define ORDER_CCW 0

// Particle types
#define PTYPE_ELECTRON 1
#define PTYPE_ION 0

static bool ORIENTATION = ORIENT_RIGHT_HANDED;
void setOrientation(bool);
bool getOrientation();

class GeometryUtils;

class InvalidPrimitiveException: public runtime_error {
public:
    InvalidPrimitiveException(string *msg): runtime_error(*msg) {}
    InvalidPrimitiveException(string &msg): runtime_error(msg) {}
};



struct Point {
    real x;
    real y;
    real z;

    Point(): x(0),y(0),z(0) {}
    Point(real _x, real _y, real _z): x(_x),y(_y),z(_z) {}
    Point(const Point &p) {x = p.x;y = p.y;z = p.z;}

    Point operator+(Vector v);
    Point operator-(Vector v);
    friend ostream& operator<<(ostream &os, const Point &p) {
        os << '(' << p.x << ',' << p.y << ',' << p.z << ')';
        return os;
    }
    bool operator==(Point a) const {
        return (x == a.x && y == a.y && z == a.z);
    }
    bool operator!=(Point b) const {
        return !(*this == b);
    }
    Point rotateAroundZ(double cos,double sin) {
        return Point(x*cos - y*sin,x*sin + y*cos,z);
    }
    Point rotateAroundY(double cos,double sin) {
        return Point(x*cos + z*sin,y,-x*sin + z*cos);
    }
    real& operator[](int i) {
        switch(i % 3) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
        }
    }
};

struct Vector: public Point {
    Vector(): Point() {}
    Vector(Point p): Point(p) {}
    Vector(real _x, real _y, real _z): Point(_x,_y,_z) {}
    Vector(Point b,Point a): Point(a.x - b.x,a.y - b.y,a.z - b.z) {}

    real operator*(Vector right) {
        return x*right.x + y*right.y + z*right.z;
    }
    Vector operator*(double k) {
        return Vector(k*x,k*y,k*z);
    }
    Vector operator/(double k) {
        return Vector(x/k,y/k,z/k);
    }
    Vector operator+(Vector v) {
        return Vector(x + v.x,y + v.y, z + v.z);
    }
    Vector operator-(Vector v) {
        return Vector(x - v.x,y - v.y, z - v.z);
    }
    real length() {
        return sqrt(x*x + y*y + z*z);
    }
    double cos(Vector right) {
        return ((*this)*right)/(this->length()*right.length());
    }
    Vector vectorProduct(Vector left) {
        return Vector(y*left.z - z*left.y, -x*left.z +
                      z*left.x, x*left.y - y*left.x);
    }
    Vector normalized() {
        double len = length();
       ////// assert(len != 0);
        return Vector(x/len,y/len,z/len);
    }
    Vector resized(real _length) {
        double coef = _length/length();
        if (std::isnan(coef))
            return Vector(x,y,z);
        return Vector(x*coef,y*coef,z*coef);
    }
    void resize(real _length) {
        double coef = _length/length();
        if (!std::isnan(coef)) {
            x *= coef;
            y *= coef;
            z *= coef;
        }
    }
};

template <unsigned int T>
struct Locus { // collection of points
    Point set[T];
    friend ostream& operator<<(ostream &os, const Locus &l) {
        for(unsigned int i = 0;i < T - 1;i++)
            os << l.set[i] << ", ";
        os << l.set[T - 1];
        return os;
    }
};

struct Line: public Locus<2> {
    Vector directionVector;
    Line(Point _a,Point _b): a(set[0]), b(set[1]) {
        set[0] = _a; set[1] = _b;
        directionVector = Vector(_a,_b);
    }
    Line(Point _a,Vector v): a(set[0]), b(set[1]) {
        directionVector = v;
        set[0] = _a;
        set[1] = _a + v;
    }
    Point& a;
    Point& b;
    Point pointByCoef(real coef) {
        return Point(set[0].x + coef*directionVector.x,
                     set[0].y + coef*directionVector.y,
                     set[0].z + coef*directionVector.z);
    }
};

struct ThreePoints : public Locus<3> {
    ThreePoints(): a(set[0]), b(set[1]), c(set[2]) {}
    ThreePoints(const ThreePoints &tP): a(set[0]), b(set[1]), c(set[2]) {
        set[0] = tP.set[0]; set[1] = tP.set[1]; set[2] = tP.set[2];
    }
    ThreePoints(Point _a,Point _b,Point _c): a(set[0]), b(set[1]), c(set[2]) {
        set[0] = _a; set[1] = _b; set[2] = _c;
    }
    ThreePoints& operator=(const ThreePoints& right) {
        if (this != &right) {
            memcpy(set, right.set, 3*sizeof(Point));
        }
        return *this;
    }

    Point& a;
    Point& b;
    Point& c;
    virtual Vector getNormal() {
        return Vector(a,b).vectorProduct(Vector(a,c) );
    }
};

struct Triangle: public ThreePoints {
    Triangle(Point _a,Point _b,Point _c): ThreePoints(_a,_b,_c) {}
    Triangle(ThreePoints &tP): ThreePoints(tP) {}
    Point centerOfMass() {
        return Point( (a.x + b.x + c.x) / 3.0,
                      (a.y + b.y + c.y) / 3.0,
                      (a.z + b.z + c.z) / 3.0);
    }
};

struct Plane: public ThreePoints {
    Plane(): ThreePoints() {}
    Plane(Point _a,Point _b,Point _c): ThreePoints(_a,_b,_c) {}
    Plane(ThreePoints &tP): ThreePoints(tP) {}
    Plane(Point, Vector);
    bool doesPointBelongPlane(Point p) {
        /// TODO fix possible error because of mashine precision
        return Vector(a,p)*getNormal() == 0;
    }
};

struct OrientedPlane: public Plane {
    Vector normal;
    OrientedPlane(): Plane(), normal() {}
    OrientedPlane(Point _a,Point _b,Point _c, bool _pointsOrder = ORDER_CCW):
        Plane(_a,_b,_c), pointsOrder(_pointsOrder) {
        initNormal();
    }
    OrientedPlane(ThreePoints &tP,  bool _pointsOrder = ORDER_CCW):
        Plane(tP), pointsOrder(_pointsOrder) {
        initNormal();
    }
    OrientedPlane(Plane p, Vector v): Plane(p), normal(v) {}
    OrientedPlane(Point p, Vector v): Plane(p,v), normal(v) {}
    OrientedPlane(const OrientedPlane &op): Plane((ThreePoints&)op),
        normal(op.normal) {}
    Vector getNormal() {
        return normal;
    }

private:
    bool pointsOrder;
    void initNormal() {
        Vector ab(a,b);
        Vector ac(a,c);
        // get cross product
        normal = ab.vectorProduct(ac);
        if (pointsOrder == ORDER_CW) {
            normal = normal*(-1);
        }
        ///////// assert (!isnan(normal.normalized().x) && !isnan(normal.normalized().y) && !isnan(normal.normalized().z));
        if (isnan(normal.normalized().x) || isnan(normal.normalized().y) || isnan(normal.normalized().z)) {
            cout << a << "  " << b << "  " << c << endl;
            assert(false);
        }
    }
};

struct Particle: public Point {
public:
    Vector step;
    real ttl;
    Particle operator+(Vector v);
    Particle operator-(Vector v);
    Particle(): Point(), step(), ttl(-1) {}
    Particle(Point p, Vector s,real ttl_ = -1): Point(p), step(s), ttl(ttl_) {}
};

struct Sphere {
    Point center;
    real radius;
    Sphere(): center(), radius(0) {}
    Sphere(Point _p, real _r): center(_p), radius(_r) {}
    Sphere(const Sphere &_s): center(_s.center), radius(_s.radius) {}
};

struct Object3D: public Sphere {
    Vector front;
    Point maxCoords, minCoords;
    Point nearestPoint, furthermostPoint; // relatively to front of the object
    vector<PlaneType> *polygons;
    velocity speed;

    Object3D(int polygonsNumber, Vector _front = Vector(100,0,0)): front(_front) {
        polygons = new vector<PlaneType>(polygonsNumber);
        init();
    }
    Object3D(vector<PlaneType> *_polygons, Vector _front = Vector(100,0,0)):
        front(_front), polygons(_polygons) {
        init();
    }

    void init();

    PlaneType& operator[](int i) {
        return polygons->at(i);
    }

    Vector step() {
        return front.normalized()*speed;
    }
};

struct GenerativeSphere: public Sphere {
private:
    velocity electronVelocityGenerator() {
        static MaxwellDistributionSpeedGenerator generator =
                // getGaussianDistributionGenerator(ELECTRON_VELOCITY,ELECTRON_VELOCITY/4.0);
                getMaxwellDistributionSpeedGenerator(ELECTRON_VELOCITY_M,ELECTRON_VELOCITY_D);
        return generator();
    }
    velocity ionVelocityGenerator() {
        static MaxwellDistributionSpeedGenerator generator =
                // getGaussianDistributionGenerator(ION_VELOCITY,ION_VELOCITY/4.0);
                getMaxwellDistributionSpeedGenerator(ION_VELOCITY_M,ION_VELOCITY_D);
        return generator();
    }
    Object3D &object;
    Vector objectStep;
    // Sphere sphereAroundObject;
    void init();

public:
    GenerativeSphere(Point _p, real _r,Object3D _object):
        Sphere(_p, _r), object(_object), objectStep(object.step()) {}
    GenerativeSphere(const Sphere &_s,Object3D _object):
        Sphere(_s), object(_object), objectStep(object.step()) {}

    Particle generateRandomParticle(int);
    //Particle generateParticleWhichIntersectsObject(int);
    ParticlePolygon generateParticleWhichIntersectsObject(int,bool);
};

#endif // TYPES_H

