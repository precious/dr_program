#include "geometry_utils.h"

using namespace std;

template <typename T>
inline char sign(T t) {
    return (t > 0)? 1: (t < 0)? -1: 0;
}

template <typename T>
inline bool inInterval(T x,T a,T b) {
    return x <= max(a,b) && x >= min(a,b);
}

bool Geometry::isPointInsideTriangle(ThreePoints &t,Point k) {
    Vector v0(Point(0,0,0)), v1(k,t.a), v2(k,t.b), v3(k,t.c);
    if (v1 == v0 || v2 == v0 || v3 == v0)
        return true;
    real cos1 = v1.cos(v2), cos2 = v1.cos(v3), cos3 = v2.cos(v3);
    switch( sign(cos1) + sign(cos2) + sign(cos3) ) {
        case -3:
            return true;
        case 3:
        case 1:
            return false;
        default:
            // delta is required to prevent mashine imprecision
            real delta = 0.0001;
            return !(acos(cos1) + acos(cos2) + acos(cos3) < M_PI*2 - delta);
    }
}

bool Geometry::isPointInsideTriangle2(ThreePoints &t,Point k) {
    Point p = getPointOnLineProjection(Line(t.a,t.b),t.c);
    if (Vector(p,t.c).cos(Vector(p,k)) < 0)
        return false;
    p = getPointOnLineProjection(Line(t.c,t.b),t.a);
    if (Vector(p,t.a).cos(Vector(p,k)) < 0)
        return false;
    p = getPointOnLineProjection(Line(t.c,t.a),t.b);
    if (Vector(p,t.b).cos(Vector(p,k)) < 0)
        return false;
    return true;
}

// коэффициент точки пересечения находим подставляя параметрические уравнения прямой
// в векторное уравнение плоскости
Point Geometry::getPlaneAndLineIntersection2(ThreePoints &plane,Line line) {
    real coef = plane.getNormal()*Vector(line.a,plane.a) /
            (plane.getNormal()*line.directionVector);
    return line.pointByCoef(coef);
}

Point Geometry::getPointOnPlaneProjection(ThreePoints& plane,Point p) {
    return getPlaneAndLineIntersection2(plane,Line(p,plane.getNormal()));
}

// this method is deprecated. use getPlaneAndLineIntersection2 instead.
// коэффициент точки пересечения находим подставляя параметрические уравнения прямой
// в каноническое уравнение плоскости
Point Geometry::getPlaneAndLineIntersection(ThreePoints &plane,Line line) {
    real coef1 = (plane.b.y - plane.a.y)*(plane.c.z - plane.a.z) -
            (plane.c.y - plane.a.y)*(plane.b.z - plane.a.z);
    real coef2 = (plane.b.x - plane.a.x)*(plane.c.z - plane.a.z) -
            (plane.c.x - plane.a.x)*(plane.b.z - plane.a.z);
    real coef3 = (plane.b.x - plane.a.x)*(plane.c.y - plane.a.y) -
            (plane.c.x - plane.a.x)*(plane.b.y - plane.a.y);
    real coef = ((line.a.x - plane.a.x)*coef1 -
                 (line.a.y - plane.a.y)*coef2 +
                 (line.a.z - plane.a.z)*coef3) /
            (line.directionVector.x*coef1 - line.directionVector.y*coef2 +
             line.directionVector.z*coef3);

    // inf means that the line is parallel towards plane
    // nan means that the line belongs to plane
    if(isinf(coef) || isnan(coef)) {
        return Point(coef,coef,coef);
    }

    return line.pointByCoef(-coef);
}

// коэффициент точки пересечения находим из условия перпендикулярности напрвляющего вектора
// прямой и вектора, образованного заданной точкой и ее проекией (в качестве координат последней
// берем параметрические уравнения прямой)
Point Geometry::getPointOnLineProjection(Line line,Point point) {
    real coef = line.directionVector*Vector(line.a,point) /
            (line.directionVector*line.directionVector);
    return line.pointByCoef(coef);
}

bool Geometry::doesLineIntersectTriangle(ThreePoints &triangle,Line line) {
    Point intersection = getPlaneAndLineIntersection2(triangle,line);
    if (isinf(intersection.x)) {
        cerr << "INF" << endl;////////////////////////////////////////////////
        /// TODO: here should be checking whether the line intersects
        /// at least one of the triangles side
        /// for this function for finding two lines intersection should be
        /// implemented
        return false;
    }
    if (isnan(intersection.x)) {
        cerr << "NAN" << endl;////////////////////////////////////////////////
        return false;
    }
    bool retVal = isPointInsideTriangle2(triangle,intersection);
    return retVal;
}

Point Geometry::getRandomPointFromSphere(Sphere s) {
    Point randPoint;
    do {
        randPoint.x = (Time::getRandom()*2 - 1)*s.radius;
        randPoint.y = (Time::getRandom()*2 - 1)*s.radius;
        randPoint.z = (Time::getRandom()*2 - 1)*s.radius;
    } while(Geometry::getDistanceBetweenPoints(randPoint,POINT_OF_ORIGIN) > s.radius);
    return s.center + Vector(POINT_OF_ORIGIN,randPoint);
}

Point Geometry::getRandomPointFromSphere2(Sphere s) {
    return s.center + Vector(POINT_OF_ORIGIN,
                             Point(Time::getRandom() - 0.5,Time::getRandom() - 0.5,Time::getRandom() - 0.5))
            * sqrt(s.radius*Time::getRandom(0,s.radius));
}

Point Geometry::getRandomPointOnSphere(Sphere s) {
    return s.center + Vector(s.center,getRandomPointFromSphere2(s)).resized(s.radius);
}

Vector Geometry::getRandomOrthogonalVector(Vector v) {
    Vector a;
    if (v.x != 0) {
        a.y = Time::getRandom();
        a.z = Time::getRandom();
        a.x = (-v.y*a.y - v.z*a.z)/v.x;
    } else if (v.y != 0) {
        a.x = Time::getRandom();
        a.z = Time::getRandom();
        a.y = (-v.x*a.x - v.z*a.z)/v.y;
    } else if (v.z != 0) {
        a.y = Time::getRandom();
        a.x = Time::getRandom();
        a.z = (-v.y*a.y - v.x*a.x)/v.z;
    }
    return a;
}

real Geometry::getDistanceBetweenPoints(Point a,Point b) {
    return sqrt(pow(a.x - b.x,2) + pow(a.y - b.y,2) + pow(a.z - b.z,2));
}

real Geometry::getDistanceBetweenPointAndPlane(ThreePoints& plane,Point p) {
    return getDistanceBetweenPoints(getPointOnPlaneProjection(plane,p),p);
}

real Geometry::getDistanceBetweenPointAndSphere(Sphere& s,Point p) {
    return max<real>(getDistanceBetweenPoints(p,s.center) - s.radius,0);
}


bool Geometry::isPointInsideParallelepiped(Point a,Point v1,Point v2) {
    return a.x <= max(v1.x,v2.x) && a.x >= min(v1.x,v2.x) &&
            a.y <= max(v1.y,v2.y) && a.y >= min(v1.y,v2.y) &&
            a.z <= max(v1.z,v2.z) && a.z >= min(v1.z,v2.z);
}

/*bool GeometryUtils::isPointInsideObject(Point p,Object3D& obj) {
    for (vector<PlaneType>::iterator it = obj.polygons->begin();it != obj.polygons->end();it++) {
        if doesLineIntersectTriangle(*it,)
    }
    return false;
}*/

bool Geometry::doesLineIntersectParallelepiped(Line l,Point p1,Point p2) {
    Plane planePendicularToOX(p1,Vector(1,0,0));
    Plane planePendicularToOY(p1,Vector(0,1,0));
    Plane planePendicularToOZ(p1,Vector(0,0,1));

    Point pX = getPlaneAndLineIntersection2(planePendicularToOX,l);
    Point pY = getPlaneAndLineIntersection2(planePendicularToOY,l);
    Point pZ = getPlaneAndLineIntersection2(planePendicularToOZ,l);

    return  (inInterval(pZ.x,p1.x,p2.x) && inInterval(pZ.y,p1.y,p2.y)) ||
            (inInterval(pX.y,p1.y,p2.y) && inInterval(pX.z,p1.z,p2.z)) ||
            (inInterval(pY.x,p1.x,p2.x) && inInterval(pY.z,p1.z,p2.z));
}

bool Geometry::doesParticlesTrajectoryIntersectObject(Particle p,Object3D &obj) {
    Line line(p,p.step);
    /*if (!doesLineIntersectParallelepiped(line,obj.maxCoords,obj.minCoords))
        return false;*/
    if ( !doesLineIntersectSphere(line,obj) )
        return false;
    for (unsigned int i = 0;i < obj.polygons->size();i++)
        if (doesLineIntersectTriangle(obj.polygons->at(i),line))
            return true;
    return false;
}

// line should intersect sphere
real Geometry::getChordLength(Sphere sphere,Line line) {
    return 2*sqrt(sphere.radius*sphere.radius -
                  pow(getDistanceBetweenPoints(sphere.center,
                                               getPointOnLineProjection(line,sphere.center)),2) );
}

Point Geometry::getRandomPointFromTriangle(ThreePoints& tp) {
    return tp.a + ( Vector(tp.a,tp.b)*Time::getRandom() + Vector(tp.a,tp.c)*Time::getRandom() )*0.5;
}

bool Geometry::doesLineIntersectSphere(Line l,Sphere s) {
    return getDistanceBetweenPoints(getPointOnLineProjection(l,s.center), s.center) <= s.radius;
}

/*Point Geometry::getNearestObject3DAndParticleTrajectoryIntersection(Object3D&,Particle) {
    return Point();
}*/

int Geometry::getIndexOfPolygonThatParicleIntersects(Object3D& obj,Particle p) {
    Line line(p,p.step);
    if ( !doesLineIntersectSphere(line,obj) )
        return -1;
    for (unsigned int i = 0;i < obj.polygons->size();i++)
        if (doesLineIntersectTriangle(obj.polygons->at(i),line))
            return i;
    return -1;
}

// see explanation at pages 9-10 of draft
Point Geometry::rotatePointAroundLine(Point p,Line l,double angle) {
    Point projection = getPointOnLineProjection(l,p);
    Vector j(projection,p);
    double length = j.length();
    Vector i = j.vectorProduct(l.directionVector).resized(length);
    return projection - i*length*sin(angle) + j*length*cos(angle);
}

