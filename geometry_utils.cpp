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

bool GeometryUtils::isPointInsideTriangle(ThreePoints &t,Point &k) {
    Vector v0(Point(0,0,0)), v1(k,t.a()), v2(k,t.b()), v3(k,t.c());
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

Point GeometryUtils::getPlaneAndLineIntersection2(ThreePoints &plane,Line line) {
    Vector planeNormal = plane.getNormal();
    Vector lineDirection = line.directionVector;
    real coef = planeNormal.x*(line.a().x - plane.a().x) +
                planeNormal.y*(line.a().y - plane.a().y) +
                planeNormal.z*(line.a().z - plane.a().z);
    coef /= planeNormal*lineDirection;
    return line.getPointByCoef(coef);
}

Point GeometryUtils::getPointOnPlaneProjection(ThreePoints& plane,Point p) {
    return getPlaneAndLineIntersection2(plane,Line(p,plane.getNormal()));
}

Point GeometryUtils::getPlaneAndLineIntersection(ThreePoints &plane,Line line) {
    real coef1 = (plane.b().y - plane.a().y)*(plane.c().z - plane.a().z) -
            (plane.c().y - plane.a().y)*(plane.b().z - plane.a().z);
    real coef2 = (plane.b().x - plane.a().x)*(plane.c().z - plane.a().z) -
            (plane.c().x - plane.a().x)*(plane.b().z - plane.a().z);
    real coef3 = (plane.b().x - plane.a().x)*(plane.c().y - plane.a().y) -
            (plane.c().x - plane.a().x)*(plane.b().y - plane.a().y);
    real coef = ((line.a().x - plane.a().x)*coef1 -
                 (line.a().y - plane.a().y)*coef2 +
                 (line.a().z - plane.a().z)*coef3) /
                (line.directionVector.x*coef1 - line.directionVector.y*coef2 +
                 line.directionVector.z*coef3);

    // inf means that the line is parallel towards plane
    // nan means that the line belongs to plane
    if(isinf(coef) || isnan(coef)) {
        return Point(coef,coef,coef);
    }

    return line.getPointByCoef(coef);
}

bool GeometryUtils::doesLineIntersectTriangle(ThreePoints &triangle,Line &line) {
    Point intersection = getPlaneAndLineIntersection2(triangle,line);
    if (isinf(intersection.x)) {
        cout << "INF" << endl;////////////////////////////////////////////////
        /// TODO: here should be checking whether the line intersects
        /// at least one of the triangles side
        /// for this function for finding two lines intersection should be
        /// implemented
        return false;
    }
    if (isnan(intersection.x)) {
        cout << "NAN" << endl;////////////////////////////////////////////////
        return false;
    }
    bool retVal = isPointInsideTriangle(triangle,intersection);
    return retVal;
}

Point GeometryUtils::getRandomPointFromSphere(Sphere s) {
    Point randPoint;
    do {
        randPoint.x = (getRandom()*2 - 1)*s.radius;
        randPoint.y = (getRandom()*2 - 1)*s.radius;
        randPoint.z = (getRandom()*2 - 1)*s.radius;
    } while(Vector(randPoint,s.center).length() > s.radius);
    return randPoint;
}

Vector GeometryUtils::getRandomOrthogonalVector(Vector v) {
    Vector a;
    if (v.x != 0) {
        a.y = getRandom();
        a.z = getRandom();
        a.x = (-v.y*a.y - v.z*a.z)/v.x;
    } else if (v.y != 0) {
        a.x = getRandom();
        a.z = getRandom();
        a.y = (-v.x*a.x - v.z*a.z)/v.y;
    } else if (v.z != 0) {
        a.y = getRandom();
        a.x = getRandom();
        a.z = (-v.y*a.y - v.x*a.x)/v.z;
    }
    return a;
}

real GeometryUtils::getDistanceBetweenPoints(Point a,Point b) {
    return sqrt(pow(a.x - b.x,2) + pow(a.y - b.y,2) + pow(a.z - b.z,2));
}

real GeometryUtils::getDistanceBetweenPointAndPlane(ThreePoints& plane,Point p) {
    return getDistanceBetweenPoints(getPointOnPlaneProjection(plane,p),p);
}

Particle GeometryUtils::getFastestParticle(vector<Particle>& vec,Point target,function<bool (Particle)> check) {
    Particle fastestParticle = Particle(Point(),Vector(),0);
    real longestProjectionLength = 0, projectionLength, cos;
    for (vector<Particle>::iterator it = vec.begin();it != vec.end();it++) {
        cos = Vector((*it),target).cos((*it).step);
        if (cos <= 0)
            continue;
        projectionLength = (*it).step.length()*cos;
        if (projectionLength > longestProjectionLength && check(*it)) {
            cout << "i'm here" << endl;//////////////////////////////////////////////////////
            longestProjectionLength = projectionLength;
            fastestParticle = *it;
        }
    }
    return fastestParticle;
}

bool GeometryUtils::isPointInsideParallelepiped(Point a,Point v1,Point v2) {
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

bool GeometryUtils::doesLineIntersectParallelepiped(Line l,Point p1,Point p2) {
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

bool GeometryUtils::doesParticlesTrajectoryIntersectObject(Particle p,Object3D &obj) {
    Line line(p,p.step);
    if (!doesLineIntersectParallelepiped(line,obj.maxCoords,obj.minCoords))
        return false;
    for (vector<PlaneType>::iterator it = obj.polygons->begin();it != obj.polygons->end();it++)
        if (doesLineIntersectTriangle(*it,line))
            return true;
    return false;
}

