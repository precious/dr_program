#ifndef GEOMETRY_UTILS_H
#define GEOMETRY_UTILS_H

#include "types.h"
#include <algorithm>
#include <iostream>
#include <cmath>

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED(func) __declspec(deprecated) func
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED(func) func
#endif

using namespace std;

struct ThreePoints;
struct Line;
struct Point;
struct Particle;
struct Plane;
struct Sphere;
struct Vector;
struct Triangle;
struct Particle;
struct Object3D;


namespace Geometry {
    // geometry utils with random
    DEPRECATED(Point getRandomPointFromSphere(Sphere));
    Point getRandomPointFromSphere2(Sphere);
    Point getRandomPointOnSphere(Sphere);
    Vector getRandomOrthogonalVector(Vector);
    Point getRandomPointFromTriangle(ThreePoints&);

    // geometry utils returning distances and lengths
    real getDistanceBetweenPoints(Point,Point);
    real getDistanceBetweenPointAndPlane(ThreePoints&,Point);
    real getChordLength(Sphere,Line);

    // geometry utils for intersections and projections calculation
    Point getPointOnPlaneProjection(ThreePoints&,Point);
    Point getPointOnLineProjection(Line,Point);
    DEPRECATED(Point getPlaneAndLineIntersection(ThreePoints&,Line));
    Point getPlaneAndLineIntersection2(ThreePoints&,Line);
    Point getNearestObject3DAndParticleTrajectoryIntersection(Object3D&,Particle);

    // geometry utils for conditions checking
    DEPRECATED(bool isPointInsideTriangle(ThreePoints&,Point&));
    bool isPointInsideTriangle2(ThreePoints&,Point&);
    bool doesLineIntersectTriangle(ThreePoints&,Line&);
    bool isPointInsideParallelepiped(Point,Point,Point);
    /*static bool isPointInsideObject(Point,Object3D&);*/
    bool doesParticlesTrajectoryIntersectObject(Particle,Object3D&);
    bool doesLineIntersectParallelepiped(Line,Point,Point);
    bool doesLineIntersectSphere(Line,Sphere);

    // other geometry utils
    Point rotatePointAroundLine(Point,Line,double);
}

#endif // GEOMETRY_UTILS_H

