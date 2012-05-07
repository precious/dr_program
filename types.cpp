#include "types.h"
#include <assert.h>

Point POINT_OF_ORIGIN = Point(0,0,0);

Point Point::operator+(Vector v) {
    return Point(x + v.x,y + v.y, z + v.z);
}

Point Point::operator-(Vector v) {
    return Point(x - v.x,y - v.y, z - v.z);
}

Particle Particle::operator+(Vector v) {
    return Particle(Point(*this) + v,step,ttl,type,polygonIndex);
}

Particle Particle::operator-(Vector v) {
    return Particle(Point(*this) - v,step,ttl,type,polygonIndex);
}

Plane::Plane(Point p, Vector v):
    ThreePoints(p,p + Geometry::getRandomOrthogonalVector(v),
                p + Geometry::getRandomOrthogonalVector(v)) {}

void GenerativeSphere::checkForIntersectionsAndSetTtl(Particle &p) {
    p.polygonIndex = Geometry::getIndexOfPolygonThatParicleIntersects(object,p);
    if (p.polygonIndex != -1) {
        p.ttl = Geometry::getDistanceBetweenPoints(p,
            Geometry::getPlaneAndLineIntersection2(object.polygons->at(p.polygonIndex),Line(p,p.step)))/p.step.length();
    } else { // see description at page 11 of the draft
        p.ttl = 2*radius*p.step.cos(Vector(p,center)) / p.step.length();
    }
}

Particle GenerativeSphere::generateParticleInSphere(int type) {
    Point initialPosition = Geometry::getRandomPointFromSphere2(*this);
    Vector step(Time::getRandom() - 0.5,Time::getRandom() - 0.5,Time::getRandom() - 0.5);
    switch(type) {
    case PTYPE_ELECTRON:
        step = step.resized(electronVelocityGenerator()) - objectStep;
        break;
    case PTYPE_ION:
        step = step.resized(ionVelocityGenerator()) - objectStep;
        break;
    }

    Particle p(initialPosition,step,0,type);
    checkForIntersectionsAndSetTtl(p);
    return p;
}

Particle GenerativeSphere::generateParticleWhichIntersectsObject(int type,bool isParticleOnSphere) {
    int polygonIndex = RAND(object.polygons->size());
    Vector n = object.polygons->at(polygonIndex).getNormal().normalized();
    velocity particleSpeed;

    switch(type) {
    case PTYPE_ELECTRON:
        particleSpeed = electronVelocityGenerator();
        break;
    case PTYPE_ION:
        particleSpeed = ionVelocityGenerator();
        break;
    }

    Point p = Geometry::getRandomPointFromTriangle(object.polygons->at(polygonIndex));
    Line auxLine(p,(object.polygons->at(polygonIndex).a != p)?
                     object.polygons->at(polygonIndex).a: object.polygons->at(polygonIndex).b);

    // see explanation at page 2 of draft
    if (particleSpeed <= -object.speed*n.cos(objectStep)) {
        particleSpeed = Time::getRandom(max<velocity>(0.1,-object.speed*n.cos(objectStep)),2.*ELECTRON_VELOCITY_M); /// TODO fix me
    }

    double cos = Time::getRandom(-1,min<velocity>(1,object.speed*n.cos(objectStep)/particleSpeed));
    // see explanation at page 8 of draft
    Vector s(p,Geometry::rotatePointAroundLine(p + n,auxLine,acos(cos)));
    s = s.resized(particleSpeed) - objectStep;

    real distanceBetweenParticleAndPolygon;
    if (isParticleOnSphere) {
        // now we should calculate point on sphere were particle will be initially placed
        // see explanation at page 1 of draft
        real a = Geometry::getDistanceBetweenPoints(center,p);
        cos = Vector(center,p).cos(s);
        distanceBetweenParticleAndPolygon = sqrt(a*a*(cos*cos - 1) + radius*radius) + a*cos;

        // now p is point on sphere
        p = p - s.resized(distanceBetweenParticleAndPolygon);

        //assert (abs(radius - GU::getDistanceBetweenPoints(p,center)) <= 0.00001);
    } else {
        // see explanation at page 5 of draft
        real distanceBetweenParticleAndPolygon = sqrt(Time::getRandom(object.radius,radius)*radius);
        p = p - s.resized(distanceBetweenParticleAndPolygon);
    }

    return Particle(p,s,distanceBetweenParticleAndPolygon/particleSpeed,type,polygonIndex);
}

Particle GenerativeSphere::generateParticleOnSphere(int type) {
    Point initialPosition = Geometry::getRandomPointOnSphere(*this);
    Vector step(initialPosition,Geometry::getRandomPointOnSphere(*this));
    switch(type) {
    case PTYPE_ELECTRON:
        step = step.resized(electronVelocityGenerator()) - objectStep;
        break;
    case PTYPE_ION:
        step = step.resized(ionVelocityGenerator()) - objectStep;
        break;
    }

    Particle p(initialPosition,step,0,type);
    checkForIntersectionsAndSetTtl(p);
    return p;
}

void GenerativeSphere::populateArray(Particle *particles,int number,int type,int FLAGS) {
    int n = 0;
    bool isOnSphere = FLAGS & GEN_ON_SPHERE;
    if (FLAGS & GEN_INTERSECT_OBJ)
        for(;n < number;++n) {
            particles[n] = generateParticleWhichIntersectsObject(type,isOnSphere);
        }
    else if (FLAGS & GEN_RANDOM)
        for(;n < number;++n) {
            particles[n] = generateParticleInSphere(type);
        }
    else if (isOnSphere)
        for(;n < number;++n) {
            particles[n] = generateParticleOnSphere(type);
        }
}

void Object3D::init() {
    charge = 0;
    speed = ORBITAL_VELOCITY;
    nearestPoint = furthermostPoint = maxCoords = minCoords = polygons->at(0).set[0];
    for(vector<PlaneType>::iterator it = polygons->begin();it != polygons->end();it++)
        for(int i = 0;i < 3;i++) {
            if (Vector(nearestPoint,(*it).set[i]).cos(front) > 0)
                nearestPoint = (*it).set[i];
            if (Vector(furthermostPoint,(*it).set[i]).cos(front) < 0)
                furthermostPoint = (*it).set[i];

            if (maxCoords.x < (*it).set[i].x) maxCoords.x = (*it).set[i].x;
            if (maxCoords.y < (*it).set[i].y) maxCoords.y = (*it).set[i].y;
            if (maxCoords.z < (*it).set[i].z) maxCoords.z = (*it).set[i].z;
            if (minCoords.x > (*it).set[i].x) minCoords.x = (*it).set[i].x;
            if (minCoords.y > (*it).set[i].y) minCoords.y = (*it).set[i].y;
            if (minCoords.z > (*it).set[i].z) minCoords.z = (*it).set[i].z;
        }
    center = Point((maxCoords.x + minCoords.x)/2,
                   (maxCoords.y + minCoords.y)/2,
                   (maxCoords.z + minCoords.z)/2);
    radius = Geometry::getDistanceBetweenPoints(center,maxCoords);
    polygonsCharges = new double[polygons->size()];
    for(unsigned int i = 0;i < polygons->size();++i) {
        polygonsCharges[i] = 0;
    }
}
