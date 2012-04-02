#include "types.h"
#include <assert.h>

Point POINT_OF_ORIGIN = Point(0,0,0);

bool getOrientation() {
    return ORIENTATION;
}

void setOrientation(bool orientation) {
    ORIENTATION = orientation;
}


Point Point::operator+(Vector v) {
    return Point(x + v.x,y + v.y, z + v.z);
}

Point Point::operator-(Vector v) {
    return Point(x - v.x,y - v.y, z - v.z);
}

Particle Particle::operator+(Vector v) {
    return Particle(Point(*this) + v,step,ttl);
}

Particle Particle::operator-(Vector v) {
    return Particle(Point(*this) - v,step,ttl);
}

Plane::Plane(Point p, Vector v):
    ThreePoints(p,p + GeometryUtils::getRandomOrthogonalVector(v),
                p + GeometryUtils::getRandomOrthogonalVector(v)) {}

Particle GenerativeSphere::generateRandomParticle(int type) {
    Point initialPosition = GeometryUtils::getRandomPointFromSphere2(*this);
    Vector step(getRandom() - 0.5,getRandom() - 0.5,getRandom() - 0.5);
    switch(type) {
    case PTYPE_ELECTRON:
        step = step.resized((*electronVelocityGenerator)()) - objectStep;
        break;
    case PTYPE_ION:
        step = step.resized((*ionVelocityGenerator)()) - objectStep;
        break;
    }
    return Particle(initialPosition,step);
}

ParticlePolygon GenerativeSphere::generateParticleWhichIntersectsObject(int type,bool isParticleOnSphere = false) {
    PlaneType *polygon = &object.polygons->at(RAND(object.polygons->size()));
    Vector n = polygon->getNormal().normalized();
    velocity particleSpeed;

    switch(type) {
    case PTYPE_ELECTRON:
        particleSpeed = ( (*electronVelocityGenerator)() );
        break;
    case PTYPE_ION:
        particleSpeed = ( (*ionVelocityGenerator)() );
        break;
    }

    Point p = GeometryUtils::getRandomPointFromTriangle(*polygon);
    Line auxLine(p,(polygon->a != p)? polygon->a: polygon->b);

    // see explanation at page 2 of draft
    if (particleSpeed < object.speed*n.cos(objectStep)) {
        particleSpeed = getRandom(object.speed*n.cos(objectStep),2.*ELECTRON_VELOCITY); /// TODO fix me
    }
    double cos = getRandom(-1,object.speed*n.cos(objectStep)/particleSpeed);
    // see explanation at page 8 of draft
    Vector s(p,GU::rotatePointAroundLine(p + n,auxLine,acos(cos)));
    s = s.resized(particleSpeed) - objectStep;
/*
    if (n.cos(s) >= 0) {
        cout << n.cos(s) << endl;/////////////////////
    }
*/
    real ttl = 0;
    if (isParticleOnSphere) {
        // now we should calculate point on sphere were particle will be initially placed
        // see explanation at page 1 of draft
        real a = GeometryUtils::getDistanceBetweenPoints(center,p);
        cos = Vector(center,p).cos(s);
        real step_length = sqrt(a*a*(cos*cos - 1) + radius*radius) + a*cos;

        // now p is point on sphere
        p = p - s.resized(step_length);

        /*if (abs(radius- GU::getDistanceBetweenPoints(p,center)) > 0.00001) {
            cout << "fail" << endl;
        }*/////////////////////////////////
    } else {
        // see explanation at page 5 of draft
        real distanceBetweenParticleAndPolygon = sqrt(getRandom(object.radius,radius)*radius);
        p = p - s.resized(distanceBetweenParticleAndPolygon);
        ttl = distanceBetweenParticleAndPolygon/particleSpeed;
        //Point p2 = GU::getRandomPointFromSphere2(*this);
        //p = p - s.resized(GU::getDistanceBetweenPoints(p,p2));
    }

    return ParticlePolygon(Particle(p,s,ttl),polygon);
}

void GenerativeSphere::init()
{
    objectStep = object.step();
    /*sphereAroundObject = Sphere(object.center(),
                                GeometryUtils::getDistanceBetweenPoints(object.center(),object.maxCoords));
    */
    // expectation and standart deviation are calculated due the 4-sigma rule
    // max and min possible velocities are 2*ELECTRON_VELOCITY and 0 respectively
    electronVelocityGenerator = getGaussianDistributionGenerator(ELECTRON_VELOCITY,ELECTRON_VELOCITY/4.0);
    ionVelocityGenerator = getGaussianDistributionGenerator(ION_VELOCITY,ION_VELOCITY/4.0);
}

void Object3D::init() {
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
    radius = GeometryUtils::getDistanceBetweenPoints(center,maxCoords);
}
