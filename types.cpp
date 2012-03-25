#include "types.h"

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
    return Particle(Point(*this) + v,step/*,weight*/);
}

Plane::Plane(Point p, Vector v):
    ThreePoints(p,p + GeometryUtils::getRandomOrthogonalVector(v),
                p + GeometryUtils::getRandomOrthogonalVector(v)) {}

Particle GenerativeSphere::generateRandomParticle(int type) {
    Point initialPosition = GeometryUtils::getRandomPointFromSphere(*this);
    Vector step(getRandom() - 0.5,getRandom() - 0.5,getRandom() - 0.5);
    switch(type) {
    case PTYPE_ELECTRON:
        step = step.normalize()*( (*electronVelocityGenerator)() ) - objectStep;
        break;
    case PTYPE_ION:
        step = step.normalize()*( (*ionVelocityGenerator)() ) - objectStep;
        break;
    }
    return Particle(initialPosition,step);
}

Particle GenerativeSphere::generateParticleWhichIntersectsObject(int type) {
    Point initialPosition = GeometryUtils::getRandomPointFromSphere(*this);
    Vector step(initialPosition,GeometryUtils::getRandomPointFromSphere(sphereAroundObject));
    switch(type) {
    case PTYPE_ELECTRON:
        step = step.normalize()*( (*electronVelocityGenerator)() ) - objectStep;
        break;
    case PTYPE_ION:
        step = step.normalize()*( (*ionVelocityGenerator)() ) - objectStep;
        break;
    }
    return Particle(initialPosition,step);
}

void GenerativeSphere::init()
{
    objectStep = object.getStep();
    sphereAroundObject = Sphere(object.center(),
                                max(GeometryUtils::getDistanceBetweenPoints(object.center(),object.maxCoords),
                                    GeometryUtils::getDistanceBetweenPoints(object.center(),object.minCoords)));
    // expectation and standart deviation are calculated due the 3-sigma rule
    // max and min possible velocities are 2*ELECTRON_VELOCITY and 0 respectively
    electronVelocityGenerator = getGaussianDistributionGenerator(ELECTRON_VELOCITY,ELECTRON_VELOCITY/3.0);
    ionVelocityGenerator = getGaussianDistributionGenerator(ION_VELOCITY,ION_VELOCITY/3.0);
}
