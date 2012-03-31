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
    return Particle(Point(*this) + v,step);
}

Particle Particle::operator-(Vector v) {
    return Particle(Point(*this) - v,step);
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

ParticlePolygon GenerativeSphere::generateParticleWhichIntersectsObject(int type) {
    PlaneType *polygon = &object.polygons->at(RAND(object.polygons->size()));
    Vector step;

    do {
        step =  Vector(getRandom() - 0.5,getRandom() - 0.5,getRandom() - 0.5);
        switch(type) {
        case PTYPE_ELECTRON:
            step = step.normalize()*( (*electronVelocityGenerator)() ) - objectStep;
            break;
        case PTYPE_ION:
            step = step.normalize()*( (*ionVelocityGenerator)() ) - objectStep;
            break;
        }
    } while(polygon->getNormal().cos(step) > 0);


    /*if (polygon->getNormal().cos(step) > 0) {
        cout << "fail" << endl;
    }*/

    Point p = GeometryUtils::getRandomPointFromTriangle(*polygon);

    // now we should calculate point on sphere were particle will be initially placed
    // see explanation at page 1 of workbook
    real a = GeometryUtils::getDistanceBetweenPoints(center,p);
    double cos = Vector(center,p).cos(step);
    real step_length = sqrt(a*a*(cos*cos - 1) + radius*radius) - a*cos;

    // now p is point on sphere
    p = p + step.normalize()*step_length;

    /*if (abs(radius- GU::getDistanceBetweenPoints(p,center)) > 0.00001) {
        cout << "fail" << endl;
    }*/

    return ParticlePolygon(Particle(p,step),polygon);
}

void GenerativeSphere::init()
{
    objectStep = object.getStep();
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
