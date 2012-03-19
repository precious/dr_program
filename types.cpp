#include "types.h"

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

Point GenerativeSurface::generatePoint() {
    Point p = center + i*(getRandom()*width)+ j*(getRandom()*height);
    ////////////cout << "w: " << width << endl;//////////////////////////////////////////////
    ////////////cout << "h: " << height << endl;//////////////////////////////////////////////
    return p;
}

Particle GenerativeSurface::generateParticle(int type) {
    Point initialPosition = generatePoint();
    Vector randDeviation(getRandom() - 0.5,getRandom() - 0.5,getRandom() - 0.5);
    Vector step = objectDirection;
    switch(type) {
    case PTYPE_ELECTRON:
        step = step + randDeviation.normalize()*( (*electronVelocityGenerator)() );
        break;
    case PTYPE_ION:
        step = step + randDeviation.normalize()*( (*ionVelocityGenerator)() );
        break;
    }
    /// TODO weight should be properly specified
    return Particle(initialPosition,step,1);
}

Particle Particle::operator+(Vector v) {
    return Particle(Point(*this) + v,step,weight);
}

Plane::Plane(Point p, Vector v):
    ThreePoints(p,p + GeometryUtils::getRandomOrthogonalVector(v),
                p + GeometryUtils::getRandomOrthogonalVector(v)) {}

