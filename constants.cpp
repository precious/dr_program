#include "constants.h"

// http://ceres.hsc.edu/homepages/classes/astronomy/spring99/Mathematics/sec10.html

map<string,velocity> *initVelocityMap() {
    map<string,velocity> *vm = new map<string,velocity>();
    (*vm)[string("orbitalVelocity")] = 7907.343098064;
    (*vm)[string("ionVelocity")] = 9090.364708796;
    (*vm)[string("electronVelocity")] = 389206.742021674;
    (*vm)[string("maxElectronVelocity")] = vm->at(string("electronVelocity"));
    (*vm)[string("minElectronVelocity")] = vm->at(string("electronVelocity"))*0.5;
    (*vm)[string("maxIonVelocity")] = vm->at(string("ionVelocity"));
    (*vm)[string("minIonVelocity")] = vm->at(string("ionVelocity"))*0.5;
    return vm;
}

map<string,velocity> *velocityMap = initVelocityMap();

velocity getVelocitySi(char* velocityName) {
    return velocityMap->at(string(velocityName));
}

velocity getVelocityGl(char* velocityName) {
    return velocityMap->at(string(velocityName))*1000;
}
